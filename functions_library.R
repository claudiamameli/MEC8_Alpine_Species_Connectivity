# Library -----------------------------------------------------------------
library(tidyverse)
library(sp)
library(colorRamps)
library(RColorBrewer)
library(paletteer)
library(terra)
library(sf)
library(rworldmap)
library(gpkg)
library(raster)
library(igraph)
library(ggplot2)
library(ggridges)
library(readxl)
library(patchwork)
library(networktools)
library(cluster)
library(cowplot)
library(caseconverter) 
library(mgcv)
library(gam)




# Functions ---------------------------------------------------------------
# Function for plotting nodes
plot_nodes <- function(df, scenario, buffer) {
  ggplot(df, aes(x = longitude, y = latitude)) +
    geom_point(size = 0.4) +
    coord_fixed() +
    theme_bw() +
    labs(title = paste(scenario, ", Buffer", buffer)) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none"
    )
}

# Function for edge-to-edge network construction based on Euclidean disctance and analysis
euclidean_network_e2e <- function(d, df_patch) {
  # matrix of distances between all pairs of patches
  mat_dist <- as.matrix(dist(dplyr::select(df_patch, longitude,latitude), method="euclidean", diag=TRUE, upper=TRUE))
  print('1/7: Preliminary matrix created')
  
  # calculate and extrapolate the radius of individual patches
  df_patch$radius <- sqrt(df_patch$area_m2 / pi)
  radius_list <- dplyr::select(df_patch, radius)
  
  # creation  of a new distance matrix considering edge to edge distance
  e2e_mat_dist <- mat_dist
  e2e_mat_dist[] <- 0 
  
  # this avoids loops, which was taking too long. 
  radius_matrix <- outer(radius_list$radius, radius_list$radius, `+`)
  e2e_mat_dist <- mat_dist - radius_matrix
  print('2/7: Edge to edge matrix created')
  
  # This is to make sure diagonal and negative values are set to 0 
  diag(e2e_mat_dist) <- 0 
  e2e_mat_dist[e2e_mat_dist < 0] <- 0
  print('3/7: Matrix cleaned up')
  
  # adjacency matrix -> 1 = link between patch i and j; 0 = no link between patch i and j
  # assign 1 if distance between patch i and j is less or equal to threshold distance
  # assign 0 if distance between patch i and j is greater than threshold distance
  m_adj <- e2e_mat_dist
  m_adj[e2e_mat_dist <= d] <- 1
  m_adj[e2e_mat_dist > d] <- 0
  diag(m_adj) <- 0 # set diagonal to 0
  print('4/7: Binary adjacency matrix created')

  # Create dataframe with connections between patches
  edge_index <- which(m_adj==1, arr.ind=TRUE)
  patch_names <- df_patch$patch
  edges <- data.frame(from = patch_names[edge_index[, "row"]],
                      to = patch_names[edge_index[, "col"]]) %>%
    rowwise() %>%
    mutate(a = min(from, to),
           b = max(from, to)) %>%
    ungroup() %>%
    dplyr::select(from = a, to = b) %>%
    distinct()

  # Create igraph object
  g <- graph_from_data_frame(d=edges, vertices=df_patch, directed=FALSE)
  
  print('5/7: Edges and graph created')
  
  # Connectance (number of realised links / number of all possible links)
  n = nrow(m_adj) # number of patches
  c = sum(m_adj) / (n*(n-1))
  
  # Modularity
  modules <- cluster_infomap(g) # switch to different function if needed
  m <- modularity(modules)
  
  # Centrality metrics
  c_betwenness <- betweenness(g, directed=FALSE)
  c_closeness <- closeness(g)
  deg_list <- degree(g)
  
  # Compute components
  comp <- components(g, mode = "weak")  # Use mode = "strong" for directed Strongly Connected Component (SCC)
  largest_size <- max(comp$csize)  
  total_nodes <- vcount(g)   
  fraction <- largest_size / total_nodes  # Fraction of nodes in largest component
  
  # Calculate percentage of nodes connected
  connected_nodes <- sum(deg_list >= 1)
  percent_connected <- (connected_nodes / nrow(m_adj)) * 100
  
  print('6/7: All metrics computed')
  print('7/7: Returning full list of results')
  
  
  return(list(
    edge_index = edge_index,
    distance = d, 
    connectance = c, 
    modularity = m,
    betweenness = c_betwenness,
    closeness = c_closeness,
    degree = deg_list,
    comp_number = comp$no,
    largest_comp = largest_size,
    fraction_comp = fraction,
    percent_connected = percent_connected,
    graph = g
  ))
}

#' Function for plotting graph 
#' g --> graph from the euclidean_network_e2e function
#' d --> distance used for the network creation
graph_plot_fun <- function(g, d, title){ 
  nodes <- igraph::as_data_frame(g, what = "vertices")

  layout_coords <- nodes %>%
    dplyr::select(name, longitude, latitude) %>%
    filter(name %in% V(g)$name) %>%
    arrange(match(name, V(g)$name)) %>%
    dplyr::select(longitude, latitude) %>% 
    as.matrix()
  plot(g,
       layout = layout_coords,
       vertex.size = log10(nodes$area_m2),
       #vertex.label = V(g)$name,
       vertex.label = "",
       vertex.color = "black",
       edge.color = "red",
       edge.width = 2,
       main = paste0(title, " Spatial Network, d = ",d))
}

# Functions to create a full dataframe with the results from the euclidean_network_e2e function
net_metrics_fun <-  function(full_results, scenario_){
  data.frame(
  distance = full_results$distance,
  percent_connected = full_results$percent_connected,
  connectance = full_results$connectance,
  modularity = full_results$modularity,
  comp_number = full_results$comp_number, 
  comp_largest = full_results$largest_comp,
  comp_fraction = full_results$fraction_comp,
  scenario = scenario_)
}

node_metrics_fun <- function(full_results, original_node_df){
  df <- data.frame(
    patch = names(full_results$betweenness),
    betweenness = full_results$betweenness,
    closeness = full_results$closeness,
    degree = full_results$degree
  )
  
  # Change NaN --> 0 
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
    col[is.nan(col)] <- 0
    return(col)
  })
  
  # Adding node areas and ranking nodes by connectivity measures
  df <- df %>%
    mutate(area_m2 = original_node_df$area_m2) %>% 
    mutate(
      betweenness_rank = ifelse( betweenness != 0, rank(-betweenness, ties.method = "min"), NA),
      closeness_rank = ifelse(closeness != 0 & closeness != 1, rank(-closeness, ties.method = "min"), NA),
      degree_rank = ifelse(degree != 0, rank(-degree, ties.method = "min"), NA),
      area_rank = rank(-area_m2)
    )
  return(df)
}

# Modularity analysis function 
random_mod_cal_fun <- function(g, num_graphs) {
  random_modularity <- numeric(num_graphs)
  for (i in 1:num_graphs) {
    g_rand <- sample_degseq(degree(g), method = "configuration")
    mod_rand <- cluster_infomap(g_rand) #can be changed to other clustering functions - for comparison better to match with same as in euclidean_network_e2e function
    random_modularity[i] <- modularity(mod_rand)
  }
  return(random_modularity)
}
