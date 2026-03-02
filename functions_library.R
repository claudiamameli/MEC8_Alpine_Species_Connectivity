# Library -----------------------------------------------------------------
library(caseconverter) 
library(colorRamps)
library(cowplot)
library(cluster)
library(gam)
library(ggplot2)
library(ggridges)
library(gpkg)
library(igraph)
library(networktools)
library(mgcv)
library(paletteer)
library(patchwork)
library(PresenceAbsence)
library(raster)
library(RColorBrewer)
library(readxl)
library(rworldmap)
library(sf)
library(sp)
library(terra)
library(tidyverse)




# Functions ---------------------------------------------------------------
# Nodes creation function
node_df_fun <- function(original_map, direction, cell_res, buffer_val = NULL){
  # Buffer addition
  if(is.null(buffer_val) == F){
    cat("Adding buffer \n")
    r_buffered <- as.numeric(buffer(original_map, width = buffer_val))
    r_buffered[r_buffered[] == 0] <- NA
    
    cat("Creating patches\n")
    # Connect neighbouring cells (direction = 8 for corners, = 4 for sides only)
    patch_data <- patches(r_buffered, directions = direction)
    
  } else {
    cat("Creating patches\n")
    patch_data <- patches(original_map, directions = direction, zeroAsNA=T)
  }
  
  # Create polygons and find centroids
  cat("Creating polygons and centroids \n")
  patch_polygons <- as.polygons(patch_data, dissolve = TRUE)
  df_centroids <- crds(centroids(patch_polygons), df = TRUE)
  
  # Calculate cell counts and area
  cat("Counting cells \n")
  patch_cell_count <- freq(patch_data, bylayer = FALSE)
  area_m2 <- patch_cell_count$count * cell_res
  
  
  cat("Returning data frame \n")
  return(data.frame(
    patch = patch_polygons$patches,
    area_m2 = area_m2,
    longitude = df_centroids$x,
    latitude = df_centroids$y
  )
  )
}


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
  c = sum(m_adj) / (n*(n-1)/2)
  
  # Calculate PC (Probability of connectivity - area-weighted) 
  total_area <- sum(df_patch$area_m2)
  groups <- components(g)$membership # finds which patches are in the same component
  group_areas <- tapply(df_patch$area_m2, groups, sum)
  PC <- sum(group_areas^2) / (total_area^2)

  
  # Modularity
  modules <- cluster_louvain(g) # switch to different function if needed
  m <- modularity(modules)
  members <- membership(modules) # returns memberships based on modularity --> functional communities
  
  # Centrality metrics
  c_betwenness <- betweenness(g, directed=FALSE)
  c_closeness <- closeness(g)
  deg_list <- degree(g)
  
  # Compute components
  comp <- components(g)
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
    membership = members,
    betweenness = c_betwenness,
    closeness = c_closeness,
    degree = deg_list,
    components = comp,
    comp_number = comp$no,
    largest_comp = largest_size,
    fraction_comp = fraction,
    percent_connected = percent_connected,
    graph = g,
    prob_conn = PC
  ))
}

#' Function for plotting graph 
#' g --> graph from the euclidean_network_e2e function
#' d --> distance used for the network creation
graph_plot_fun <- function(g, d, title, colour = "black"){ 
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
       vertex.color = colour,
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
  prob_connectivity = full_results$prob_conn,
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
    mod_rand <- cluster_louvain(g_rand) #can be changed to other clustering functions - for comparison better to match with same as in euclidean_network_e2e function
    random_modularity[i] <- modularity(mod_rand)
  }
  return(random_modularity)
}
