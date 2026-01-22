# Library -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
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



# TO DEAL WITH ------------------------------------------------------------



quartz()

ggplot(data=g_supinum_pres_buff3, 
       aes(x=longitude, y=latitude, 
           #size=area_m2
       )) +
  geom_point() +
  #geom_text(aes(label=patch_ID), hjust=-0.1, vjust=-0.1, col="darkgrey", size=4) +
  coord_fixed() +
  #labs(size="Patch area\n(m2)") +
  theme_bw() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_line(colour="grey"),
        legend.position="bottom")


quartz()
ggplot(data=merged_patch_info, 
       aes(x=longitude, y=latitude, 
           #size=area_m2
       )) +
  geom_point() +
  #geom_text(aes(label=patch_ID), hjust=-0.1, vjust=-0.1, col="darkgrey", size=4) +
  coord_fixed() +
  labs(size="Patch area\n(m2)") +
  theme_bw() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_line(colour="grey"),
        legend.position="bottom")







ggplot(data=fut_patch_info, 
       aes(x=longitude, y=latitude, 
           #size=area_m2
       )) +
  geom_point() +
  #geom_text(aes(label=patch_ID), hjust=-0.1, vjust=-0.1, col="darkgrey", size=4) +
  coord_fixed() +
  labs(size="Patch area\n(m2)") +
  theme_bw() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_line(colour="grey"),
        legend.position="bottom")


ggplot(data=fut_merged_patch_info, 
       aes(x=longitude, y=latitude, 
           #size=area_m2
       )) +
  geom_point() +
  #geom_text(aes(label=patch_ID), hjust=-0.1, vjust=-0.1, col="darkgrey", size=4) +
  coord_fixed() +
  labs(size="Patch area\n(m2)") +
  theme_bw() +
  theme(axis.title=element_blank(), axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_line(colour="grey"),
        legend.position="bottom")


# edge to edge function  --------------------------------------------------
df_patch <- merged_patch_info
d <- 2 

euclidean_network_e2e <- function(d, df_patch) {
  # matrix of distances between all pairs of patches
  mat_dist <- as.matrix(dist(dplyr::select(df_patch, longitude,latitude), method="euclidean", diag=TRUE, upper=TRUE))
  
  # Calculate and extrapolates the radius of individual patches
  df_patch$radius <- sqrt(df_patch$area_m2 / pi)
  radius_list <- dplyr::select(df_patch, radius)
  
  
  # creation  of a new distance matrix considering edge to edge distance
  e2e_mat_dist <- mat_dist
  e2e_mat_dist[] <- 0 
  print('preliminary matrix created')
  
  #this takes a long time <<<<<<<<<<<<<<<<<<< Changed into the version below
  # for (i in 1:nrow(mat_dist)) {
  #   for (j in 1:ncol(mat_dist)) {
  #     e2e_mat_dist[i, j] <- mat_dist[i, j] - (radius_list$radius[i] + radius_list$radius[j])
  #   }
  # }
  
  # this avoids loops, which was taking too long. 
  radius_matrix <- outer(radius_list$radius, radius_list$radius, `+`)
  e2e_mat_dist <- mat_dist - radius_matrix
  print('edge to edge created')
  
  
  # This is to make sure diagonal and negative values are set to 0 
  diag(e2e_mat_dist) <- 0 
  # sum(is.na(e2e_mat_dist))
  # sum(e2e_mat_dist < 0)
  e2e_mat_dist[e2e_mat_dist < 0] <- 0
  # sum(e2e_mat_dist == 0.1)
  
  # adjacency matrix - 1 = link between patch i and j; 0 = no link between patch i and j
  # assign 1 if distance between patch i and j is less or equal to threshold distance
  # assign 0 if distance between patch i and j is greater than threshold distance
  m_adj <- e2e_mat_dist
  m_adj[e2e_mat_dist <= d] <- 1
  m_adj[e2e_mat_dist > d] <- 0
  diag(m_adj) <- 0 # set diagonal to 0
  # sum(m_adj == 1)
  # sum(m_adj == 0)
  # sum(diag(m_adj))
  print('Matrix cleaned up')
  
  # Connectance (number of realised links / number of all possible links)
  n = nrow(m_adj) # number of patches
  c = sum(m_adj) / (n*(n-1))
  
  
  # create dataframe with connections between patches
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
  
  # rename site to name (for igraph)
  nodes <- df_patch %>% rename(name = patch)
  
  # create igraph object
  g <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  
  # Modularity
  modules <- cluster_infomap(g) #switch to different function if needed
  m <- modularity(modules)
  
  # Centrality metrics
  c_betwenness <- betweenness(g, directed=FALSE)
  c_closeness <- closeness(g)
  # c_bridge <- bridge(g, communities = modules) #takes a long time <<<<<<<<<
  deg_list <- degree(g)
  
  # Compute components
  comp <- components(g, mode = "weak")  # Use mode = "strong" for directed Strongly Connected Component (SCC)
  largest_size <- max(comp$csize)  
  total_nodes <- vcount(g)   
  fraction <- largest_size / total_nodes  # Fraction of nodes in largest component
  
  # Calculate percentage of nodes connected
  connected_nodes <- sum(deg_list >= 1)
  percent_connected <- (connected_nodes / nrow(m_adj)) * 100
  
  return(list(
    edge_index = edge_index,
    distance = d, 
    connectance = c, 
    modularity = m,
    betweenness = c_betwenness,
    closeness = c_closeness,
    #bridge = c_bridge,
    degree = deg_list,
    comp_number = comp$no,
    largest_comp = largest_size,
    fraction_comp = fraction,
    percent_connected = percent_connected,
    graph = g
  ))
}

object.size(merged_patch_info) %>% format(units = "Mb")


trial_net <- euclidean_network_e2e(2, merged_patch_info)  
fut_trial_net <- euclidean_network_e2e(2, fut_merged_patch_info)  


trial_net$percent_connected
fut_trial_net$percent_connected

# plot igraph object
# define node coordinates (igraph layout)

g <- trial_net$graph
nodes <- igraph::as_data_frame(g, what = "vertices")

layout_coords <- nodes %>%
  dplyr::select(name, longitude, latitude) %>%
  filter(name %in% V(g)$name) %>%
  arrange(match(name, V(g)$name)) %>%
  dplyr::select(longitude, latitude) %>% 
  as.matrix()

plot(g,
     layout = layout_coords,
     vertex.size = 0.3,
     #vertex.label = V(g)$name,
     vertex.label = "",
     vertex.color = "black",
     edge.color = "red",
     edge.width = 2,
     main = paste0("Pres Spatial Network, d = ",d) )


g2 <- fut_trial_net$graph
nodes <- igraph::as_data_frame(g2, what = "vertices")

layout_coords <- nodes %>%
  dplyr::select(name, longitude, latitude) %>%
  filter(name %in% V(g2)$name) %>%
  arrange(match(name, V(g2)$name)) %>%
  dplyr::select(longitude, latitude) %>% 
  as.matrix()

plot(g2,
     layout = layout_coords,
     vertex.size = 0.3,
     #vertex.label = V(g2)$name,
     vertex.label = "",
     vertex.color = "black",
     edge.color = "red",
     edge.width = 2,
     main = paste0("Fut Spatial Network, d = ",d) )


pres_dist_metrics <-  data.frame(
  distance = trial_net$distance,
  percent_connected = trial_net$percent_connected,
  connectance = trial_net$connectance,
  modularity = trial_net$modularity,
  comp_number = trial_net$comp_number, 
  comp_largest = trial_net$largest_comp,
  comp_fraction = trial_net$fraction_comp
)
fut_dist_metrics <-  data.frame(
  distance = fut_trial_net$distance,
  percent_connected = fut_trial_net$percent_connected,
  connectance = fut_trial_net$connectance,
  modularity = fut_trial_net$modularity,
  comp_number = fut_trial_net$comp_number, 
  comp_largest = fut_trial_net$largest_comp,
  comp_fraction = fut_trial_net$fraction_comp
)

pres_dist_metrics
fut_dist_metrics

