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

# Importing functions 
source("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/functions_script.R")


# Import node dataframes ------------------------------------------------------------
graphs_path <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/species_specific_df/"
graph_df <- list.files(graphs_path, pattern = "buff3.csv$", full.names = TRUE)

graph_df_names <- tools::file_path_sans_ext(basename(graph_df))


list2env(setNames(lapply(graph_df, read.csv),
                  graph_df_names),
         envir = .GlobalEnv)




# Node Visualisation and Distribution of area values  -------------------------------------------------
plot_nodes(g_supinum_pres_buff3, "Present", "3")
plot_nodes(g_supinum_fut245_buff3, "Future 2-4.5", "3")
plot_nodes(g_supinum_fut585_buff3, "Future 5-8.5", "3")

g_supinum_buff3 <- bind_rows(
  g_supinum_pres_buff3 %>% mutate(scenario = "Present", buffer = "3"),
  g_supinum_fut245_buff3 %>% mutate(scenario = "Future 2–4.5", buffer = "3"),
  g_supinum_fut585_buff3 %>% mutate(scenario = "Future 5–8.5", buffer = "3")
)

# Order scenarios
g_supinum_buff3 <- g_supinum_buff3 %>%
  mutate(scenario = factor(scenario, levels = c("Present", "Future 2–4.5", "Future 5–8.5")))


ggplot(subset(g_supinum_buff3, buffer == "3"),
       aes(x = scenario, y = log(area_m2))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Patch area - log(m²)",
    title = "Patch area distribution"
  )


# Network creation & Metric calculation  --------------------------------------------------
# See "functions_script" for full function


# G. supinum networks ------------------------------------------------------
gnap_pres_net <- euclidean_network_e2e(50, g_supinum_pres_buff3)  
gnap_fut245_net <- euclidean_network_e2e(50, g_supinum_fut245_buff3)  
gnap_fut585_net <- euclidean_network_e2e(50, g_supinum_fut585_buff3)  


gnap_pres_net$percent_connected
gnap_fut245_net$percent_connected
gnap_fut585_net$percent_connected


g <- trial_net$graph
nodes <- igraph::as_data_frame(g, what = "vertices")

plot(g)

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
     main = paste0("Pres Spatial Network, d = ",d) )


g2 <- trial_net_fut$graph
nodes <- igraph::as_data_frame(g2, what = "vertices")

layout_coords <- nodes %>%
  dplyr::select(name, longitude, latitude) %>%
  filter(name %in% V(g2)$name) %>%
  arrange(match(name, V(g2)$name)) %>%
  dplyr::select(longitude, latitude) %>% 
  as.matrix()

plot(g2,
     layout = layout_coords,
     vertex.size = log10(nodes$area_m2),
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
  distance = trial_net_fut$distance,
  percent_connected = trial_net_fut$percent_connected,
  connectance = trial_net_fut$connectance,
  modularity = trial_net_fut$modularity,
  comp_number = trial_net_fut$comp_number, 
  comp_largest = trial_net_fut$largest_comp,
  comp_fraction = trial_net_fut$fraction_comp
)

pres_dist_metrics
fut_dist_metrics

summary_table_scenarios
