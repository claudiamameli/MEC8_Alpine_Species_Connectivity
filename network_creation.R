# Library -----------------------------------------------------------------
# Importing functions and full library
source("functions_library.R")
set.seed(54367) # Setting seed for random calculations

# 1. Import maps ------------------------------------------------------
files_rast <- list.files("Data/map_predictions" ,pattern = "tif$", full.names=T)

list2env(setNames(lapply(files_rast, rast),
                  tools::file_path_sans_ext(basename(files_rast))), envir = .GlobalEnv)

plot(gnaphalium_supinum_current)
plot(gnaphalium_supinum_ssp245)
plot(gnaphalium_supinum_ssp585)

# 2. Graph creation - Nodes ----------------------------------------------------------
## 2.1 Species - specific graphs ----------------------------------------------------------
## > G_supinum  -------------------------------------------------------
# ### a. present ----------------------------------------------------------
# g_supinum_pres <- node_df_fun(gnaphalium_supinum_current,
#                                      direction = 8,
#                                      cell_res = 1,
#                                      buffer_val = 3)
# 
# ### b. Future SSP2-4.5 -------------------------------------------------------
# g_supinum_fut245 <- node_df_fun(gnaphalium_supinum_ssp245,
#                                     direction = 8,
#                                     cell_res = 1,
#                                     buffer_val = 3)
# 
# ### c. Future SSP5-8.5 -------------------------------------------------------
# g_supinum_fut585 <- node_df_fun(gnaphalium_supinum_ssp585,
#                                       direction = 8,
#                                       cell_res = 1,
#                                       buffer_val = 3)
# 
# 3. Node visualisation --------------------------------
# Load all previously saved dataframes
nodes_df <- list.files("node_df", pattern = "\\.csv$", full.names = TRUE)
df_names <- tools::file_path_sans_ext(basename(nodes_df))

list2env(setNames(lapply(nodes_df, read.csv),
                  df_names),
         envir = .GlobalEnv)


plot_nodes(g_supinum_pres, "Present", "Gnaphalium Supinum")
plot_nodes(g_supinum_fut245, "SSP2-4.5", "Gnaphalium Supinum")
plot_nodes(g_supinum_fut585, "SSP5-8.5", "Gnaphalium Supinum")

# Distribution of area values
g_supinum_all <- bind_rows(
  g_supinum_pres %>% mutate(scenario = "Present"), 
  g_supinum_fut245 %>% mutate(scenario = "ssp245"),
  g_supinum_fut585 %>% mutate(scenario = "ssp585") %>% 
    mutate(scenario = factor(scenario, levels = c("Present", "ssp245", "ssp585")))
)


ggplot(g_supinum_all, aes(x = scenario, y = log(area_m2))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Scenario",
    y = "Patch area - log(m²)",
    title = "Patch area distribution"
  )


# 4. Network creation & Metric calculation  --------------------------------------------------
## 4.1 G. supinum networks ------------------------------------------------------
gnap_pres_net <- euclidean_network_e2e(50, g_supinum_pres)
gnap_ssp245_net <- euclidean_network_e2e(50, g_supinum_fut245)
gnap_ssp585_net <- euclidean_network_e2e(50, g_supinum_fut585)



### a. Network Visualisation  --------------------------------------------
graph_plot_fun(gnap_pres_net$graph, 50, "Present")
graph_plot_fun(gnap_ssp245_net$graph, 50, "Future spp2-4.5")
graph_plot_fun(gnap_ssp585_net$graph, 50, "Future ssp5-8.5")


### b. Network Metrics ---------------------------------------------------------
gnap_net_metrics <- data.frame(
  rbind(net_metrics_fun(gnap_pres_net, "pres"),
        net_metrics_fun(gnap_ssp245_net, "ssp245"),
        net_metrics_fun(gnap_ssp585_net, "ssp585")
  ))

### c. Node Metrics ---------------------------------------------------------
gnap_pres_node_metrics <- node_metrics_fun(gnap_pres_net, g_supinum_pres) 
gnap_fut245_node_metrics <- node_metrics_fun(gnap_ssp245_net, g_supinum_fut245) 
gnap_fut585_node_metrics <- node_metrics_fun(gnap_ssp585_net, g_supinum_fut585) 


#### >>> Save files -----------------------------------------------------
node_df_list <- list(g_supinum_pres = g_supinum_pres,
                       g_supinum_fut245 = g_supinum_fut245,
                       g_supinum_fut585 = g_supinum_fut585)

for(i in names(node_df_list)){
  write.csv(node_df_list[[i]],
            paste0("graph_df/", i,".csv"),
            row.names = FALSE)
}

net_list <- list(gnap_pres_net = gnap_pres_net, 
                 gnap_ssp245_net = gnap_ssp245_net, 
                 gnap_ssp585_net = gnap_ssp585_net)

for(i in names(net_list)){
  write_rds(net_list[[i]],
            paste0("graph_files/", i,".rds"))
}

graph_files_list <- list(gnap_net_metrics = gnap_net_metrics, 
                         gnap_pres_node_metrics = gnap_pres_node_metrics, 
                         gnap_fut245_node_metrics = gnap_fut245_node_metrics, 
                         gnap_fut585_node_metrics = gnap_fut585_node_metrics)

for(i in names(graph_files_list)){
  write.csv(graph_files_list[[i]],
            paste0("graph_df/", i,".csv"),
            row.names = FALSE)
}
