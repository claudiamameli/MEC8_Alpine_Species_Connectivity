# Library -----------------------------------------------------------------
# Importing functions and full library
source("functions_library.R")
set.seed(54367) # Setting seed for random calculation to return the same values 

# Import node dataframes ------------------------------------------------------------
node_df <- list.files("node_df", full.names = TRUE)
df_names <- tools::file_path_sans_ext(basename(node_df))


list2env(setNames(lapply(node_df, read.csv),
                  df_names),
         envir = .GlobalEnv)


# 1. Node Visualisation and Distribution of area values  -------------------------------------------------
plot_nodes(g_supinum_pres, "Present", "Gnaphalium Supinum")
plot_nodes(g_supinum_fut245, "SSP2-4.5", "Gnaphalium Supinum")
plot_nodes(g_supinum_fut585, "SSP5-8.5", "Gnaphalium Supinum")



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


# 2. Network creation & Metric calculation  --------------------------------------------------
## 2.1 G. supinum networks ------------------------------------------------------
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
        net_metrics_fun(gnap_ssp245_net, "fut_ssp245"),
        net_metrics_fun(gnap_ssp585_net, "fut_ssp585")
))

### c. Node Metrics ---------------------------------------------------------
gnap_pres_node_metrics <- node_metrics_fun(gnap_pres_net, g_supinum_pres) 
gnap_fut245_node_metrics <- node_metrics_fun(gnap_ssp245_net, g_supinum_fut245) 
gnap_fut585_node_metrics <- node_metrics_fun(gnap_ssp585_net, g_supinum_fut585) 


# 3. Network Analysis -----------------------------------------------------
## 3.1 Modularity ----------------------------------------------------
### a. Present -----------------------------------------------------------
pres_random_mod <- random_mod_cal_fun(gnap_pres_net$graph, 1000)
pres_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "pres"]

hist(pres_random_mod, breaks = 10, col = "#A4C9EC",
     main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
     xlab = "Modularity",
     xlim = c(round(min(pres_random_mod),1), na.rm = TRUE),
              max(c(max(pres_random_mod), pres_mod), na.rm = TRUE))

abline(v = pres_mod, col = "#FF4000", lwd = 2) # <<<<<<<<< subset scenario 
legend("topright", legend = "Real Network Modularity", col = "#FF4000", lwd = 2)

# Calculate Z-score
pres_z_score <- (pres_mod - mean(pres_random_mod, na.rm = TRUE)) / sd(pres_random_mod, na.rm = TRUE)
print(paste("Present Z-score (Degree controlled):", round(pres_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_pres_net$graph, 50, "Present", gnap_pres_net$membership)

### b. Future ssp245 -----------------------------------------------------------
fut245_random_mod <- random_mod_cal_fun(gnap_ssp245_net$graph, 1000)
fut245_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "fut_ssp245"]

hist(fut245_random_mod, breaks = 10, col = "#A4C9EC",
     main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
     xlab = "Modularity",
     xlim = c(round(min(fut245_random_mod),1), na.rm = TRUE),
     max(c(max(fut245_random_mod), fut245_mod), na.rm = TRUE))

abline(v = fut245_mod, col = "#FF4000", lwd = 2) 
legend("topright", legend = "Real Network Modularity", col = "#FF4000", lwd = 2)

### Calculate Z-score
fut245_z_score <- (fut245_mod - mean(fut245_random_mod, na.rm = TRUE)) / sd(fut245_random_mod, na.rm = TRUE)
print(paste("Fut 245 Z-score (Degree controlled):", round(fut245_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp245_net$graph, 50, "ssp245", gnap_ssp245_net$membership)

### b. Future ssp585 -----------------------------------------------------------
fut585_random_mod <- random_mod_cal_fun(gnap_ssp585_net$graph, 1000)
fut585_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "fut_ssp585"]

hist(fut585_random_mod, breaks = 10, col = "#A4C9EC",
     main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
     xlab = "Modularity",
     xlim = c(round(min(fut585_random_mod),1), na.rm = TRUE),
     max(c(max(fut585_random_mod), fut585_mod), na.rm = TRUE))

abline(v = fut585_mod, col = "#FF4000", lwd = 2) # <<<<<<<<< subset scenario 
legend("topright", legend = "Real Network Modularity", col = "#FF4000", lwd = 2)

### Calculate Z-score
fut585_z_score <- (fut585_mod - mean(fut585_random_mod, na.rm = TRUE)) / sd(fut585_random_mod, na.rm = TRUE)
print(paste("Fut 585 Z-score (Degree controlled):", round(fut585_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp585_net$graph, 50, "ssp585", gnap_ssp585_net$membership)

## 3.2 Degree distribution ------------------------------------
# This section aims to understand which degree distribution our network has to classify if the networks are closer to a random, scale-free or small-world network.

## 3.2.1 Present -----------------------------------------------------------------
### > Original network --------------------------------------------------------
g_gnap_pres <- gnap_pres_net$graph

vcount(g_gnap_pres)
ecount(g_gnap_pres)
transitivity(g_gnap_pres) # clustering coefficient
transitivity(g_gnap_pres, "local")  # Per-patch
mean(transitivity(g_gnap_pres, "local"), na.rm=TRUE)

components(g_gnap_pres)
mean(igraph::degree(g_gnap_pres))

hist(igraph::degree(g_gnap_pres), xlim = c(0, max(igraph::degree(g_gnap_pres))), breaks = 100)
text(x = max(igraph::degree(g_gnap_pres)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)


### > Random network ----------------------------------------------------------
random_net <- sample_correlated_gnp(g_gnap_pres, p = edge_density(g_gnap_pres), corr = 0.6)

vcount(random_net)
ecount(random_net)
transitivity(random_net) 
components(random_net)
mean(igraph::degree(random_net))



#Visualisation 
hist(igraph::degree(random_net), xlim = c(0, max(igraph::degree(random_net))), breaks = 100)
abline(v = mean(igraph::degree(g_gnap_pres)), col = "#FF4000", lwd = 2, lty = 2)
text(x = 45,
     y = par("usr")[4] * 0.9,
     labels = paste("Original Network \n Mean = ", round(mean(igraph::degree(g_gnap_pres)), 2)),
     col = "#FF4000",
     cex = 1)
text(x = max(igraph::degree(random_net)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)


### > Scale-free network ------------------------------------------------------
m_estimate <- round(ecount(g_gnap_pres) / vcount(g_gnap_pres))


scalefree_net <- sample_pa(n = vcount(g_gnap_pres), m = m_estimate, directed = FALSE, zero.appeal = 0)

vcount(scalefree_net)
ecount(scalefree_net)
transitivity(scalefree_net) 
components(scalefree_net)
mean(igraph::degree(scalefree_net))

#Visualisation 

hist(igraph::degree(scalefree_net), xlim = c(0, max(igraph::degree(scalefree_net))), breaks = 100)
abline(v = mean(igraph::degree(g_gnap_pres)), col = "#FF4000", lwd = 2, lty = 2)
text(x = 20,
     y = par("usr")[4] * 0.9,
     labels = paste0("Original Network \n Mean = ", round(mean(igraph::degree(g_gnap_pres)), 2)),
     col = "#FF4000",
     cex = 1)
text(x = max(igraph::degree(scalefree_net)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)

# Create table with transitivity etc <<<<<< 