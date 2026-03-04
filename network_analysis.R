# Library -----------------------------------------------------------------
# Importing functions and full library
source("functions_library.R")
set.seed(54367) # Setting seed for random calculations

# Import data ------------------------------------------------------------
read_file_fun("graph_df/", read.csv)
read_file_fun("graph_files/", readRDS)



# 3. Network Analysis -----------------------------------------------------
## 3.1 Modularity ----------------------------------------------------
### a. Present -----------------------------------------------------------
pres_random_mod <- random_mod_cal_fun(gnap_pres_net$graph, 1000)
pres_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "pres"]

hist_rand_degree(pres_random_mod, pres_mod, 10, "Present")

# Calculate Z-score
pres_z_score <- mod_z_score(pres_random_mod, pres_mod)
print(paste("Present Z-score (Degree controlled):", round(pres_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_pres_net$graph, 50, "Present", gnap_pres_net$membership)

### b. Future ssp245 -----------------------------------------------------------
fut245_random_mod <- random_mod_cal_fun(gnap_ssp245_net$graph, 1000)
fut245_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "ssp245"]

hist_rand_degree(fut245_random_mod, fut245_mod, 10, "ssp245")

# Calculate Z-score
fut245_z_score <- mod_z_score(fut245_random_mod, fut245_mod)
print(paste("Fut 245 Z-score (Degree controlled):", round(fut245_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp245_net$graph, 50, "ssp245", gnap_ssp245_net$membership)

### b. Future ssp585 -----------------------------------------------------------
fut585_random_mod <- random_mod_cal_fun(gnap_ssp585_net$graph, 1000)
fut585_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "ssp585"]

hist_rand_degree(fut585_random_mod, fut585_mod, 10, "ssp585")

# Calculate Z-score
fut585_z_score <- mod_z_score(fut585_random_mod, fut585_mod)
print(paste("Fut 585 Z-score (Degree controlled):", round(fut585_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp585_net$graph, 50, "ssp585", gnap_ssp585_net$membership)

## 3.2 Degree distribution ------------------------------------
# This section aims to understand which degree distribution our network has to classify if the networks are closer to a random, scale-free or small-world network.

### a. Present -----------------------------------------------------------------
#### > Original network --------------------------------------------------------
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


#### > Random network ----------------------------------------------------------
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


#### > Scale-free network ------------------------------------------------------
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