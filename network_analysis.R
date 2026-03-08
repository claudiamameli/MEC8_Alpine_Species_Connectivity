# Library -----------------------------------------------------------------
# Importing functions and full library
source("functions_library.R")
set.seed(54367) # Setting seed for random calculations

# Import data ------------------------------------------------------------
read_file_fun("graph_df/", read.csv)
read_file_fun("graph_files/", readRDS)

graph_plot_fun(gnap_pres_net$graph, 50, "Present", gnap_pres_net$components$membership)
graph_plot_fun(gnap_ssp245_net$graph, 50, "SSP2-4.5", gnap_ssp245_net$components$membership)
graph_plot_fun(gnap_ssp585_net$graph, 50, "SSP5-8.5", gnap_ssp585_net$components$membership)

# 1.Area vs node degree -----------------------------------------------------
ggplot(gnap_pres_node_metrics, aes(y = log10(degree + 1), x = log10(area_m2))) +
  geom_point()+
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Patch Area log10(m²)", y = "Node Degree log10(n + 1)")+
  theme_bw()

dfs <- list(
  present = gnap_pres_node_metrics,
  fut245  = gnap_fut245_node_metrics,
  fut585  = gnap_fut585_node_metrics
)

for (name in names(dfs)) {
  cat("\n=== Model for", name, "===\n")
  model <- lm(degree ~ area_m2, data = dfs[[name]])
  print(summary(model))
}

for (name in names(dfs)) {
  p <- ggplot(dfs[[name]], aes(y = log10(degree + 1), x = log10(area_m2))) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = name,
      x = "Patch area log10(m²)",
      y = "Node degree log10(n + 1)"
    ) +
    theme_bw()
  
  print(p)
}

# 2. Modularity ----------------------------------------------------
### a. Present -----------------------------------------------------------
pres_random_mod <- random_mod_cal_fun(gnap_pres_net$graph, 1000)
pres_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "pres"]

hist_rand_mod(pres_random_mod, pres_mod, 10, "Present")

# Calculate Z-score
pres_z_score <- mod_z_score(pres_random_mod, pres_mod)
print(paste("Present Z-score (Degree controlled):", round(pres_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_pres_net$graph, 50, "Present", gnap_pres_net$membership)

### b. Future ssp245 -----------------------------------------------------------
fut245_random_mod <- random_mod_cal_fun(gnap_ssp245_net$graph, 1000)
fut245_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "ssp245"]

hist_rand_mod(fut245_random_mod, fut245_mod, 10, "ssp245")

# Calculate Z-score
fut245_z_score <- mod_z_score(fut245_random_mod, fut245_mod)
print(paste("Fut 245 Z-score (Degree controlled):", round(fut245_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp245_net$graph, 50, "ssp245", gnap_ssp245_net$membership)

### c. Future ssp585 -----------------------------------------------------------
fut585_random_mod <- random_mod_cal_fun(gnap_ssp585_net$graph, 1000)
fut585_mod <- gnap_net_metrics$modularity[gnap_net_metrics$scenario == "ssp585"]

hist_rand_mod(fut585_random_mod, fut585_mod, 10, "ssp585")

# Calculate Z-score
fut585_z_score <- mod_z_score(fut585_random_mod, fut585_mod)
print(paste("Fut 585 Z-score (Degree controlled):", round(fut585_z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_ssp585_net$graph, 50, "ssp585", gnap_ssp585_net$membership)

### >>> Table -----------------------------------------------------------
modularity_summary <- data.frame(
  scenario = c("Present", "ssp2-4.5", "ssp5-8.5"),
  nodes      = c(vcount(gnap_pres_net$graph),
                   vcount(gnap_ssp245_net$graph),
                   vcount(gnap_ssp585_net$graph)),
  edges      = c(ecount(gnap_pres_net$graph),
                   ecount(gnap_ssp245_net$graph),
                   ecount(gnap_ssp585_net$graph)),
  modules    = c(length(unique(gnap_pres_net$membership)),
                   length(unique(gnap_ssp245_net$membership)),
                   length(unique(gnap_ssp585_net$membership))),
  modularity = round(c(pres_mod, fut245_mod, fut585_mod), 4),
  z_score = round(c(pres_z_score, fut245_z_score, fut585_z_score),4)
)

modularity_summary

## 3. Degree distribution ------------------------------------
# This section aims to understand which degree distribution our network has to classify if the networks are closer to a random, scale-free or small-world network.

### a. Present -----------------------------------------------------------------
#### > Original network --------------------------------------------------------
g_gnap_pres <- gnap_pres_net$graph

vcount(g_gnap_pres)
ecount(g_gnap_pres)
gnap_pres_net$global_trans
components(g_gnap_pres)
mean(igraph::degree(g_gnap_pres))

# Visualisation
hist_degree(g_gnap_pres, "Present Original Network")

#### > Random network ----------------------------------------------------------
random_net <- sample_correlated_gnp(g_gnap_pres, p = edge_density(g_gnap_pres), corr = 0.6)

vcount(random_net)
ecount(random_net)
transitivity(random_net) 
components(random_net)
mean(igraph::degree(random_net))

# Visualisation 
hist_degree(random_net, "Present Random Network")

#### > Scale-free network ------------------------------------------------------
m_estimate <- round(ecount(g_gnap_pres) / vcount(g_gnap_pres))
scalefree_net <- sample_pa(n = vcount(g_gnap_pres), m = m_estimate, directed = FALSE, zero.appeal = 0)

vcount(scalefree_net)
ecount(scalefree_net)
transitivity(scalefree_net) 
components(scalefree_net)
mean(igraph::degree(scalefree_net))

#Visualisation 
hist_degree(scalefree_net, "Present Scale-free Network")

### b. Future ssp245 -----------------------------------------------------------------
#### > Original network --------------------------------------------------------
g_gnap_ssp245 <- gnap_ssp245_net$graph

vcount(g_gnap_ssp245)
ecount(g_gnap_ssp245)
gnap_ssp245_net$global_trans
components(g_gnap_ssp245)
mean(igraph::degree(g_gnap_ssp245))

# Visualisation
hist_degree(g_gnap_ssp245, "ssp245 Original Network")

#### > Random network ----------------------------------------------------------
random_net_ssp245 <- sample_correlated_gnp(g_gnap_ssp245, p = edge_density(g_gnap_ssp245), corr = 0.6)

vcount(random_net_ssp245)
ecount(random_net_ssp245)
transitivity(random_net_ssp245) 
components(random_net_ssp245)
mean(igraph::degree(random_net_ssp245))

# Visualisation 
hist_degree(random_net_ssp245, "ssp245 Random Network")

#### > Scale-free network ------------------------------------------------------
m_estimate <- round(ecount(g_gnap_ssp245) / vcount(g_gnap_ssp245))
scalefree_net_ssp245 <- sample_pa(n = vcount(g_gnap_ssp245), m = m_estimate, directed = FALSE, zero.appeal = 0)

vcount(scalefree_net_ssp245)
ecount(scalefree_net_ssp245)
transitivity(scalefree_net_ssp245) 
components(scalefree_net_ssp245)
mean(igraph::degree(scalefree_net_ssp245))

#Visualisation 
hist_degree(scalefree_net_ssp245, "ssp245 Scale-free Network")

### c. Future ssp585 -----------------------------------------------------------
#### > Original network --------------------------------------------------------
g_gnap_ssp585 <- gnap_ssp585_net$graph

vcount(g_gnap_ssp585)
ecount(g_gnap_ssp585)
gnap_ssp585_net$global_trans
components(g_gnap_ssp585)
mean(igraph::degree(g_gnap_ssp585))

# Visualisation
hist_degree(g_gnap_ssp585, "ssp585 Original Network")

#### > Random network ----------------------------------------------------------
random_net_ssp585 <- sample_correlated_gnp(g_gnap_ssp585, p = edge_density(g_gnap_ssp585), corr = 0.6)

vcount(random_net_ssp585)
ecount(random_net_ssp585)
transitivity(random_net_ssp585) 
components(random_net_ssp585)
mean(igraph::degree(random_net_ssp585))

# Visualisation 
hist_degree(random_net_ssp585, "ssp585 Random Network")

#### > Scale-free network ------------------------------------------------------
m_estimate <- round(ecount(g_gnap_ssp585) / vcount(g_gnap_ssp585))
scalefree_net_ssp585 <- sample_pa(n = vcount(g_gnap_ssp585), m = m_estimate, directed = FALSE, zero.appeal = 0)

vcount(scalefree_net_ssp585)
ecount(scalefree_net_ssp585)
transitivity(scalefree_net_ssp585) 
components(scalefree_net_ssp585)
mean(igraph::degree(scalefree_net_ssp585))

#Visualisation 
hist_degree(scalefree_net_ssp585, "ssp585 Scale-free Network")

### >>> Table ---------------------------------------------------------------
degree_summary <- data.frame(
  scenario = rep(c("Present", "ssp2-4.5", "ssp5-8.5"), each = 3),
  network_type = rep(c("Original", "Random", "Scale-free"), 3),
  transitivity = c(gnap_pres_net$global_trans, 
                   transitivity(random_net), 
                   transitivity(scalefree_net),
                   gnap_ssp245_net$global_trans, 
                   transitivity(random_net_ssp245), 
                   transitivity(scalefree_net_ssp245),
                   gnap_ssp585_net$global_trans, 
                   transitivity(random_net_ssp585), 
                   transitivity(scalefree_net_ssp585)),
  mean_degree = c(round(mean(igraph::degree(g_gnap_pres)), 2), 
                  round(mean(igraph::degree(random_net)), 2), 
                  round(mean(igraph::degree(scalefree_net)), 2),
                  round(mean(igraph::degree(g_gnap_ssp245)), 2), 
                  round(mean(igraph::degree(random_net_ssp245)), 2), 
                  round(mean(igraph::degree(scalefree_net_ssp245)), 2),
                  round(mean(igraph::degree(g_gnap_ssp585)), 2),
                  round(mean(igraph::degree(random_net_ssp585)), 2), 
                  round(mean(igraph::degree(scalefree_net_ssp585)), 2)),
  max_degree = c(max(igraph::degree(g_gnap_pres)), 
                 max(igraph::degree(random_net)), 
                 max(igraph::degree(scalefree_net)),
                 max(igraph::degree(g_gnap_ssp245)), 
                 max(igraph::degree(random_net_ssp245)), 
                 max(igraph::degree(scalefree_net_ssp245)),
                 max(igraph::degree(g_gnap_ssp585)),
                 max(igraph::degree(random_net_ssp585)), 
                 max(igraph::degree(scalefree_net_ssp585))))

degree_summary



