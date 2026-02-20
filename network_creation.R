# Library -----------------------------------------------------------------
# Importing functions and full library
source("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/functions_library.R")


# Import node dataframes ------------------------------------------------------------
graphs_path <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/species_specific_df/"

# Selecting 3 m buffer maps only
graph_df <- list.files(graphs_path, pattern = "buff3.csv$", full.names = TRUE)
graph_df_names <- tools::file_path_sans_ext(basename(graph_df))


list2env(setNames(lapply(graph_df, read.csv),
                  graph_df_names),
         envir = .GlobalEnv)


# 1. Node Visualisation and Distribution of area values  -------------------------------------------------
plot_nodes(g_supinum_pres_buff3, "Present", "3")
plot_nodes(g_supinum_fut245_buff3, "Future 2-4.5", "3")
plot_nodes(g_supinum_fut585_buff3, "Future 5-8.5", "3")

g_supinum_buff3 <- bind_rows(
  g_supinum_pres_buff3 %>% mutate(scenario = "Present", buffer = "3"),
  g_supinum_fut245_buff3 %>% mutate(scenario = "Future 2–4.5", buffer = "3"),
  g_supinum_fut585_buff3 %>% mutate(scenario = "Future 5–8.5", buffer = "3")
)

# Order scenarios for visualisation
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


# 2. Network creation & Metric calculation  --------------------------------------------------
# See "functions_script" for full function


## 2.1 G. supinum networks ------------------------------------------------------
gnap_pres_net <- euclidean_network_e2e(50, g_supinum_pres_buff3)      # present
gnap_fut245_net <- euclidean_network_e2e(50, g_supinum_fut245_buff3)  # future 2-4.5
gnap_fut585_net <- euclidean_network_e2e(50, g_supinum_fut585_buff3)  # future 5-5.8



### a. Network Visualisation  --------------------------------------------
graph_plot_fun(gnap_pres_net$graph, 50, "Present")
graph_plot_fun(gnap_fut245_net$graph, 50, "Future 2-4.5")
graph_plot_fun(gnap_fut585_net$graph, 50, "Future 5-5.8")


### b. Network Metrics ---------------------------------------------------------
gnap_net_metrics <- data.frame(
  rbind(net_metrics_fun(gnap_pres_net, "Present"),
  net_metrics_fun(gnap_fut245_net, "Future 2-4.5"),
  net_metrics_fun(gnap_fut585_net, "Future 5-8.5"))
)

### c. Node Metrics ---------------------------------------------------------
gnap_pres_node_metrics <- node_metrics_fun(gnap_pres_net, g_supinum_pres_buff3) 
gnap_fut245_node_metrics <- node_metrics_fun(gnap_fut245_net, g_supinum_fut245_buff3) 
gnap_fut585_node_metrics <- node_metrics_fun(gnap_fut585_net, g_supinum_fut585_buff3) 



# 3. Network Analysis -----------------------------------------------------


## 3.1 Modularity   ----------------------------------------------------
trial_mod <- random_mod_cal_fun(gnap_pres_net$graph, 1000)

hist(trial_mod, breaks = 20, col = "#A4C9EC",
     main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
     xlab = "Modularity",
     xlim = c(min(0.2, na.rm = TRUE),
              max(c(trial_mod, gnap_pres_net$modularity), na.rm = TRUE))
)
abline(v = gnap_pres_net$modularity, col = "#FF4000", lwd = 2)
legend("topright", legend = "Real Network Modularity", col = "#FF4000", lwd = 2)

# Calculate Z-score
z_score <- (gnap_pres_net$modularity - mean(trial_mod, na.rm = TRUE)) / sd(trial_mod, na.rm = TRUE)
print(paste("Z-score (Degree controlled):", round(z_score, 5)))


## 3.2 Degree distribution ------------------------------------
# This section aims to understand which degree distribution our network has to classify if the networks are closer to a random, scale-free or small-world network.
### > Original network --------------------------------------------------------
g_gnap_pres <- gnap_pres_net$graph

g_poll_plot()
vcount(g_gnap_pres)
ecount(g_gnap_pres)
transitivity(g_gnap_pres) # clustering coefficient
components(g_gnap_pres)
mean(degree(g_gnap_pres))

hist(degree(g_gnap_pres), xlim = c(0, 300), breaks = 60)
text(x = max(degree(g_gnap_pres)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)

### > Random network ----------------------------------------------------------
set.seed(54367) # always run this before creating random network so that we get same results
random_net <- sample_correlated_gnp(g_gnap_pres, p = edge_density(g_gnap_pres), corr = 0.6)

vcount(random_net)
ecount(random_net)
transitivity(random_net) 
components(random_net)
mean(degree(random_net))



#Visualisation 
hist(degree(random_net), xlim = c(0, 300), breaks = 60)
abline(v = mean(degree(g_gnap_pres)), col = "#FF4000", lwd = 2, lty = 2)
text(x = 45,
     y = par("usr")[4] * 0.9,
     labels = paste("Original Network \n Mean = ", round(mean(degree(g_gnap_pres)), 2)),
     col = "#FF4000",
     cex = 1)
text(x = max(degree(random_net)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)


### > Scale-free network ------------------------------------------------------
m_estimate <- ceiling(ecount(g_gnap_pres) / vcount(g_gnap_pres))

set.seed(54367)
scalefree_net <- sample_pa(n = vcount(g_gnap_pres), m = m_estimate, directed = FALSE, zero.appeal = 0)
vcount(scalefree_net)
ecount(scalefree_net)
transitivity(scalefree_net) 
components(scalefree_net)
mean(degree(scalefree_net))
?sample_pa
#Visualisation 

hist(degree(scalefree_net), xlim = c(0, max(degree(scalefree_net))))
abline(v = mean(degree(g_gnap_pres)), col = "#FF4000", lwd = 2, lty = 2)
text(x = 20,
     y = par("usr")[4] * 0.9,
     labels = paste("Original Network \n Mean = ", round(mean(degree(g_gnap_pres)), 2)),
     col = "#FF4000",
     cex = 1)
text(x = max(degree(scalefree_net)),
     y = 1,
     labels = "*",
     col = "black",
     cex = 4)


### > Small-world network -----------------------------------------------------
sw_net <- sample_smallworld(dim = 1, size = vcount(g_gnap_pres), nei = floor(mean(degree(g_gnap_pres))/2), p = 0.1)

hist(degree(sw_net))

#' Problem is that when we are attaching edges to the scale-free it is impossible to attach less than 1 edge per node, so our histograms might become extremely skewed.
#' Generally, this type of analysis is easier when you have well connected networks



