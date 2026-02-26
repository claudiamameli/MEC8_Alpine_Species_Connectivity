# Library -----------------------------------------------------------------
# Importing functions and full library
source("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/functions_library.R")
set.seed(54367) # Setting seed for random calculation to return the same values 

# Import node dataframes ------------------------------------------------------------
graphs_path <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/species_specific_df/"

# Selecting 3 m buffer maps only
graph_df <- list.files(graphs_path, pattern = "buff3.csv$", full.names = TRUE)
graph_df_names <- tools::file_path_sans_ext(basename(graph_df))


list2env(setNames(lapply(graph_df, read.csv),
                  graph_df_names),
         envir = .GlobalEnv)


# 1. Node Visualisation and Distribution of area values  -------------------------------------------------
plot_nodes(gnap_pres_buff3, "Present", "3")
# plot_nodes(gnap_fut_buff3, "Future, "3") # to be created



g_supinum_buff3 <- bind_rows(
  gnap_pres_buff3 %>% mutate(scenario = "Present", buffer = "3")
  # , gnap_fut_buff3 %>% mutate(scenario = "Future", buffer = "3")
)

# Order scenarios for visualisation
g_supinum_buff3 <- g_supinum_buff3 %>%
  mutate(scenario = factor(scenario, levels = c("Present"
                                                # , "Future"
                                                )))


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
gnap_pres_net <- euclidean_network_e2e(50, gnap_pres_buff3)
# future need to be added <<<<<<<<<



### a. Network Visualisation  --------------------------------------------
graph_plot_fun(gnap_pres_net$graph, 50, "Present")


### b. Network Metrics ---------------------------------------------------------
gnap_net_metrics <- data.frame(
  rbind(net_metrics_fun(gnap_pres_net, "Present") # add future when you have it <<<<<<
))

### c. Node Metrics ---------------------------------------------------------
gnap_pres_node_metrics <- node_metrics_fun(gnap_pres_net, gnap_pres_buff3) 
# add future when you have it <<<<<<



# 3. Network Analysis -----------------------------------------------------
## 3.1 Modularity ----------------------------------------------------

# a. Present -----------------------------------------------------------
pres_random_mod <- random_mod_cal_fun(gnap_pres_net$graph, 1000)

hist(pres_random_mod, breaks = 10, col = "#A4C9EC",
     main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
     xlab = "Modularity",
     xlim = c(round(min(pres_random_mod),1), na.rm = TRUE),
              max(c(max(pres_random_mod), gnap_net_metrics$modularity), na.rm = TRUE))

abline(v = gnap_net_metrics$modularity, col = "#FF4000", lwd = 2)
legend("topright", legend = "Real Network Modularity", col = "#FF4000", lwd = 2)

# Calculate Z-score
z_score <- (gnap_net_metrics$modularity - mean(pres_random_mod, na.rm = TRUE)) / sd(pres_random_mod, na.rm = TRUE)
print(paste("Z-score (Degree controlled):", round(z_score, 3)))

# Visualise memberships
graph_plot_fun(gnap_pres_net$graph, 50, "Present", gnap_pres_net$membership)

# b. Future -----------------------------------------------------------
# add future when you have it <<<<<<

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
mean(degree(g_gnap_pres))

hist(degree(g_gnap_pres), xlim = c(0, max(degree(g_gnap_pres))), breaks = 30)
text(x = max(degree(g_gnap_pres)),
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
mean(degree(random_net))



#Visualisation 
hist(degree(random_net), xlim = c(0, max(degree(random_net))), breaks = 30)
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
m_estimate <- round(ecount(g_gnap_pres) / vcount(g_gnap_pres))


scalefree_net <- sample_pa(n = vcount(g_gnap_pres), m = m_estimate, directed = FALSE, zero.appeal = 0)

vcount(scalefree_net)
ecount(scalefree_net)
transitivity(scalefree_net) 
components(scalefree_net)
mean(degree(scalefree_net))

#Visualisation 

hist(degree(scalefree_net), xlim = c(0, max(degree(scalefree_net))), breaks = 30)
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

# Create table with transitivity etc <<<<<< 