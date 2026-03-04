# Library -----------------------------------------------------------------
library(caseconverter) 
library(colorRamps)
library(cowplot)
library(cluster)
library(foreach)
library(furrr)
library(gam)
library(ggplot2)
library(ggridges)
library(gpkg)
library(igraph)
library(irr)
library(networktools)
library(mgcv)
library(paletteer)
library(patchwork)
library(PresenceAbsence)
library(pROC)
library(randomForest)
library(raster)
library(RColorBrewer)
library(readxl)
library(rsample)
library(rworldmap)
library(sf)
library(sp)
library(tictoc)
library(terra)
library(tidymodels)
library(tidyverse)


# Read file functions -----------------------------------------------------
# read_function = eg. read.csv or readRDS
read_file_fun <- function(path, read_function){
  list_ <- list.files(path, full.names = TRUE)
  names <- tools::file_path_sans_ext(basename(list_))
  
  list2env(setNames(lapply(list_, read_function),
                    names),
           envir = .GlobalEnv)
}

# SDMs Functions ----------------------------------------------------------
#' The followings fit_* models are teh basis for the full model evaluation function
#' Make sure the predictors are well defined within the 3 functions and their nomenclature matches that of the data frame provided
#' df is the dataframe with all presence/absence values (the response variable) and values for the predictors (the explanatory variables)

fit_gam <- function(df) {
  gam_model <- bam(
    value ~ s(temperature) +
      s(snow) +
      s(carbon) +
      s(SLP_1) + 
      s(HSD_1) +
      s(TWI_25) +
      s(TPI_3),
    data = df, 
    family = binomial
  )
  return(gam_model)
}

fit_glm <- function(df) {
  glm(
    value ~
      poly(temperature, 2) +
      poly(snow, 2) +
      poly(carbon, 2) +
      poly(SLP_1, 2) +
      poly(HSD_1, 2) +
      poly(TWI_25, 2) +
      poly(TPI_3, 2),
    data = df,
    family = binomial
  )
}

fit_rf <- function(df) {
  rf_model <- randomForest(
    value ~ temperature + 
      snow + 
      carbon +
      SLP_1 + 
      HSD_1 + 
      TWI_25 + 
      TPI_3,
    data = df |> mutate(value = factor(value)), 
    family = binomial()
  )
  return(rf_model)
}

# Compute TSS (True Skill Statistic) values - measures the performance of the models
get_tss <- function(m, df) {
  if (class(m)[2] == "randomForest") {
    prediction <- data.frame(predict(m, df, type = "prob"))$X1
  } else {
    prediction <- predict(m, df, type = "response")
  }
  
  auc <- as.numeric(roc(df$value, prediction, quiet = T)$auc)
  vals <- coords(roc(df$value, prediction, quiet = T), "best")[1, ]
  predicted_binary <- ifelse(prediction > vals$threshold[[1]], 1, 0)
  accuracy <- mean(predicted_binary == df$value)
  tss <- vals$sensitivity[1] + vals$specificity[1] - 1
  
  
  stats <- tibble(
    predicted = prediction,
    observed = df$value,
    TSS = tss,
    sensitivity = vals$sensitivity,
    specificity = vals$specificity,
    threshold = vals$threshold,
    accuracy = accuracy,
    auc = auc
  )
  return(stats)
}



# Compute final TSS for all models
get_final_tss <- function(df) {
  auc <- as.numeric(roc(df$observed, df$predicted, quiet = T)$auc)
  vals <- coords(roc(df$observed, df$predicted, quiet = T), "best")
  predicted_binary <- ifelse(df$predicted > vals$threshold[[1]], 1, 0)
  accuracy <- mean(predicted_binary == df$observed)
  tss <- vals$sensitivity[1] + vals$specificity[1] - 1
  
  out <- tibble(
    threshold = vals$threshold[[1]],
    TSS = tss, 
    sensitivity = vals$sensitivity,
    specificity = vals$specificity,
    auc = auc,
    accuracy = accuracy
  )
  return(out)
}

# Outputs statistics for the models calculating thresholds for the current data
# input_data = the dataframe with presence/absence data and extracted values for all predictors at specific GEOM point
get_model_stats <- function(input_data) {
  
  vfold_cv(input_data, strata = value, v = 5) |>
    mutate(output = map(splits, \(mdf){
      list(
        get_tss(fit_glm(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "glm"),
        get_tss(fit_gam(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "gam"),
        get_tss(fit_rf(analysis(mdf)), assessment(mdf)) |>
          mutate(model = "rf")
      ) |>
        bind_rows()
    })) |>
    select(output) |>
    unnest() -> out
  
  out |>
    group_by(model) |>
    mutate(id = row_number()) |>
    ungroup() |>
    group_by(id) |>
    summarise(
      predicted = weighted.mean(predicted, w = TSS),
      observed = mean(observed)
    ) |>
    get_final_tss() |>
    mutate(model = "ensemble") -> final_TSS
  
  out |>
    summarise(across(TSS:auc, mean), .by = model) |>
    bind_rows(final_TSS) -> stats
  
  stats
}


# Computes ensemble prediction for current and future conditions and scenarios
# input_data = same as above
# rs_df = raster dataframe - with all predictor values
# stats = stats produced from the get_model_stats() output

get_ensemble_prediction <- function(input_data, rs_df, stats) {

  # glm
  glm_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_glm(input_data), rs_df, type = "response"))
  glm_prob <- rast(glm_prediction, crs = "EPSG:3035") |> setNames(paste("GLM"))
  glm_class <- ifel(glm_prob > stats$threshold[stats$model == "glm"], 1, 0)
  
  # gam
  m_gam <- fit_gam(input_data)
  prediction <- rs_df |>
    select(-x, -y) |>
    group_by(chunk = ntile(row_number(), 32)) |> #
    group_split() |>
    future_map(\(x) predict(m_gam, x, type = "response")) |>
    unlist() |>
    as_tibble()
  
  gam_prob <- bind_cols(rs_df[c("x", "y")], prediction) |>
    rast(crs = "EPSG:3035") |>
    setNames("GAM")
  gam_class <- ifel(gam_prob > stats$threshold[stats$model == "gam"], 1, 0)
  
  # rf
  rf_prediction <- rs_df[c("x", "y")] |> bind_cols(predict(fit_rf(input_data), rs_df, type = "prob")[, 2])
  rf_prob <- rast(rf_prediction, crs = "EPSG:3035") |> setNames("RF")
  rf_class <- ifel(rf_prob > stats$threshold[stats$model == "rf"], 1, 0)
  
  ensemble_prediction <- terra::weighted.mean(
    x = c(
      glm_prob,
      gam_prob,
      rf_prob
    ),
    w = c(
      stats$TSS[stats$model == "glm"],
      stats$TSS[stats$model == "gam"],
      stats$TSS[stats$model == "rf"]
    )
  )
  
  ensemble_prediction_class <- ifel(ensemble_prediction > stats$threshold[stats$model == "ensemble"], 1, 0) |>
    setNames("ensemble_prediction")
  
  return(ensemble_prediction_class)
}


# Graph Functions ---------------------------------------------------------------
# Nodes creation function
node_df_fun <- function(original_map, direction, cell_res, buffer_val = NULL){
  original_map[original_map[] == 0] <- NA
  
  # Buffer addition
  if(is.null(buffer_val) == F){
    cat("Adding buffer \n")
    r_buffered <- as.numeric(terra::buffer(original_map, width = buffer_val))
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
plot_nodes <- function(df, scenario, species) {
  ggplot(df, aes(x = longitude, y = latitude)) +
    geom_point(size = 0.4) +
    coord_fixed() +
    theme_bw() +
    labs(title = paste(species, " - ", scenario, " Scenario")) +
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
  deg_list <- igraph::degree(g)
  
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
    g_rand <- sample_degseq(igraph::degree(g), method = "configuration")
    mod_rand <- cluster_louvain(g_rand) #can be changed to other clustering functions - for comparison better to match with same as in euclidean_network_e2e function
    random_modularity[i] <- modularity(mod_rand)
  }
  return(random_modularity)
}

# Modularity Z-score
mod_z_score <- function(random_mod, real_mod){
  (real_mod - mean(random_mod, na.rm = TRUE)) / sd(random_mod, na.rm = TRUE)
  }

# Extra Visualisation Functions ---------------------------------------------------------------
plot_changes <- function(current, future, scenario, species) {
  change <- current + 2 * future
  # Map values to categories
  # 0: 0+0*2 = 0 → No presence
  # 1: 1+0*2 = 1 → Lost
  # 2: 0+1*2 = 2 → Gained
  # 3: 1+1*2 = 3 → Stable
  categories <- c("No presence", "Lost", "Gained", "Stable")
  change_cat <- classify(change, cbind(0:3, 0:3))
  values(change_cat) <- categories[values(change_cat) + 1]
  df <- as.data.frame(change_cat, xy = TRUE)
  colnames(df) <- c("x", "y", "change")
  dummy <- tibble(
    x = NA, y = NA,
    change = factor(categories, levels = categories)
  )
  bind_rows(df, dummy) |>
    ggplot(aes(x = x, y = y, fill = factor(change, levels = categories))) +
    geom_raster() +
    annotate(
      "text",
      x = min(df$x), # left
      y = max(df$y), # top
      label = scenario,
      hjust = 0, # align left
      vjust = 1, # align top
      size = 5
    ) +
    scale_fill_manual(values = c(
      "No presence" = "grey90",
      "Lost"        = "red",
      "Gained"      = "green",
      "Stable"      = "blue"
    ), drop = F) +
    coord_equal() +
    theme_bw() +
    labs(
      subtitle = paste0("TSS: ", round(stats_out$TSS[stats_out$model == "ensemble"], 2)),
      fill = "Change category",
      title = paste0(species)
    ) +
    theme(
      legend.position = c(1, 1), legend.justification = c(1, 1),
      legend.background = element_blank(),
      plot.title = element_text(face = "italic"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) ->p
  return(p)
}

hist_rand_degree <- function(random_mod, real_mod, bin, scenario){
  hist(random_mod, breaks = bin, col = "#A4C9EC",
       main = "Distribution of Modularity\n (Degree-controlled Random Graphs)",
       xlab = "Modularity",
       xlim = c(round(min(random_mod),1), na.rm = TRUE),
       max(c(max(random_mod), real_mod), na.rm = TRUE))
  abline(v = real_mod, col = "#FF4000", lwd = 2)
  legend("topright", legend = paste0(scenario," Network Modularity"), col = "#FF4000", lwd = 2)
  }
