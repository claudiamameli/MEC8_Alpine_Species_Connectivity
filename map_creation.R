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
library(ggpubr)


# 1.  Import data frames and maps ------------------------------------------------------
path_df <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/prelim_dfs"
files_df <- list.files(path_df, pattern = "\\.csv$", full.names = TRUE)

list2env(setNames(lapply(files_df, read.csv),
    tools::file_path_sans_ext(basename(files_df))), envir = .GlobalEnv)

path_maps <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/final_overlay_maps"
files_rast <- list.files(path_maps ,pattern = "tif$", full.names=T)

list2env(setNames(lapply(files_rast, rast),
                  tools::file_path_sans_ext(basename(files_rast))), envir = .GlobalEnv)

plot(g_supinum_pres_bin)
plot(g_supinum_fut245_bin)
plot(g_supinum_fut585_bin)


# 2. Graph creation ----------------------------------------------------------
## 2.1 Nodes creation function  ----------------------------------------------------------------

node_df_fun <- function(original_map, direction, cell_res, buffer_val = NULL){
  
  # Connect neighbouring cells (direction = 8 for corners, = 4 for sides only)
  cat("Creating patches \n")
  patch_data <- patches(original_map, directions = direction, zeroAsNA=T)
  
  # Calculate cell counts and area
  cat("Counting cells \n")
  patch_cell_count <- freq(patch_data, bylayer = FALSE)  
  area_m2 <- patch_cell_count$count * cell_res
  
  # Create polygons and find centroids
  cat("Creating polygons and centroids \n")
  patch_polygons <- as.polygons(patch_data, dissolve = TRUE)
  df_centroids <- crds(centroids(patch_polygons), df = TRUE)
  
  # Buffer addition
  if(is.null(buffer_val) == F){
    cat("Adding buffer \n")
    patch_polys_buffered <- buffer(patch_polygons, width = buffer_val)
    
    r_buffered <- rasterize(patch_polys_buffered, g_supinum_pres_bin, field = 1)
    
    buffer_patches <- patches(r_buffered, directions = direction)
    patch_polygons <- as.polygons(buffer_patches, dissolve = TRUE)
    
    area_m2 <- expanse(patch_polygons, unit = "m")
    df_centroids <- crds(centroids(patch_polygons), df = TRUE)
  }

  
  cat("Returning data frame \n")
  return(data.frame(
      patch = patch_polygons$patches,
      area_m2 = area_m2,
      longitude = df_centroids$x,
      latitude = df_centroids$y
    )
  )
}



## 2.2 Species - specific graphs ----------------------------------------------------------


### > G_supinum  -------------------------------------------------------
#### a. present ----------------------------------------------------------
# g_supinum_pres_buff0 <- node_df_fun(g_supinum_pres_bin, 
#                                      direction = 8, 
#                                      cell_res = 1, 
#                                      buffer_val = NULL)
# g_supinum_pres_buff1 <- node_df_fun(g_supinum_pres_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 1)
# g_supinum_pres_buff2 <- node_df_fun(g_supinum_pres_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 2)
# g_supinum_pres_buff3 <- node_df_fun(g_supinum_pres_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 3)
# 
#### b. future 245 -------------------------------------------------------
# g_supinum_fut245_buff0 <- node_df_fun(g_supinum_fut245_bin, 
#                                      direction = 8, 
#                                      cell_res = 1, 
#                                      buffer_val = NULL)
# g_supinum_fut245_buff1 <- node_df_fun(g_supinum_fut245_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 1)
# g_supinum_fut245_buff2 <- node_df_fun(g_supinum_fut245_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 2)
# g_supinum_fut245_buff3 <- node_df_fun(g_supinum_fut245_bin, 
#                                     direction = 8, 
#                                     cell_res = 1, 
#                                     buffer_val = 3)
# 
#### c. future 585 -------------------------------------------------------
# g_supinum_fut585_buff0 <- node_df_fun(g_supinum_fut585_bin, 
#                                        direction = 8, 
#                                        cell_res = 1, 
#                                        buffer_val = NULL)
# g_supinum_fut585_buff1 <- node_df_fun(g_supinum_fut585_bin, 
#                                       direction = 8, 
#                                       cell_res = 1, 
#                                       buffer_val = 1)
# g_supinum_fut585_buff2 <- node_df_fun(g_supinum_fut585_bin, 
#                                       direction = 8, 
#                                       cell_res = 1, 
#                                       buffer_val = 2)
# g_supinum_fut585_buff3 <- node_df_fun(g_supinum_fut585_bin, 
#                                       direction = 8, 
#                                       cell_res = 1, 
#                                       buffer_val = 3)
# 
# 
# 
# 
#### >>> Save files --------------------------------------------------------------

df_graphs_list <- list(g_supinum_pres_buff0 = g_supinum_pres_buff0,
                g_supinum_pres_buff1 = g_supinum_pres_buff1,
                g_supinum_pres_buff2 = g_supinum_pres_buff2,
                g_supinum_pres_buff3 = g_supinum_pres_buff3,
                g_supinum_fut245_buff0 = g_supinum_fut245_buff0,
                g_supinum_fut245_buff1 = g_supinum_fut245_buff1,
                g_supinum_fut245_buff2 = g_supinum_fut245_buff2,
                g_supinum_fut245_buff3 = g_supinum_fut245_buff3,
                g_supinum_fut585_buff0 = g_supinum_fut585_buff0,
                g_supinum_fut585_buff1 = g_supinum_fut585_buff1,
                g_supinum_fut585_buff2 = g_supinum_fut585_buff2,
                g_supinum_fut585_buff3 = g_supinum_fut585_buff3)

graphs_path <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/species_specific_df/"


for(i in names(df_graphs_list)){
  write.csv(df_graphs_list[[i]], paste0(graphs_path, i,".csv"))
}


# Visual comparison of different scenarios --------------------------------
summary_table_scenarios <- read.csv("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/species_specific_df/summary_table_scenarios.csv")

summary_table_scenarios


# Change order to change appearance in the Legends.
summary_table_scenarios$scenario <- factor(summary_table_scenarios$scenario, levels = c('Present', 'Future 245', 'Future 585'))


total_patches_plot <- ggplot(data = summary_table_scenarios, 
       aes(x = buffer, y = tot_patches, 
           group = scenario, 
           colour = scenario 
       )) +
  geom_line() +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() + 
  scale_color_manual(values=c("#00AA00", "#FFA500", "#FF4000")) +
  labs(x = "Buffer (m)", y = "Total Number",
         title = "a) Number of patches") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  ) +
  guides(colour=guide_legend(title="Climate Scenarios"))


smallest_area_plot <- ggplot(data = summary_table_scenarios, 
       aes(x = buffer, y = log(min_area_m2), 
           group = scenario, 
           colour = scenario 
       )) +
  geom_line() +
  theme_bw()+
  scale_color_manual(values=c("#00AA00", "#FFA500", "#FF4000")) +
  labs(x = "Buffer (m)", y = "Area log(m2)",
       title ="b) Smallest patch") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  ) +
  guides(colour=guide_legend(title="Climate Scenarios"))


largest_area_plot <- ggplot(data = summary_table_scenarios, 
       aes(x = buffer, y = log(max_area_m2), 
           group = scenario, 
           colour = scenario 
       )) +
  geom_line() +
  theme_bw()+
  scale_color_manual(values=c("#00AA00", "#FFA500", "#FF4000")) +
  labs(x = "Buffer (m)", y = "Area log(m2)",
       title ="c) Largest patch") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  ) +
  guides(colour = guide_legend(title="Climate Scenarios"))

# Comaparison plot with main features of the created graphs
ggarrange(total_patches_plot, smallest_area_plot, largest_area_plot, 
          align = "v", ncol = 1, 
          common.legend = T,
          legend = "right"
          )
