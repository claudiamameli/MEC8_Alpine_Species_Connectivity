# Library -----------------------------------------------------------------
# Importing functions and full library
source("functions_library.R")

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
# #### >>> Save files --------------------------------------------------------------
# 
# node_df_list <- list(g_supinum_pres = g_supinum_pres,
#                        g_supinum_fut245 = g_supinum_fut245,
#                        g_supinum_fut585 = g_supinum_fut585)
# 
# 
# for(i in names(node_df_list)){
#   write.csv(node_df_list[[i]],
#             paste0("node_df/", i,".csv"),
#             row.names = FALSE)
# }


# 3. Node visualisation --------------------------------
# Load all previously saved dataframes
nodes_df <- list.files("node_df", pattern = "\\.csv$", full.names = TRUE)
df_names <- tools::file_path_sans_ext(basename(nodes_df))

list2env(setNames(lapply(nodes_df, read.csv),
                  df_names),
         envir = .GlobalEnv)


plot_nodes(g_supinum_pres, "Present", "3")
plot_nodes(g_supinum_fut245, "Future 2-4.5", "3")
plot_nodes(g_supinum_fut585, "Future 5-8.5", "3")
