# Library -----------------------------------------------------------------
# Importing functions and full library
source('functions_library.R')
set.seed(41032) # Setting seed for random calculation to return the same values 

# 1. Data preparation and functions ------------------------------------------

## 1.1 SPECIES DATA - function ---------------------------------------------------------------------
#' Function for formatting specific species data 
#' When used in a loop, different species can be implemented to calculate specific maps

get_species_data <- function(species) {
  spe_data <- read_csv("Data/species_data/species_long.csv") |>
    select(logger_ID, name, value) |>
    left_join(read_csv("Data/species_data/T5_logger_coordinates_EPSG25832.csv")) |>
    st_as_sf(coords = c("X", "Y"), crs = 25832) |>
    filter(name == species)
  return(spe_data)
}


## 1.2 PRESENT CONDITION EXTRACTION ---------------------------------------------------------------------

### a. Days of snowcover - present --------------------------------------------

tibble(files = list.files("Data/predictors/snow_past", full.names = T)) |>
  filter(grepl("1_7", files)) |>
  filter(grepl("2017|2018|2019|2020|2021", files)) |>
  pull(files) |> 
  rast() |> 
  app("mean") |>
  setNames('snow') -> snow_current


### b. Soil temperature - present ----------------------------------------------

tibble(files = list.files("Data/predictors/temp_past", full.names = T)) |>
  filter(grepl("2017|2018|2019|2020|2021", files)) |>
  pull(files) |> 
  rast() |> 
  app("mean") |>
  setNames('temperature') -> temperature_current


### c. Static soil and topographic predictors ----------------------------------
# Make sure the predictors are implemented into functions_library.csv for glm, gam and rf models.

other_variables <- rast(c("Data/predictors/carbon/carbon_crop.tif", 
                          "Data/predictors/topography/HSD_1.tif", 
                          "Data/predictors/topography/SLP_1.tif", 
                          "Data/predictors/topography/TPI_3.tif", 
                          "Data/predictors/topography/TWI_25.tif"))

## 1.2 FUTURE PREDICTORS EXTRACTION - functions --------------------------------------------
# These functions are used to extrapolate future conditions under different climate scenarios.

### a. Days of snowcover - future ----------------------------------------------

get_future_snow <- function(scenario) {
  tibble(files = list.files(paste0("Data/predictors/snow_future/", scenario), full.names = T)) |>
    filter(grepl("1_7", files)) |>
    filter(grepl("2095|2096|2097|2098|2099", files)) |>
    pull(files) |>
    rast() |>
    app("mean") |>
    setNames("snow") -> snow_future
  return(snow_future)
}


### b. Soil temperature - future --------------------------------------------

get_future_temperature <- function(scenario) {
  tibble(files = list.files(paste0("Data/predictors/temp_future/", scenario), full.names = T)) |>
    filter(grepl("2095|2096|2097|2098|2099", files)) |>
    pull(files) |>
    rast() |>
    app("mean") |>
    setNames("temperature") -> temperature_future
  return(temperature_future)
}


# 2. Models and predictions -----------------------------------------------
# Put together dataframes for model training with all chosen predictors
training_df <- c(other_variables, snow_current, temperature_current)
rs_df <- as.data.frame(training_df, xy = T)

# Uncomment to re-run
# Create predicted suitability binary rasters for all current and future scenario for chosen species
# for (target_species in c("Gnaphalium supinum"
#                          # "Luzula alpino-pilosa", 
#                          # "Salix herbacea", 
#                          # "Sibbaldia procumbens"
# )) {
#   print(target_species)  
#   clean_name <- tolower(gsub(" ", "_", target_species)) # for well formatted outputs
#   
#   spe_data <- get_species_data(target_species)
#   
#   df <- bind_cols(as_tibble(spe_data)["value"], as_tibble(terra::extract(training_df, spe_data))[-1])
#   
#   stats_out <- get_model_stats(df)
#   write.csv(stats_out, paste0("Data/stats/", clean_name, "_model_stats.csv"), row.names = FALSE)
#   print("Model statistics computed")
#   
#   ensemble_prediction_current <- get_ensemble_prediction(df, rs_df, stats_out)
#   print("Current ensamble prediction finished")
#   
#   ensemble_prediction_current |>
#     writeRaster(paste0("Data/map_predictions/", clean_name, "_current.tif"), overwrite = T)
#   
#   print("Starting future predictions")
#   for (scenario in list.files("Data/predictors/temp_future")) {
#     print(scenario)
#     rs_df_future <- as.data.frame(c(other_variables, 
#                                     get_future_snow(scenario),
#                                     get_future_temperature(scenario)), xy = T)
#     get_ensemble_prediction(df, rs_df_future, stats_out) |>
#       setNames(scenario) |>
#       writeRaster(paste0("Data/map_predictions/", clean_name, "_", scenario, ".tif"), overwrite = T)
#   }
# }



# 3. Visualisation of habitat gain/loss -----------------------------------
# Reload maps created 
files_rast <- list.files("Data/map_predictions" ,pattern = "tif$", full.names=T)
stats_out <- read_csv("Data/stats/gnaphalium_supinum_model_stats.csv") 

list2env(setNames(lapply(files_rast, rast),
                  tools::file_path_sans_ext(basename(files_rast))), envir = .GlobalEnv)


ssp245_change <- plot_changes(gnaphalium_supinum_current, gnaphalium_supinum_ssp245, "SSP2-4.5", "Gnaphalium supinum", "A) ")
ssp585_change <- plot_changes(gnaphalium_supinum_current, gnaphalium_supinum_ssp585, "SSP5-8.5", "Gnaphalium supinum", "B) ")

(ssp245_change | ssp585_change) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Gnaphalium supinum",
    subtitle = paste0("TSS: ", round(stats_out$TSS[stats_out$model == "ensemble"], 2))
  ) &
  theme(
    plot.title = element_text(face = "italic", size = 20),
    plot.subtitle = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
  )

