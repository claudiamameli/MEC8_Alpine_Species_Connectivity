# Library -----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(gpkg)
library(terra)
library(sf)

# Specied Data Import & Cleanup ------------------------------------------------------------
df_species_tot <- read.csv("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/species_wide.csv")
str(df_species_tot)

# Gnaphalium supinum, Luzula alpino-pilosa, Salix herbacea, Sibbaldia procumbens
colnames(df_species_tot)[grepl("Gnaphalium", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Luzula", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Salix", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Sibbaldia", colnames(df_species_tot))]

df_species <- df_species_tot %>% 
  select("logger_ID", "Gnaphalium.supinum", "Luzula.alpino.pilosa", "Salix.herbacea", "Sibbaldia.procumbens")

write.csv(df_species, "~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/df_species.csv")



# Geo Data  Import & Cleanup ----------------------------------------------------------------
st_layers("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/site_data.gpkg")
geo_data <- read_sf("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/site_data.gpkg")

plot(geo_data$geom)


all_data <- df_species %>% 
  cbind(geo_data[-1])


g_supinum_data <- df_species %>%
  select("logger_ID", "Gnaphalium.supinum") %>% 
  filter(Gnaphalium.supinum != 0) %>% 
  left_join(geo_data, by = "logger_ID")

l_alpinospinosa_data <- df_species %>%
  select("logger_ID", "Luzula.alpino.pilosa") %>% 
  filter(Luzula.alpino.pilosa != 0) %>% 
  left_join(geo_data, by = "logger_ID")

s_herbacea_data <- df_species %>%
  select("logger_ID", "Salix.herbacea") %>% 
  filter(Salix.herbacea != 0) %>% 
  left_join(geo_data, by = "logger_ID")

s_procumbens_data <- df_species %>%
  select("logger_ID", "Sibbaldia.procumbens") %>% 
  filter(Sibbaldia.procumbens != 0) %>% 
  left_join(geo_data, by = "logger_ID")
  
