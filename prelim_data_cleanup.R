# Library -----------------------------------------------------------------
# Importing functions and full library
source("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/functions_library.R")

# 1. Species Data  ------------------------------------------------------------
## 1.1 Presence/absence ----------------------------------------------------
df_species_tot <- read.csv("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/species_wide.csv")
str(df_species_tot)

# Gnaphalium supinum, Luzula alpino-pilosa, Salix herbacea, Sibbaldia procumbens
colnames(df_species_tot)[grepl("Gnaphalium", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Luzula", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Salix", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Sibbaldia", colnames(df_species_tot))]

# Change from abundance into pres/abs data
df_species_bin <- df_species_tot
df_species_bin[ , -1] <- as.numeric(ifelse(df_species_bin[ , -1] > 0, 1, 0))

# Select target species
target_species <- df_species_bin %>% 
  dplyr::select("logger_ID", "Gnaphalium.supinum", "Luzula.alpino.pilosa", "Salix.herbacea", "Sibbaldia.procumbens")



# # Uncomment if need to save the file again
# write.csv(target_species, "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/target_species.csv")

## 1.2 Geo Data - Schrankogel ----------------------------------------------------------------
geo_data_full <- read_sf("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/site_data.gpkg")

target_species_data <- target_species %>% 
  left_join(dplyr::select(geo_data_full, logger_ID, geom), by = "logger_ID")


## 1.3 Kernel data ---------------------------------------------------------
df_species_kernels <- read.csv("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/kernels_clean.csv")
head(df_species_kernels)
str(df_species_kernels)

# Gnaphalium supinum, Luzula alpino-pilosa, Salix herbacea, Sibbaldia procumbens
df_species_kernels$species[grep("gnaphalium", df_species_kernels$species)]
df_species_kernels$species[grep("luzula", df_species_kernels$species)]
df_species_kernels$species[grep("salix", df_species_kernels$species)]
df_species_kernels$species[grep("sibbaldia", df_species_kernels$species)]

# Select target species
target_species_names <- c("gnaphalium_supinum", "luzula_alpino-pilosa", "salix_herbacea", "sibbaldia_procumbens")

kernels_target_sp <- df_species_kernels %>%
  filter(species %in% target_species_names)

kernels_target_sp

# 2. Climate data  --------------------------------------------------------
## 2.1 Snowcover data (Jan - Jul)------------------------------------------------------
snowbeds_list <- list.files("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/snow_cover_maps" ,pattern = "tif$", full.names=T)

snowbed_rast <-  rast(snowbeds_list)
#plot(snowbed_rast)
names(snowbed_rast)

years <- stringr::str_extract(names(snowbed_rast), "\\d{4}") %>% as.numeric()
months <- stringr::str_extract(names(snowbed_rast), "\\d_\\d_")


### a. Snowcover - Present -----------------------------------------------------
snowbed_pres_1_7 <- snowbed_rast[[years >= 2015 & years <= 2025 & months == "1_7_"]]
#plot(snowbed_pres_1_7)

snowbed_pres_1_7_mean <- app(snowbed_pres_1_7, mean)
plot(snowbed_pres_1_7_mean)

### b. Snowcover - Future ------------------------------------------------------
snowbed_fut_1_7  <- snowbed_rast[[years >= 2090 & years <= 2100 & months == "1_7_"]]
#plot(snowbed_fut_1_7) 

snowbed_fut_1_7_mean <- app(snowbed_fut_1_7, mean)
plot(snowbed_fut_1_7_mean) 




### /// Example visualisation - Snowcover v Presence ----------------------------
# Uncomment if needed 

# plot(snowbed_pres_1_7_mean)
# points(g_supinum_data$geom[g_supinum_data$pres == 1], col="red", pch = 16, cex=1)
# points(l_alpinospinosa_data$geom[l_alpinospinosa_data$pres == 1], col="purple", pch = 16, cex=1)
# points(s_herbacea_data$geom[s_herbacea_data$pres == 1], col="blue", pch = 16, cex=1)
# points(s_procumbens_data$geom[s_procumbens_data$pres ==1], col="green", pch = 16, cex=1)







## 2.2 Temperature ---------------------------------------------------------
### a. Temperature - Present ---------------------------------------------------
pres_temp_list <- list.files("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/temp_pres" ,pattern = "tif$", full.names=T)

pres_temp_list

pres_temp_rast <-  rast(pres_temp_list)
#plot(pres_temp_rast)
names(pres_temp_rast)

years <- stringr::str_extract(names(pres_temp_rast), "\\d{4}") %>% as.numeric()

temp_pres <- pres_temp_rast[[years >= 2015 & years <= 2023]]

#plot(temp_pres)

temp_pres_mean <- app(temp_pres, mean)
plot(temp_pres_mean)
temp_pres_mean_crop <- crop(temp_pres_mean, snowbed_pres_1_7_mean)
plot(temp_pres_mean_crop)

### b. Temperature - Future_245 ------------------------------------------------
fut245_temp_list <- list.files("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/temp_fut_245" ,pattern = "tif$", full.names=T)

fut245_temp_list

fut245_temp_rast <-  rast(fut245_temp_list)

#plot(fut245_temp_rast)
names(fut245_temp_rast)

years <- stringr::str_extract(names(fut245_temp_rast), "\\d{4}") %>% as.numeric()

temp_fut245 <- fut245_temp_rast[[years >= 2090 & years <= 2100]]

#plot(temp_fut245)

temp_fut245_mean <- app(temp_fut245, mean)
#plot(temp_fut245_mean)

temp_fut245_mean_crop <- crop(temp_fut245_mean, snowbed_pres_1_7_mean)
plot(temp_fut245_mean_crop)

### c. Temperature - Future_585 ------------------------------------------------
fut585_temp_list <- list.files("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/temp_fut_585/" ,pattern = "tif$", full.names=T)

fut585_temp_list

fut585_temp_rast <-  rast(fut585_temp_list)

#plot(fut585_temp_rast)
names(fut585_temp_rast)

years <- stringr::str_extract(names(fut585_temp_rast), "\\d{4}") %>% as.numeric()

temp_fut585 <- fut585_temp_rast[[years >= 2090 & years <= 2100]]

#plot(temp_fut585)

temp_fut585_mean <- app(temp_fut585, mean)
#plot(temp_fut585_mean)

temp_fut585_mean_crop <- crop(temp_fut585_mean, snowbed_pres_1_7_mean)
plot(temp_fut585_mean_crop)



## 2.3 Carbon --------------------------------------------------------------
carbon <- rast("~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/carbon/carbon.tif")
plot(carbon)


# 3. Map Values Extractions ---------------------------------------
## 3.1 Snowcover ---------------------------------------------------------------
snow_pres_1_7_extract <- terra::extract(snowbed_pres_1_7_mean, geo_data_full, ID = FALSE)
names(snow_pres_1_7_extract) <- "snow_pres_1_7"

snow_pres_1_7_2021_extract <- terra::extract(snowbed_pres_1_7$snowcover1_7_2021, geo_data_full, ID = FALSE)
names(snow_pres_1_7_2021_extract) <- "snow_pres_1_7_2021"

## 3.2 Temperature  ------------------------------------------------------------
temp_pres_extract <- terra::extract(temp_pres_mean_crop, geo_data_full, ID = FALSE)
names(temp_pres_extract) <- "temp_pres"

temp_pres_2021_extract <- terra::extract(temp_pres$MSummer_2021, geo_data_full, ID = FALSE)
names(temp_pres_2021_extract) <- "temp_pres_2021"

## 3.3 Carbon ------------------------------------------------------------------
carbon_extract <- terra::extract(carbon, geo_data_full, ID = FALSE)
names(carbon_extract) <- "carbon"

# 4. Dataframe constructions ----------------------------------------------------
## 4.1 Final table for all species and all climate data -------------------------------------------------------------
full_data_df <- data.frame(target_species_data, 
                           snow_pres_1_7_extract,
                           snow_pres_1_7_2021_extract,
                           temp_pres_extract,
                           temp_pres_2021_extract,
                           carbon_extract)
str(full_data_df)


## 4.2 Summary tables ----------------------------------------------------------
### a. g_supinum ---------------------------------------------------------------
summary_g_supinum <- full_data_df %>%
  filter(Gnaphalium.supinum == 1) %>%
  pivot_longer(
    cols = c(
      snow_pres_1_7,
      temp_pres,
      carbon
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  summarise(
    q25    = round(quantile(value, 0.25, na.rm = TRUE)),
    mean   = round(mean(value, na.rm = TRUE)),
    sd     = round(sd(value, na.rm = TRUE)),
    median = round(median(value, na.rm = TRUE)),
    q75    = round(quantile(value, 0.75, na.rm = TRUE)),
    q95    = round(quantile(value, 0.95, na.rm = TRUE)),
    .by = metric
  )


### b. l_alpino_pilosa ---------------------------------------------------------
summary_l_alpino_pilosa <- full_data_df %>%
  filter(Luzula.alpino.pilosa == 1) %>%
  pivot_longer(
    cols = c(
      snow_pres_1_7,
      temp_pres,
      carbon
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  summarise(
    q25    = round(quantile(value, 0.25, na.rm = TRUE)),
    mean   = round(mean(value, na.rm = TRUE)),
    sd     = round(sd(value, na.rm = TRUE)),
    median = round(median(value, na.rm = TRUE)),
    q75    = round(quantile(value, 0.75, na.rm = TRUE)),
    q95    = round(quantile(value, 0.95, na.rm = TRUE)),
    .by = metric
  )


### c. s_herbacea --------------------------------------------------------------
summary_s_herbacea <- full_data_df %>%
  filter(Salix.herbacea == 1) %>%
  pivot_longer(
    cols = c(
      snow_pres_1_7,
      temp_pres,
      carbon
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  summarise(
    q25    = round(quantile(value, 0.25, na.rm = TRUE)),
    mean   = round(mean(value, na.rm = TRUE)),
    sd     = round(sd(value, na.rm = TRUE)),
    median = round(median(value, na.rm = TRUE)),
    q75    = round(quantile(value, 0.75, na.rm = TRUE)),
    q95    = round(quantile(value, 0.95, na.rm = TRUE)),
    .by = metric
  )


### d. s_procumbens ------------------------------------------------------------
summary_s_procumbens <- full_data_df %>%
  filter(Sibbaldia.procumbens == 1) %>%
  pivot_longer(
    cols = c(
      snow_pres_1_7,
      temp_pres,
      carbon
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  summarise(
    q25    = round(quantile(value, 0.25, na.rm = TRUE)),
    mean   = round(mean(value, na.rm = TRUE)),
    sd     = round(sd(value, na.rm = TRUE)),
    median = round(median(value, na.rm = TRUE)),
    q75    = round(quantile(value, 0.75, na.rm = TRUE)),
    q95    = round(quantile(value, 0.95, na.rm = TRUE)),
    .by = metric
  )

# 5. FINAL DATAFRAMES ---------------------------------------------------------
full_data_df

summary_g_supinum
summary_l_alpino_pilosa
summary_s_herbacea
summary_s_procumbens



df_list <- list(full_data_df = full_data_df,
                summary_g_supinum = summary_g_supinum, 
                summary_l_alpino_pilosa = summary_l_alpino_pilosa, 
                summary_s_herbacea = summary_s_herbacea, 
                summary_s_procumbens = summary_s_procumbens)

path <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/Data/prelim_dfs/"


for(i in names(df_list)){
  write.csv(df_list[[i]], paste0(path, i,".csv"))
}



# 6. FINAL MEAN MAPS ------------------------------------------------------
plot(snowbed_pres_1_7_mean)
plot(snowbed_fut_1_7_mean)
plot(temp_pres_mean_crop)
plot(temp_fut245_mean_crop)
plot(temp_fut585_mean_crop)

tifs_list <- list(snowbed_pres_1_7_mean = snowbed_pres_1_7_mean,
                  snowbed_fut_1_7_mean = snowbed_fut_1_7_mean,
                  temp_pres_mean_crop = temp_pres_mean_crop,
                  temp_fut245_mean_crop = temp_fut245_mean_crop,
                  temp_fut585_mean_crop = temp_fut585_mean_crop)

path_tifs <- "~/Desktop/Repositories/MEC8_Snowbed_Alpine_Species/gis_map_overlay/"

#for(i in names(tifs_list)){
#  writeRaster(tifs_list[[i]], paste0(path_tifs, i,".tif"))
#}



# Regression line  --------------------------------------------------------
full_data_df_clean <- na.omit(full_data_df)

gnap_glm <- glm(data = full_data_df_clean, 
                Gnaphalium.supinum ~ snow_pres_1_7 + temp_pres + carbon, 
                family = binomial)
summary(gnap_glm)

sum_gnap_glm <- summary(gnap_glm)
d_sq_glm <- ((sum_gnap_glm$null.deviance - sum_gnap_glm$deviance)/sum_gnap_glm$null.deviance)

# also tried gams and poly but they all almost give the same results - the regression does not show significance with these parameters, so does it make sense using them? I am a little confused personally. 

z_values <- c(snow = 0.912, temp = 0.080, carbon = 0.631)
percentages <- round((abs(z_values) / sum(abs(z_values))) * 100, 1)
print(percentages)



gnap_glm_poly <- glm(data = full_data_df_clean, 
                Gnaphalium.supinum ~ poly(snow_pres_1_7,2) + poly(temp_pres,2) + carbon, 
                family = binomial)
summary(gnap_glm_poly)

gnap_gam_s <- gam(data = full_data_df_clean, 
                Gnaphalium.supinum ~ s(snow_pres_1_7) + s(temp_pres) + carbon, 
                family = binomial)
summary(gnap_gam_s)

gnap_gam_s4 <- gam(data = full_data_df_clean, 
                Gnaphalium.supinum ~ s(snow_pres_1_7,4) + s(temp_pres,4) + carbon, 
                family = binomial)
summary(gnap_gam_s4)
