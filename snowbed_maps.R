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


# 1. DATA IMPORT & CLEANUP---------------------------------------------------------
## 1.1 Species Data  ------------------------------------------------------------
df_species_tot <- read.csv("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/species_wide.csv")
str(df_species_tot)

# Gnaphalium supinum, Luzula alpino-pilosa, Salix herbacea, Sibbaldia procumbens
colnames(df_species_tot)[grepl("Gnaphalium", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Luzula", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Salix", colnames(df_species_tot))]
colnames(df_species_tot)[grepl("Sibbaldia", colnames(df_species_tot))]

df_species_bin <- df_species_tot
df_species_bin[ , -1] <- as.factor(ifelse(df_species_bin[ , -1] > 0, 1, 0))

target_species <- df_species_bin %>% 
  select("logger_ID", "Gnaphalium.supinum", "Luzula.alpino.pilosa", "Salix.herbacea", "Sibbaldia.procumbens")



# # Uncomment if need to save the file again
# write.csv(target_species, "~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/target_species.csv")

## 1.2 Geo Data  ----------------------------------------------------------------
st_layers("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/site_data.gpkg")
geo_data_full <- read_sf("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/site_data.gpkg")

plot(geo_data$geom)

all_data <- target_species %>% 
  left_join(select(geo_data_full, logger_ID, geom), by = "logger_ID")


g_supinum_data <- all_data %>%
  mutate(pres = Gnaphalium.supinum) %>% 
  select("logger_ID", "pres", "geom")

l_alpinospinosa_data <- all_data %>%
  mutate(pres = Luzula.alpino.pilosa) %>% 
  select("logger_ID", "pres", "geom")

s_herbacea_data <- all_data %>%
  mutate(pres = Salix.herbacea) %>% 
  select("logger_ID", "pres" , "geom")

s_procumbens_data <- all_data %>%
  mutate(pres = Sibbaldia.procumbens) %>% 
  select("logger_ID", "pres", "geom")
## 1.3 Snowcover Data -------------------------------------------------------------
snowbeds_tif <- list.files("~/Desktop/GitHub/MEC8_Snowbed_Alpine_Species/Data/snow_cover_maps" ,pattern = "tif$", full.names=T)

snowbed_list <-  rast(snowbeds_tif)
plot(snowbed_list)
names(snowbed_list)

years <- stringr::str_extract(names(snowbed_list), "\\d{4}") %>% as.numeric()
months <- stringr::str_extract(names(snowbed_list), "\\d_\\d_")

snowbed_pres_1_7 <- snowbed_list[[years >= 2015 & years <= 2025 & months == "1_7_"]]
snowbed_fut_1_7  <- snowbed_list[[years >= 2090 & years <= 2100 & months == "1_7_"]]
plot(snowbed_pres_1_7)
plot(snowbed_fut_1_7) 

snowbed_pres_1_7_sum <- app(snowbed_pres_1_7, mean)
snowbed_fut_1_7_sum <- app(snowbed_fut_1_7, mean)
plot(snowbed_pres_1_7_sum)
plot(snowbed_fut_1_7_sum) 


snowbed_pres_1_5 <- snowbed_list[[years >= 2015 & years <= 2025 & months == "1_5_"]]
snowbed_fut_1_5  <- snowbed_list[[years >= 2090 & years <= 2100 & months == "1_5_"]]
plot(snowbed_pres_1_5)
plot(snowbed_fut_1_5) 

snowbed_pres_1_5_sum <- app(snowbed_pres_1_5, mean)
snowbed_fut_1_5_sum <- app(snowbed_fut_1_5, mean)
plot(snowbed_pres_1_5_sum)
plot(snowbed_fut_1_5_sum) 


plot(snowbed_pres_1_7_sum)
points(g_supinum_data$geom[g_supinum_data$pres == 1], col="red", pch = 16, cex=1)
points(l_alpinospinosa_data$geom[l_alpinospinosa_data$pres == 1], col="purple", pch = 16, cex=1)
points(s_herbacea_data$geom[s_herbacea_data$pres == 1], col="blue", pch = 16, cex=1)
points(s_procumbens_data$geom[s_procumbens_data$pres ==1], col="green", pch = 16, cex=1)


pres_1_7_mean <- terra::extract(snowbed_pres_1_7_sum, geo_data, ID = FALSE)
names(pres_1_7_mean) <- "pres_1_7"
fut_1_7_mean <- terra::extract(snowbed_fut_1_7_sum, geo_data, ID = FALSE)
names(fut_1_7_mean) <- "fut_1_7"
pres_1_5_mean <- terra::extract(snowbed_pres_1_5_sum, geo_data, ID = FALSE)
names(pres_1_5_mean) <- "pres_1_5"
fut_1_5_mean <- terra::extract(snowbed_fut_1_5_sum, geo_data, ID = FALSE)
names(fut_1_5_mean) <- "fut_1_5"

g_supinum_data <- data.frame(g_supinum_data, pres_1_7_mean, fut_1_7_mean, pres_1_5_mean, fut_1_5_mean)
l_alpinospinosa_data <- data.frame(l_alpinospinosa_data, pres_1_7_mean, fut_1_7_mean, pres_1_5_mean, fut_1_5_mean)
s_herbacea_data <- data.frame(s_herbacea_data, pres_1_7_mean, fut_1_7_mean, pres_1_5_mean, fut_1_5_mean)
s_procumbens_data <- data.frame(s_procumbens_data, pres_1_7_mean, fut_1_7_mean, pres_1_5_mean, fut_1_5_mean)

# 3. VISUALISATION & EXTRAPOLATION-----------------------------------------
## 3.1 Present January - July Boxplots -------------------------------------------------
ggplot(g_supinum_data, aes(y = pres_1_7, x = pres))+
  geom_boxplot()+
  labs(title = "Gnaphalium.supinum - Present_1_7", 
       x = "Presence/Absence", 
       y = "Days of Snowcover")+
  theme_bw()

ggplot(l_alpinospinosa_data, aes(y = pres_1_7, x = pres))+
  geom_boxplot()+
  labs(title = "Luzula.alpino.pilosa - Present_1_7",
       x = "Presence/Absence", 
       y = "Days of Snowcover")+
  theme_bw()

ggplot(s_herbacea_data, aes(y = pres_1_7, x = pres))+
  geom_boxplot()+
  labs(title = "Salix.herbacea - Present_1_7",
       x = "Presence/Absence", 
       y = "Days of Snowcover")+
  theme_bw()

ggplot(s_procumbens_data, aes(y = pres_1_7, x = pres))+
  geom_boxplot()+
  labs(title = "Sibbaldia.procumbens - Present_1_7",
       x = "Presence/Absence", 
       y = "Days of Snowcover")+
  theme_bw()



## 3.2 Summary of Snowcover days - Jan-Jul -----------------------------------------
### > Gnaphalium.supinum --------------------------------------------------------------
summary_g_supinum_1_7 <- g_supinum_data %>%
  filter(pres == 1) %>% 
  summarise(
    species = "g_supinum",
    q25 = round(quantile(pres_1_7, 0.25, na.rm = TRUE)),
    mean = round(mean(pres_1_7, na.rm = TRUE)),
    sd = round(sd(pres_1_7, na.rm = TRUE)),
    median = round(median(pres_1_7, na.rm = TRUE)),
    q75= round(quantile(pres_1_7, 0.75, na.rm = TRUE)),
    q95 = round(quantile(pres_1_7, 0.95, na.rm = TRUE))
  )

### > Luzula.alpino.pilosa -------------------------------------------------------------
summary_l_alpinospinosa_pres_1_7 <- l_alpinospinosa_data %>%
  filter(pres == 1) %>% 
  summarise(
    species = "l_alpinospinosa",
    q25 = round(quantile(pres_1_7, 0.25, na.rm = TRUE)),
    mean = round(mean(pres_1_7, na.rm = TRUE)),
    sd = round(sd(pres_1_7, na.rm = TRUE)),
    median = round(median(pres_1_7, na.rm = TRUE)),
    q75= round(quantile(pres_1_7, 0.75, na.rm = TRUE)),
    q95 = round(quantile(pres_1_7, 0.95, na.rm = TRUE))
  )

### > Salix.herbacea -------------------------------------------------------------
summary_s_herbacea_pres_1_7 <- s_herbacea_data %>%
  filter(pres == 1) %>% 
  summarise(
    species = "s_herbacea",
    q25 = round(quantile(pres_1_7, 0.25, na.rm = TRUE)),
    mean = round(mean(pres_1_7, na.rm = TRUE)),
    sd = round(sd(pres_1_7, na.rm = TRUE)),
    median = round(median(pres_1_7, na.rm = TRUE)),
    q75= round(quantile(pres_1_7, 0.75, na.rm = TRUE)),
    q95 = round(quantile(pres_1_7, 0.95, na.rm = TRUE))
  )

### > Sibbaldia.procumbens -------------------------------------------------------------
summary_s_procumbens_1_7 <- s_procumbens_data %>%
  filter(pres == 1) %>% 
  summarise(
    species = "s_procumbens",
    q25 = round(quantile(pres_1_7, 0.25, na.rm = TRUE)),
    mean = round(mean(pres_1_7, na.rm = TRUE)),
    sd = round(sd(pres_1_7, na.rm = TRUE)),
    median = round(median(pres_1_7, na.rm = TRUE)),
    q75= round(quantile(pres_1_7, 0.75, na.rm = TRUE)),
    q95 = round(quantile(pres_1_7, 0.95, na.rm = TRUE))
  )


### 3.2.1 Total summary of Snowcover days -----------------------------------------
summary_snowcover_1_7 <- summary_g_supinum_1_7 %>% 
  rbind(summary_l_alpinospinosa_pres_1_7, summary_s_herbacea_pres_1_7, summary_s_procumbens_1_7)

summary_snowcover_1_7



# D. BINARY MAPS TRIALS  --------------------------------------------------------

trials <- snowbed_pres_1_7_sum
trials_fut <- snowbed_fut_1_7_sum

g_supinum_mean <- (summary_snowcover_1_7$mean[summary_snowcover_1_7$species == "g_supinum"])

tol <- 0.5

trials <- ifel(
  abs(trials - g_supinum_mean) <= tol,
  1,
  0)

trials_fut <- ifel(
  abs(trials_fut - g_supinum_mean) <= tol,
  1,
  0)


plot(trials)
plot(trials_fut)

