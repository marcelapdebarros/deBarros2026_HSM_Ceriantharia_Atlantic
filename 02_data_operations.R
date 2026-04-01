# Packages
# INFO ------------------------------------------------------------------------#
#                                                                              #
# Script name: 02 [HSM] - Data operations                                      #
#                                                                              #
# Author:      Marcela Aparecida de Barros                                     #
# Copyright    Copyright 2026 DE BARROS, M. A.                                 #
# Email:       marcela.barros@unesp.br                                         #
#                                                                              #
# Date:        01-04-2026                                                      #
#                                                                              #
#Here we are going to:                                                         #
# - Thin the occurrences and define the final dataset;                         #                               
# - VIF analysis and define the final environmental layers to each species;    #
# - Buffer delimitation (100 km);                                              #
# - Calibration data with ENMeval (Muscarella et al. 2014), which:             #
#   - < 30 occurrences: Jackknife + unique feature classes (L, Q, H, P, T)     #
#   - > 30 ocurrencs:  Block + unique and combined (LQ, LQH, LQHP, LQHPT)      #
#                                                                              #
# All analyses follow best practices for scientific reproducibility.           #
#                                                                              #
#------------------------------------------------------------------------------#

# 00 Packages and reproducibility ----------------------------------------------
library(ENMeval)
library(terra)
library(tidyverse)
library(sf)

set.seed(5315)

# 01 Load and thin species occurrence ------------------------------------------
occ <- read.csv("01_occurrences/complete/species.csv")

occ_thin <- spThin::thin(
  occ,
  lat.col = "latitude",
  long.col = "longitude",
  spec.col = "species",
  thin.par = 5.5,
  reps = 100,
  locs.thinned.list.return = TRUE,
  write.files = FALSE,
  write.log.file = FALSE
)

occ_thinned <- occ_thin[[1]]

mapview::mapview(
  occ_thinned,
  xcol = "Longitude",
  ycol = "Latitude",
  crs = 4326,
  grid = FALSE
)

occ <- occ_thinned

# 02 Load environmental data ---------------------------------------------------
env <- rast(
  list.files(
    path = "c:02_vars/wgs84/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# 02.1 Environmental buffer (100 km)
occ_sf <- st_as_sf(occ, coords = c("Longitude", "Latitude"), crs = 4326)
occ_buffer <- sf::st_buffer(occ_sf, dist = 100000) |> 
  sf::st_union() |> 
  sf::st_sf() |>
  sf::st_transform(crs = crs(env))
env_buffered <- crop(env, occ_buffer)
env_buffered <- mask(env_buffered, occ_buffer)
plot(env_buffered)

# 03 Defining final occurrences dataset (i. e., with env values)
occ_values <- terra::extract(env_buffered, occ, ID = FALSE)
head(occ_values)
anyNA(occ_values)

# 03.1 Check which points are NA
occ_na <- which(rowSums(is.na(occ_values)) > 0)
occ_valid <- occ[-occ_na, ]
names(occ_valid)[names(occ_valid) == "Longitude"] <- "longitude"
names(occ_valid)[names(occ_valid) == "Latitude"]  <- "latitude"

# Save the final occurrence dataset
write.csv(occ_valid, "final/species.csv", row.names = FALSE)

# 04 VIF Analysis (colinearity problems) ---------------------------------------
# 04.1 Convert SpatRaster to dataframe (sample points)
env_values <- terra::extract(env_buffered, occ_valid)
names(env)

# 04.2 View the data
summary(env_values)
str(env_values)

# 04.3 Calculate VIF
vif_calc <- usdm::vifstep(env_values, th = 5)
print(vif_calc)

# 04.4 Select uncorrelated variables
env_vif <- usdm::exclude(env_buffered, vif_calc)
summary(env_vif)

# 04.5 Export cropped rasters
for (i in 1:nlyr(env_vif)) {
  writeRaster(
    env_vif[[i]],
    filename = paste0("02_vars/vif/species/", names(env_vif)[i], ".tif"),
    overwrite = TRUE
  )
}

# 05 Calibration data with ENMEVAL ---------------------------------------------
# Create the background data, size defined by Tong et al. (2023) on 1:10
bg_data <- spatSample(
  env_vif,
  size = 24800, 
  na.rm = TRUE, 
  values = FALSE,
  xy = TRUE)  %>% 
  as.data.frame() %>%
  rename(
    longitude = x,
    latitude = y
  )

# 05.1 Partition data 1: Jackknife (< 30 occurrences)
jackk <- get.jackknife(occ_valid, bg_data)
evalplot.grps(pts = occ_valid, pts.grp = jackk$occs.grp, envs = env_vif)
evalplot.grps(pts = bg_data, pts.grp = jackk$bg.grp, envs = env_vif)
table(jackk$occs.grp)

# 05.2 Partition data 2: Block (> 30 occurrences)
blk <- get.block(occ_valid, bg_data, orientation = "lat_lon")
evalplot.grps(pts = occ_valid, pts.grp = blk$occs.grp, envs = env_vif)
evalplot.grps(pts = bg_data, pts.grp = blk$bg.grp, envs = env_vif)
table(blk$occs.grp)

# Evaluating the data's relationship (change the partitions and FCs)
test <- ENMevaluate(
  occs = occ_valid,
  envs = env_vif,
  bg = bg_data,
  algorithm = "maxent.jar",
  partitions = "block",
  tune.args = list(
    fc = c("L", "Q", "H", "P", "T", "LQ", "LQH", "LQHP", "LQHPT"),
    rm = seq(1, 5, by = 0.5)
  )
)

# To include: "LQ", "LQH", "LQHP", "LQHPT"

# Exploratory analysis to define the optimal configuration to modelling
model_results <- eval.results(test)
model_results %>%
  DT::datatable(
    options = list(
      pageLength = 100,
      scrollX = TRUE
    ),
    filter = "top"
  )

evalplot.stats(
  e = maxent_eval,
  stats = "or.10p",
  color = "fc",
  x.var = "rm"
)

eval_results <- eval.results(maxent_eval)

opt_aicc <- eval_results |>
  filter(delta.AICc < 1)
opt_aicc

model_selection <- eval.models(maxent_eval)[[opt_aicc$tune.args]]
model_selection$betas
plot(model_selection, type = "cloglog")
dev.off()

prediction_preview <- eval.predictions(maxent_eval)[[as.character(opt_aicc$tune.args)]]
plot(prediction_preview, col = terrain.colors(100))