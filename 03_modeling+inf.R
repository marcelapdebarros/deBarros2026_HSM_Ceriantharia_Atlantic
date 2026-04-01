# INFO BOX --- ----------------------------------------------------------------#
#                                                                              #
# Script name: 03 [HSM] - Modeling+infinite                                    #
#                                                                              #
# Author:      Marcela Aparecida de Barros                                     #
# Copyright    Copyright 2026 DE BARROS, M. A.                                 #
# Email:       marcela.barros@unesp.br                                         #
#                                                                              #
# Date:        01-04-2026                                                      #
#                                                                              #
# Description: This script is for the modeling process. After loading the      #
# occurrences and environmental variables, the models will be calibrated,      #
# trained, and projected, followed by the creation of the final ensemble.      #
#                                                                              #
# All analyses follow best practices for scientific reproducibility.           #
#------------------------------------------------------------------------------#

# 00 Packages and reproducibility ----------------------------------------------
library(dplyr)
library(terra)
library(tidyverse)
library(sdm)

set.seed(5315)

# 01 Load datasets --------------------------------------------------------------
occ <- read.csv("01_occurrences/final/species.csv")

env <- rast(
  list.files(
    path = "02_vars/vif/species/",
    pattern = ".tif$",
    full.names = TRUE
  )
  
)

# 02 Modeling core  -----------------------------------------------------------------
occ$species <- 1

# Create the background data
bg_data <- background(
  env,
  n = 880,
  method = "gRandom"
)

model_bg_data <- sdmData(
  formula = species ~. + coords(longitude + latitude),
  train = occ,
  predictors = env,
  bg = bg_data
)

model_bg_data

# Train the model
maxent_model <- sdm(
  formula = multiplispeciescatus ~ . + coords(longitude + latitude),
  data = model_bg_data,
  methods = "maxent",
  replication = "cv",
  cv.folds = 5,
  n = 10,
  test.percent = 30,
  modelSettings = list(
    maxent = list(
      beta = 2,
      feat = c("linear", "quadratic", "hinge")
    )
  )
)

maxent_model

# Variables relative importance
getVarImp(maxent_model)

write.sdm(
  maxent_model,
  "03_models/species"
)  

# 03 Projection ---------------------------------------------------------------
# Load, crop and mask the data
occ <- read.csv("01_occurrences/final/species.csv")
maxent_model <- read.sdm("03_models/species.sdm")

env <- rast(
  list.files(
    path = "02_vars/projs/",
    pattern = ".tif$",
    full.names = TRUE
  )
)

# Study zone (Atlantic Ocean)
atlantic <- terra::vect("Shapefiles/atlantic_ocean/iho.shp")
env_atl <- terra::project(atlantic, env)
env_atl <- terra::crop(env, atlantic)
env_atl <- terra::mask(env_atl, atlantic)

plot(env_atl)

# Projection
env_bio <- subset(env_atl, maxent_model@data@features.name)
plot(env_bio)

# Projection
maxent_proj <- predict(
  maxent_model,
  newdata = env_bio,
  filename = "04_projs/species.tif",
  overwrite = TRUE,
  parallelSetting = list(
    ncore = 10
  )
)

plot(maxent_proj)

# Ensemble
maxent_ens <- ensemble(
  maxent_model,
  newdata = maxent_proj,
  setting = list(
    method = "weighted",
    stat = "TSS",
    opt = 2
  )
)

plot(maxent_ens)

# Save
terra::writeRaster(
  maxent_ens,
  "05_ensemble/species.tif",
  overwrite = TRUE
)