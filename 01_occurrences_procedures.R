# INFO BOX --------------------------------------------------------------------#
#                                                                              #
# Script name: 01 [HSM] - Gathering and processing occurrences                 #
#                                                                              #
# Author:      Marcela Aparecida de Barros                                     #
# Copyright    Copyright 2026 DE BARROS, M. A.                                 #
# Email:       marcela.barros@unesp.br                                         #
#                                                                              #
# Date:        01-04-2026                                                      #
#                                                                              #
# Description: In this script, you will download the occurrences from the GBIF #
# and OBIS datatases, save the raw data, clean the occurrences and prepare the #
# complete dataset to the thinning process later.                              #
#                                                                              #
# All analyses follow best practices for scientific reproducibility.           #
#------------------------------------------------------------------------------#

# 00 Packages and reproducibility ----------------------------------------------
library(robis)
library(rgbif)
library(dplyr)
library(readr)
library(tidyverse)

set.seed(5315)

# Species name
species <- "Genus species"

# 01 OBIS Occurrences ----------------------------------------------------------
# Raw data
obis_raw <- occurrence(species)

# Save
write_csv(obis_raw, "01_occurrences/solitarius/raw_obis.csv")

# Clean data
obis_clean <- obis_raw %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  select(
    latitude = decimalLatitude,
    longitude = decimalLongitude,
    species = scientificName
  ) %>%
  distinct(latitude, longitude, .keep_all = TRUE)


# Save
write_csv(obis_clean, "01_occurrences/americana/clean_obis.csv")

# 02 GBIF Occurrences ------------------------------------------------------------
# Search taxonKey
taxon_key <- name_backbone(name = species)$usageKey

# Raw data
gbif_raw <- occ_search(
  taxonKey = taxon_key,
  hasCoordinate = TRUE,
  limit = 200000
)$data

# Save
write_csv(gbif_raw, "01_occurrences/solitarius/raw_gbif.csv")

# Filtering and processing
gbif_clean <- gbif_raw %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  select(
    latitude = decimalLatitude,
    longitude = decimalLongitude,
    species
  ) %>%
  distinct(latitude, longitude, .keep_all = TRUE)

# Save
write_csv(gbif_clean, "01_occurrences/solitarius/clean_gbif.csv")

# Bind datasets
obis_gbif <- bind_rows(obis_clean, gbif_clean) %>%
  distinct(latitude, longitude, .keep_all = TRUE)

# Save complete data
write_csv(obis_gbif, "01_occurrences/complete/species.csv")

# 03 Check data integrity (Visualization) ------------------------------------------
# occ <- read.csv("01_occs/01_2 complete/species.csv")

# Static map (using ggplot2)
world <- annotation_borders("world", colour = "gray50", fill = "gray50") # world map borders

ggplot() +
  world +
  geom_point(
    data = obis_gbif,
    aes(
      longitude,
      latitude
    )
  ) +
  coord_fixed() +
  theme_minimal()


# Interactive map (with mapview) 
mapview::mapview(
  obis_gbif,
  xcol = "longitude",
  ycol = "latitude",
  crs = 4326,
  grid = FALSE
)