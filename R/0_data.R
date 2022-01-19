# ----------------------------------------------------
# - Script for editing the data 
# Official data file: data/sild.csv
# Updates to the data should be saved to data/herring.Rdata 
# ----------------------------------------------------
library(tidyverse)

# -- formatting --
# data$dato <- as.Date(data$dato, format = "%Y-%m-%d")

# -- load data --
# data <- read.csv("data/hercond_final_addedMissingAges.csv")
data <- readRDS("data/herring_final.rds")
#load(file = "data/herring.Rdata")

# -----------------------------------------
# -------- make changes below here --------
# -----------------------------------------

#..add day of the year as column yday..
data <- mutate(data, yday = lubridate::yday(data$dato))


# -----------------------------------------
# -- update data --
saveRDS(data, file = "data/herring_final.rds")

