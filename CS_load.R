library(sf)
library(rgdal)
library(tidyverse)
library(doParallel)
library(raster)
library(randomForest)
library(rasterVis)
library(lattice)
library(latticeExtra)
library(viridis)
library(RColorBrewer)
library(sp)
library(spdplyr)
library(ROCR)
library(caret)
library(gridExtra)
# library(tree)
# library(gbm)
# library(ranger)       # a faster implementation of randomForest
# library(caret)

# projection used throughout: https://epsg.io/5070

# List of national extents CDL rasters downloaded to machine:
cdl_list <- Sys.glob("C:/Users/eburchf/Desktop/Data/CDL/*.img")
reference_raster <- raster(cdl_list[1])

state_list <- c("ME", "NH", "VT", "MA", "RI", "CT", "NY", "PA", "NJ", "DE", "MD", "WV", "VA", "NC", "SC", "KY",
                "OH", "MI", "WI", "IA", "MN", "MO", "TN", "GA", "AL", "MS", "LA", "AR", "KS", "OK", 
                "TX", "ND", "NE", "SD", "IL", "IN", "FL", "DC")
state_fips <- read.csv("C:/Users/eburchf/OneDrive - Emory University/Data/SHP/National/state_fips.csv")
state_fips$STATEFP <- str_pad(state_fips$ï..STATE, width = 2, side = "left", pad = "0")
state_fips <- state_fips %>% filter(STUSAB %in% state_list)
state_fips <- unique(state_fips$STATEFP)

# Shapefile of states in the contiguous US
states <- st_read("C:/Users/eburchf/OneDrive - Emory University/Data/SHP/National/states.shp", stringsAsFactors = F, quiet=T)
states <- st_transform(states, reference_raster@crs) 
states <- states %>% filter(!STATE_NAME %in% c("Hawaii", "Alaska"))
states <- states %>% filter(STATE_FIPS %in% state_fips)  # eliminate filters if/as scaling to national

# reproject for visualizations
r <-raster("./out/hist_grids/prob_ALFALFA_BIO.tif")
state_sp <- as(states, "Spatial")
state_sp <- spTransform(state_sp, r@crs)

# Shapefile of counties in the contiguous US
counties <- st_read("C:/Users/eburchf/OneDrive - Emory University/Data/SHP/National/county.shp", quiet=T) 
counties <- st_transform(counties, reference_raster@crs)
counties <- counties %>% filter(STATEFP %in% state_fips) # filter counties based on state_list

# Entire region of interest
roi <- st_read("./out/roi.shp", quiet=T)
roi <- st_transform(roi, reference_raster@crs)

# Crop information
crop_list <- c(1, 2, 5, 23, 24, 36, 37) # consider expanding to double-cropping
clean_crop_names <- c("Corn", "Cotton", "Soy", "Spring wheat", "Wheat", "Alfalfa", "Hay")
crop_names <- c("CORN", "COTTON", "SOY", "SWHEAT" , "WWHEAT",  "ALFALFA", "HAY" )
crop_dc_list <- list(c(1, 12, 13, 225, 226, 228, 237, 241), # corn  
                  c(2, 232, 238, 239), # cotton
                  c(5, 239, 240, 241, 26, 254), # soy
                  c(23), # swheat
                  c(24, 225, 236, 238, 26), # wwheat
                  c(36), # alfalfa
                  c(37)) # hay

#https://www.nature.com/articles/nmeth.1618
wong_colors <- c("#E69F00", "#e34234", "#009E73", "#F0E442", "#0072B2", "brown", "black")
names(wong_colors) <- clean_crop_names
colors_f <- scale_fill_manual(name = "CROP", values = wong_colors)
colors_c <- scale_color_manual(name = "CROP", values = wong_colors)

# AP points
ap <- readRDS("./out/ap.RDS")  # same projection as reference_raster
apm <- readRDS("./out/apm.RDS")
apmc <- readRDS("./out/apmc.RDS")
# apmc <- readRDS("./out/apmc_geoid.RDS")

project_theme <- theme_bw() +
  theme(text = element_text(size = 12))

map_theme <- project_theme +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        legend.position = "bottom")

# clean column names for results
cn <- c("AP", "Mean annual temp.", "Mean diurnal range", "Isothermality", "Temp. seasonality", "Max warm temp.", "Max cold temp.", 
        "Min cold temp.", "Temp. range")

# cleaned Census data
census_sf <- readRDS("./out/census_clean.RDS")
census_sf <- st_transform(census_sf, st_crs(apm))

# subset categories
topo <- c("SLOPE", "ELEVATION")
climate_sub <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY")
climate_full <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE", "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", "B7_TEMP_RANGE",
             "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD", "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", 
             "B16_PPT_WETQ", "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")
soil_full <- c("AWC_CLASS", "DRAINAGE", "REF_DEPTH", "S_PH_H2O", "T_BS", "T_CACO3", "T_CEC_CLAY", "T_CEC_SOIL", "T_ESP", "T_GRAVEL", "T_OC",
               "T_REF_BULK_DENSITY", "T_SILT", "T_TEB", "T_TEXTURE", "T_USDA_TEX_CLASS", "T_ECE", "T_SAND", "T_CLAY")
soil_sub <- c("T_GRAVEL", "T_SAND", "T_CLAY", "T_SILT", "T_CEC_CLAY", "T_CEC_SOIL", "DRAINAGE")
irr <- c("IRR")
inputs <- c("fert", "chem", "labor_expense", "machinery")
farm_res <- c("income", "gvt_prog", "insur_acres")
farm_char <- c("exp", "occup", "tenant", "acres_per_op")

climate_s <- c("B1_MEAN_ANNUAL_TEMP", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", "B2_MEAN_DIURNAL_RANGE",
               "B12_TOTAL_PPT", "B15_PPT_SEASONALITY", "B13_PPT_WET", "B14_PPT_DRY")

varlist <- list()
# varlist[["Climate and soil - FULL"]] <-  c(topo, climate_full, soil_full)
varlist[["Biophysical"]] <- c(topo, climate_sub, soil_sub)
# varlist[["All predictors - FULL"]] <- c(topo, climate_full, soil_full, irr, "PERC_IRR", farm_char, farm_res, inputs)
varlist[["Full"]] <- c(topo, climate_sub, soil_sub, irr, farm_char, farm_res, inputs)
# varlist[["All predictors - % irrigated"]] <- c(topo, climate_sub, soil_sub, "PERC_IRR", farm_char, farm_res, inputs)
# varlist[["All predictors - MIRAD irrigation"]] <- c(topo, climate_full, soil_full, irr, farm_char, farm_res, inputs)
varlist[["Agricultural"]] <- c(irr, farm_char, farm_res, inputs)
varlist[["Climate"]] <- c(climate_sub)

frr <- readRDS("./data/frr_shp.RDS")


# varlist <- list()
# varlist[["Climate"]] <- c(topo, temp, precip)
# varlist[["Biophysical"]] <- c(biophys)
# varlist[["Irrigation"]] <- c(biophys, irr)
# varlist[["Inputs"]] <- c(biophys, inputs)
# varlist[["Farm resources"]] <- c(biophys, farm_res)
# varlist[["Farm(er) characteristics"]] <- c(biophys, farm_char)
# varlist[["All predictors"]] <- c(biophys, irr, inputs, farm_res, farm_char)
# varlist[["All climate"]] <- c(topo, "B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE", "B3_ISOTHERMALITY", 
#                                   "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", "B7_TEMP_RANGE",
#                                   "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
#                                   "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ",
#                                   "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")

