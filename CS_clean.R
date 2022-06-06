source("CS_load.R") 
source("CS_func.R")

################################################################################################################
#  Descriptive statistic extraction
################################################################################################################

# re-run this

all_states <- Sys.glob(paste0("C:/Users/eburchf/Desktop/Data/CDL/States/2019", "*", ".grd"))
ptm <- proc.time()
num_cores <- detectCores() - 8
clus <- makeCluster(num_cores)
registerDoParallel(clus)
out_data <- extract_from_state(num_cores, all_states)
stopCluster(clus)
out_df <- do.call(rbind.data.frame, out_data)
saveRDS(out_df, "./out/state_ds.RDS")

#################################################################################################################
#  Build crop frequency shapefiles
################################################################################################################

# generate random points in ROI
set.seed(300)
rpt <- st_as_sf(spsample(roi, n = 2500000, "random")) # big n needed to ensure 10K AP for each crop

state_names <- unique(states$STATE_NAME)
num_cores <- detectCores() - 3
clus <- makeCluster(num_cores)
registerDoParallel(clus)
out <- build_freq_points(num_cores, state_names)
stopCluster(clus)

# frequency count through time
ap <- st_set_geometry(out %>% select(contains("x")), NULL)

fdf <- list()
for (i in 1:length(crop_dc_list)) {
  # fdf[[i]] <- as.vector(rowSums(ap == crop_list[i], na.rm = T))
  # for double cropping:
  vls <- crop_dc_list[[i]]
  makec <- function(r) paste('ap==', r, '|')
  stmt <- paste((unlist(lapply(vls, FUN = makec))), collapse='')
  stmt <- substr(stmt, 1, nchar(stmt)-1)
  fdf[[i]] <- as.vector(rowSums(eval(parse(text=stmt)), na.rm=T))
}
fdf <- do.call("cbind.data.frame", fdf)
colnames(fdf) <- crop_names

# bind freq data and recoerce to sf object
final <- st_as_sf(cbind.data.frame(out, fdf))

# remove weak presence rows
corn <- final %>% filter(!Corn %in% c(1:2))
cotton <- final %>% filter(!Cotton %in% c(1:2))
soy <- final %>% filter(!Soy %in% c(1:2))
swheat <- final %>% filter(!Swheat %in% c(1:2))
wwheat <- final %>% filter(!Wwheat %in% c(1:2))
alfalfa <- final %>% filter(!Alfalfa %in% c(1:2))
hay <- final %>% filter(!Hay %in% c(1:2))
freq <- list(corn, cotton, soy, swheat, wwheat, alfalfa, hay)
# for (i in 1:length(freq)) {
#   df <- as.data.frame(freq[i])
#   v <- table(df %>% select(crop_names[i]))
#   print(crop_names[i])
#   print(v[1])
#   print(sum(v[2:nrow(v)]))
# }
# select 10K presence/absence
n <- 10000
for (i in 1:length(freq)) {
  df <- as.data.frame(freq[i])
  a <- df %>% filter(!!as.symbol(crop_names[i]) == 0)
  a <- a[sample(nrow(a), n), ]
  a$TYPE <- 0
  p <- df %>% filter(!!as.symbol(crop_names[i]) > 0)
  p <- p[sample(nrow(p), n), ]
  p$TYPE <- 1
  f <- rbind.data.frame(a, p)
  print(paste0(crop_names[i], nrow(f)))
  saveRDS(f, paste0("./out/", crop_names[i], "_AP.RDS")) # write out AP files
}

#################################################################################################################
#  Merge into single shapefile
################################################################################################################

corn <- clean_AP_data("Corn")
corn$CROP <- "CORN"
cotton <- clean_AP_data("Cotton")
cotton$CROP <- "COTTON"
soy <- clean_AP_data("Soy")
soy$CROP <- "SOY"
swheat <- clean_AP_data("Swheat")
swheat$CROP <- "SWHEAT"
wwheat <- clean_AP_data("Wwheat")
wwheat$CROP <- "WWHEAT"
alfalfa <- clean_AP_data("Alfalfa")
alfalfa$CROP <- "ALFALFA"
hay <- clean_AP_data("Hay")
hay$CROP <- "HAY"

ap <- rbind.data.frame(corn, cotton, soy, swheat, wwheat, alfalfa, hay)
ap$FID <- as.numeric(1:nrow(ap))
saveRDS(ap, "./out/ap.RDS")

#################################################################################################################
#  Extract slope and elevation
################################################################################################################

# s <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM_SLOPE.tif")
# d <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM.tif")
# topo <- stack(s, d)
# aps <- st_transform(ap, s@crs)
# 
# topo_ext <- extract(topo, aps, df=T)
# topo_ext$FID <- aps$FID
# colnames(topo_ext) <- c("ID", "SLOPE", "ELEVATION", "FID")
# topo_ext$ELEVATION[is.na(topo_ext$ELEVATION)] <- 0
# topo_ext <- topo_ext %>% dplyr::select(-c(ID))
# saveRDS(topo_ext, "./out/topo_extract.RDS")

#################################################################################################################
#  Extract MIRAD to points as frequency
################################################################################################################

# # 2012 as "baseline"
# mirad_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/MIRAD/mirad1km*.tif")
# mirad <- stack(mirad_list)
# apm <- st_transform(ap, mirad@crs)
# 
# mirad_ext <- raster::extract(mirad, apm, df=T)
# mirad_ext$FID <- apm$FID
# colnames(mirad_ext) <- c("ID", "IRR_02", "IRR_07", "IRR_12", "IRR_17", "FID")
# mirad_ext <- mirad_ext %>% select(-c("ID"))
# mirad_ext$IRR_AVG <- rowMeans(mirad_ext[,c("IRR_02", "IRR_07", "IRR_12", "IRR_17")])
# mirad_ext$IRR_ZERO <- ifelse(mirad_ext$IRR_AVG > 0, 1, 0)
# mirad_ext$IRR_50 <- ifelse(mirad_ext$IRR_AVG > 50, 1, 0)
# saveRDS(mirad_ext, "./out/mirad_extract.RDS")

#################################################################################################################
#  Extract DayMet
################################################################################################################

# daymet_download.R

# library(dismo)
# library(sf)
# 
# years <- 1999:2019
# 
# pptf <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/prcp_monttl_*")
# tminf <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/tmin_monavg*")
# tmaxf <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/tmax_monavg*")
# 
# for (i in 1:length(years)) {
#   print(years[i])
#   dir.create(file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",years[i])))
#   rasterOptions(tmpdir=file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",years[i])))
#   
#   ppt <- stack(pptf[i])
#   tmin <- stack(tminf[i])
#   tmax <- stack(tmaxf[i])
#   
#   roic <- spTransform(roi, tmin@crs)
#   
#   ppt <- crop(ppt, roic)
#   tmin <- crop(tmin, roic)
#   tmax <- crop(tmax, roic)
#   
#   bc <- biovars(ppt, tmin, tmax) # takes a sec, 1:07
#   writeRaster(bc, paste0("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/biovars_", years[i],".grd"), overwrite=T) 
#   unlink(paste0("C:/Users/eburchf/Desktop/TRASH/",years[i]), recursive = TRUE, force = TRUE) 
# }
# 
# 
# # average biovars rasters

# bv_list <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/biovars_*.grd"))
# r <- stack(bv_list[1])
# 
# for (i in 1:length(bv_list)) {
#   
#   r <- stack(bv_list[i])
#   print(r@crs)
#   
# }
# 
# # average of each slice of stack through time
# 
# for (i in 1:dim(r)[3]) {
# 
#   x <- stack()
#   
#   for (f in 1:length(bv_list)) {
# 
#     bvr <- stack(bv_list[f])
#     
#     bvrs <- bvr[[i]]
#     bvrs@crs <- r@crs
#     x <- stack(x, bvrs)
# 
#   }
# 
#   mn <- calc(x, fun = mean, na.rm=T)
#   i <- str_pad(i, 2, "left", "0")
#   writeRaster(mn, paste0("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/mn_99_19/bv", i, ".grd"), overwrite=T)
#   remove(x, mn)
# 
# }
# 
# # extract to points
# biovar_list <- sort(Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/mn_99_19/*.grd"))
# biovar <- stack(biovar_list)
# # projection information found here: https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html
# raster::projection(biovar) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314706705 +units=m +no_defs"
# apm <- st_transform(ap, biovar@crs)
#                    
# bv_ext <- raster::extract(biovar, apm, df=T)
# bv_ext$FID <- apm$FID
# colnames(bv_ext) <- c("ID", "BV1", "BV2", "BV3", "BV4", "BV5", "BV6", "BV7", "BV8",
#                       "BV9", "BV10", "BV11", "BV12", "BV13", "BV14", "BV15", "BV16", "BV17", "BV18", "BV19", "FID")
# bv_ext <- bv_ext %>% select(-c(ID))
# saveRDS(bv_ext, "./out/biovars_extract.RDS")


#################################################################################################################
#  HWSD madness
################################################################################################################

# # check out hwsd.R for the wizardry. 
# soil_data <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/rasters/*.tif")
# soil_ras <- stack(soil_data)
# 
# aps <- st_transform(ap, soil_ras@crs)
# 
# soil_ext <- raster::extract(soil_ras, aps, df=T)
# soil_ext$FID <- aps$FID
# soil_ext$T_ECE <- soil_ext$TS_ECE
# soil_ext$T_SAND <- soil_ext$TS_SAND
# soil_ext <- soil_ext %>% dplyr::select(-c(ID, TS_ECE, TS_SAND))
# 
# saveRDS(soil_ext, "./out/soil_extract.RDS")

#################################################################################################################
#  Merge all extracted data
################################################################################################################

bv <- readRDS("./out/extracts/biovars_extract.RDS")
colnames(bv) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                    "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                    "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                    "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                    "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ", "FID")
irr <- readRDS("./out/extracts/mirad_extract.RDS")
topo <- readRDS("./out/extracts/topo_extract.RDS")
soil <- readRDS("./out/extracts/soil_extract.RDS")

apm <- ap
apm <- merge(apm, bv, by = "FID", all=T)
apm <- merge(apm, irr, by = "FID", all=T)
apm <- merge(apm, topo, by = "FID", all=T)
apm <- merge(apm, soil, by = "FID", all=T)
apm$IRR <- ifelse(apm$IRR_AVG > 0, 1, 0)

saveRDS(apm, "./out/apm.RDS")

#################################################################################################################
#  Clean up Census data
################################################################################################################

# compute average over period of interest, 2007, 2012 and 2017
census <- readRDS("./data/BSDS_standatt_operated_v3.RDS") %>% 
  filter(year %in% c(2007, 2012, 2017)) %>%
  select(-c(year, county, perc_t)) %>% # remove null model attr
  group_by(GEOID) %>% 
  summarize_all(funs(mean), na.rm=T)
census_sf <- merge(counties, census, by = "GEOID") %>% select(-c(STATEFP, COUNTYFP, COUNTYNS, AFFGEOID,
                                                                 LSAD, ALAND, AWATER, cty_cl, cty_clp, cty_al, cty_pe,
                                                                 ph_corn, perc_cl, perc_clp, perc_p, perc_pe, perc_o, perc_wle,
                                                                 cons_wet_acres, manure_acres, male, insur_op, labor_n,
                                                                 fert_acres, comm_sales, cons_wetlands, full_owner, herb_acres,
                                                                 income_farm, insect_acres, irrig, part_owner, labor_n, female,
                                                                 crop_sales, NAME)) %>% na.omit()
saveRDS(census_sf, "./out/census_clean.RDS")

apm_census <- st_intersection(apm, census_sf) %>% select(-c(GEOID))
saveRDS(apm_census, "./out/apmc.RDS")

# add more detailed irrigation
irr <- readRDS("C:/Users/eburchf/OneDrive - Emory University/WF/Irrigation/out/clean_merged_irrigation.RDS") %>% group_by(GEOID) %>%
  summarize(PERC_IRR = mean(PERC_IRR, na.rm=T))
apmc <- merge(apmc, irr, by = "GEOID")

saveRDS(apmc, "./out/apmc.RDS")


#################################################################################################################
# Build big raster brick for historical prediction (DayMet 1K resolution just like RF training data)
################################################################################################################

prj <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314706705 +units=m +no_defs"

biovar_list <- sort(Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/mn_99_19/*.grd"))
biovar <- stack(biovar_list)
projection(biovar) <- prj  # the original projection read in R seems to be wrong - clips out Maine and someother areas
roi_bv <- st_transform(roi, biovar@crs)
biovar <- raster::crop(biovar, roi_bv)
biovar <- raster::mask(biovar, mask = roi_bv)
names(biovar) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                   "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                   "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                   "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                   "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")

mirad_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/MIRAD/mirad1km*.tif")
mirad <- stack(mirad_list)
mirad_sum <- sum(mirad)
mirad_sum[mirad_sum > 0] <- 1  # anywhere irrigation has occured since 2002
names(mirad_sum) <- "IRR"
mirad <- projectRaster(mirad_sum, biovar[[1]], method = "ngb") # categorical
roi_mirad <- st_transform(roi, biovar@crs)
mirad <- raster::crop(mirad, roi_mirad)
mirad <- raster::mask(mirad, mask = roi_mirad)

s <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM_SLOPE.tif")
d <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM.tif")
topo <- stack(s, d)
topo <- projectRaster(topo, biovar[[1]], method = "bilinear")
roi_topo <- st_transform(roi, biovar@crs)
topo <- raster::crop(topo, roi_topo)
topo <- raster::mask(topo, mask = roi_topo)
names(topo) <- c("SLOPE", "ELEVATION")

# unprojected WGS64
soil_data <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/rasters/*.tif")
soil_in <- stack(soil_data)
soil_ras <- projectRaster(soil_in, biovar[[1]], method = "ngb") # some categorical values here need to preserve
writeRaster(soil_ras, "C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_DayMet_2.tif", overwrite=T)
soil_ras <- stack("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_DayMet_2.tif")
roi_soil <- st_transform(roi, biovar@crs)
soil_crop <- raster::crop(soil_ras, roi_soil)
soil_crop <- raster::mask(soil_crop, roi_soil)
names(soil_crop) <- names(soil_in)

all_preds <- stack(biovar, topo, mirad, soil_crop) 
writeRaster(all_preds, "./out/attr_grids/hist_hires.grd", bandorder="BIL", overwrite=T) # doesn't include Census

#################################################################################################################
#  Add Census data to historical stack
################################################################################################################

library(sf)
library(fasterize)
library(raster)

stk <- stack("./out/attr_grids/hist_hires.grd")
census_sf <- readRDS("./out/census_clean.RDS") %>% dplyr::select(-c(GEOID))
census_sf <- st_transform(census_sf, stk@crs)

census_stack <- list()

for (c in 1:(length(colnames(census_sf))-1)) {
  
  census_stack[[c]] <- fasterize(census_sf, stk[[1]], colnames(census_sf)[[c]])
  
}

cstack <- stack(census_stack)
names(cstack) <- colnames(census_sf)[!colnames(census_sf) %in% c("geometry")]

allstack <- stack(stk, cstack)
writeRaster(allstack, "./out/attr_grids/hist_hires_census.grd", bandorder="BIL", overwrite=T)


#################################################################################################################
# Build big raster brick for future prediction (future at lower ~5K resolution)
################################################################################################################

# reproject and resample to WorldClim projection
dir <- "C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/"  # direct download
biovar_list  <- list.files(dir, recursive=T)

fut <- brick(paste0(dir, biovar_list[[1]]))
roi_f <- st_transform(roi, fut@crs)
fut <- raster::crop(fut, roi_f)
biovar <- raster::mask(fut, roi_f)
prj <- biovar@crs

mirad_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/MIRAD/mirad1km*.tif")
mirad <- stack(mirad_list)
mirad_sum <- sum(mirad)
mirad_sum[mirad_sum > 0] <- 1  # anywhere irrigation has occured since 2002
names(mirad_sum) <- "IRR"
mirad <- projectRaster(mirad_sum, biovar[[1]], method = "ngb")
roi_mirad <- st_transform(roi, prj)
mirad <- raster::crop(mirad, roi_mirad)
mirad <- raster::mask(mirad, mask = roi_mirad)

s <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM_SLOPE.tif")
d <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM.tif")
topo <- stack(s, d)
topo <- projectRaster(topo, biovar[[1]], method = "bilinear")
roi_topo <- st_transform(roi, prj)
topo <- raster::crop(topo, roi_topo)
topo <- raster::mask(topo, mask = roi_topo)
names(topo) <- c("SLOPE", "ELEVATION")

# unprojected WGS64
soil_data <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/rasters/*.tif")
soil_in <- stack(soil_data)
# soil_ras <- projectRaster(soil_in, biovar[[1]], method = "ngb") # change to biliear? NO! very categorical!
# writeRaster(soil_ras, "C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_lores.tif", overwrite=T)
soil_ras <- stack("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_lores.tif")
roi_soil <- st_transform(roi, prj)
soil_crop <- raster::crop(soil_ras, roi_soil)
soil_crop <- raster::mask(soil_crop, roi_soil)
names(soil_crop) <- names(soil_in)

future_attr <- stack(mirad, topo, soil_crop)
# writeRaster(future_attr, "./out/attr_grids/future_attr.grd", bandorder="BIL", overwrite=T)

for (i in 1:length(biovar_list)) {
  
  print(i)
  fut <- brick(paste0(dir, biovar_list[[i]]))
  roi_f <- st_transform(roi, fut@crs)
  fut <- raster::crop(fut, roi_f)
  fut <- raster::mask(fut, roi_f)
  names(fut) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                     "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                     "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                     "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                     "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")
  fn <- strsplit(biovar_list[[i]], "/")[[1]][1]
  writeRaster(fut, paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_grids/", fn, ".grd"), bandorder = "BIL", overwrite=T)

}


#################################################################################################################
#  DayMet resample (for error)
################################################################################################################

prj <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314706705 +units=m +no_defs"

biovar_list <- sort(Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/DayMet/mn_99_19/*.grd"))
biovar <- stack(biovar_list)
projection(biovar) <- prj  # the original projection read in R seems to be wrong - clips out Maine and someother areas
roi_bv <- st_transform(roi, biovar@crs)
biovar <- raster::crop(biovar, roi_bv)
biovar <- raster::mask(biovar, mask = roi_bv)
names(biovar) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                   "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                   "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                   "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                   "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")

futlist <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_grids/*.grd")
futref <- stack(futlist[[1]])

biovarp <- projectRaster(biovar, futref, method = "bilinear") # this does the resample

# build above
fut_attr <- stack("./out/attr_grids/future_attr.grd")
allras <- stack(biovarp, fut_attr)

census_sf <- readRDS("./out/census_clean.RDS") %>% dplyr::select(-c(GEOID))
census_sf <- st_transform(census_sf, allras@crs)

census_stack <- list()

for (c in 1:(length(colnames(census_sf))-1)) {
  
  census_stack[[c]] <- fasterize(census_sf, allras[[1]], colnames(census_sf)[[c]])
  
}

cstack <- stack(census_stack)
names(cstack) <- colnames(census_sf)[!colnames(census_sf) %in% c("geometry")]

allstack <- stack(allras, cstack)
writeRaster(allstack, "./out/attr_grids/5k_hist_attr_DayMet.grd", bandorder="BIL", overwrite=T)

#################################################################################################################
#  Historical WorldClim resample
################################################################################################################

biovar_list <- sort(Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/hist/*.tif"))
biovar <- stack(biovar_list)
roi_bv <- st_transform(roi, biovar@crs)
biovar <- raster::crop(biovar, roi_bv)
biovar <- raster::mask(biovar, mask = roi_bv)
names(biovar) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                   "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                   "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                   "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                   "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")


mirad_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/MIRAD/mirad1km*.tif")
mirad <- stack(mirad_list)
mirad_sum <- sum(mirad)
mirad_sum[mirad_sum > 0] <- 1  # anywhere irrigation has occured since 2002
names(mirad_sum) <- "IRR"
mirad <- projectRaster(mirad_sum, biovar[[1]], method = "ngb")
roi_mirad <- st_transform(roi, prj)
mirad <- raster::crop(mirad, roi_mirad)
mirad <- raster::mask(mirad, mask = roi_mirad)

s <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM_SLOPE.tif")
d <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM.tif")
topo <- stack(s, d)
topo <- projectRaster(topo, biovar[[1]], method = "bilinear")
roi_topo <- st_transform(roi, prj)
topo <- raster::crop(topo, roi_topo)
topo <- raster::mask(topo, mask = roi_topo)
names(topo) <- c("SLOPE", "ELEVATION")

# unprojected WGS64
soil_data <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/rasters/*.tif")
soil_in <- stack(soil_data)
soil_ras <- projectRaster(soil_in, biovar[[1]], method = "ngb") # change to biliear? NO! very categorical!
writeRaster(soil_ras, "C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_WCres.tif", overwrite=T)
soil_ras <- stack("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/CS_reproj_lores.tif")
roi_soil <- st_transform(roi, prj)
soil_crop <- raster::crop(soil_ras, roi_soil)
soil_crop <- raster::mask(soil_crop, roi_soil)
names(soil_crop) <- names(soil_in)

future_attr <- stack(biovar, mirad, topo, soil_crop)

census_sf <- readRDS("./out/census_clean.RDS") %>% dplyr::select(-c(GEOID))
census_sf <- st_transform(census_sf, allras@crs)

census_stack <- list()

for (c in 1:(length(colnames(census_sf))-1)) {
  
  census_stack[[c]] <- fasterize(census_sf, future_attr[[1]], colnames(census_sf)[[c]])
  
}

cstack <- stack(census_stack)
names(cstack) <- colnames(census_sf)[!colnames(census_sf) %in% c("geometry")]

allstack <- stack(future_attr, cstack)
writeRaster(allstack, "./out/attr_grids/5k_hist_attr_WC.grd", bandorder="BIL", overwrite=T)











# build above
fut_attr <- stack("./out/attr_grids/future_attr.grd")
allras <- stack(biovarp, fut_attr)

census_sf <- readRDS("./out/census_clean.RDS") %>% dplyr::select(-c(GEOID))
census_sf <- st_transform(census_sf, allras@crs)

census_stack <- list()

for (c in 1:(length(colnames(census_sf))-1)) {
  
  census_stack[[c]] <- fasterize(census_sf, allras[[1]], colnames(census_sf)[[c]])
  
}

cstack <- stack(census_stack)
names(cstack) <- colnames(census_sf)[!colnames(census_sf) %in% c("geometry")]

allstack <- stack(allras, cstack)
writeRaster(allstack, "./out/attr_grids/5k_hist_attr_DayMet.grd", bandorder="BIL", overwrite=T)


#################################################################################################################
#  Historical WorldClim resample
################################################################################################################

attr <- brick("./out/attr_grids/5k_hist_attr_WC.grd")
names <- c("Alfalfa", "Corn", "Cotton", "Hay", "Soy", "Spring wheat", "Winter wheat")
hs <- stack(Sys.glob("./out/hist_grids/prob*_BIO_DayMet.tif"))
# hs[is.na(hs)] <- 0 # in fill water with zero, but check this
names(hs) <- names
hs <- dropLayer(hs, "Spring.wheat") 

hsr <- resample(hs, attr, method = "bilinear")
names(hsr) <- hs
writeRaster(hsr, "./out/hist_grids/prob_BIO_DayMet_RS.tif", overwrite=T)




# Archive

#################################################################################################################
# Build big raster brick for historical prediction (past using WorldClim data and resolution)
################################################################################################################

biovar_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/hist/*.tif")


fut <- stack(paste0(dir, biovar_list))
roi_f <- st_transform(roi, fut@crs)
fut <- raster::crop(fut, roi_f)
fut <- raster::mask(fut, roi_f)
prj <- fut@crs
names(fut) <- c("B1_MEAN_ANNUAL_TEMP", "B2_MEAN_DIURNAL_RANGE",
                   "B3_ISOTHERMALITY", "B4_TEMP_SEASONALITY", "B5_MAXTEMP_WARM", "B6_MINTEMP_COLD", 
                   "B7_TEMP_RANGE", "B8_MEANTEMP_WET", "B9_MEANTEMP_DRY", "B10_MEANTEMP_WARM", "B11_MEANTEMP_COLD",
                   "B12_TOTAL_PPT", "B13_PPT_WET", "B14_PPT_DRY", "B15_PPT_SEASONALITY", "B16_PPT_WETQ", 
                   "B17_PPT_DRYQ", "B18_PPT_WARMQ", "B19_PPT_COLDQ")

mirad_list <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/MIRAD/mirad1km*.tif")
mirad <- stack(mirad_list)
mirad_sum <- sum(mirad)
mirad_sum[mirad_sum > 0] <- 1  # anywhere irrigation has occured since 2002
names(mirad_sum) <- "IRR"
mirad <- projectRaster(mirad_sum, fut[[1]], method = "ngb")
roi_mirad <- st_transform(roi, prj)
mirad <- raster::crop(mirad, roi_mirad)
mirad <- raster::mask(mirad, mask = roi_mirad)

s <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM_SLOPE.tif")
d <- raster("C:/Users/eburchf/OneDrive - Emory University/Data/DEM/NA_1KM.tif")
topo <- stack(s, d)
topo <- projectRaster(topo, fut[[1]], method = "bilinear")
roi_topo <- st_transform(roi, prj)
topo <- raster::crop(topo, roi_topo)
topo <- raster::mask(topo, mask = roi_topo)
names(topo) <- c("SLOPE", "ELEVATION")

# unprojected WGS64
soil_data <- Sys.glob("C:/Users/eburchf/OneDrive - Emory University/Data/HWSD/rasters/*.tif")
soil_in <- stack(soil_data)
soil_ras <- projectRaster(soil_in, fut[[1]], method = "ngb")
roi_soil <- st_transform(roi, prj)
soil_crop <- raster::crop(soil_ras, roi_soil)
soil_crop <- raster::mask(soil_crop, roi_soil)
names(soil_crop) <- names(soil_in)

future_attr <- stack(mirad, topo, soil_crop, fut)
writeRaster(future_attr, "./out/attr_grids/allhist_WC.grd", bandorder="BIL", overwrite=T)


