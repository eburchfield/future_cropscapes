build_freq_points <- function(cores, state_names) {
  foreach(i=1:length(state_names),
          .packages= c("velox","foreach", "sf", "tidyverse", "raster", "sp"),
          .combine = "rbind.data.frame",
          .export = c("state_names", "states", "rpt")) %dopar% { 
            
            dir.create(file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i])))
            rasterOptions(tmpdir=file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i])))
            
            # build state raster
            sras_files <- ifelse(state_names[i] == "Virginia",
                                 paste0("C:/Users/eburchf/Desktop/Data/CDL/States/*_Virginia.grd"),
                                 paste0("C:/Users/eburchf/Desktop/Data/CDL/States/*", state_names[i], ".grd"))
            sras_stack <- stack(Sys.glob(sras_files)) 
            
            # build state shapefile
            state_shp <- states %>% 
              filter(STATE_NAME == state_names[i]) %>% 
              dplyr::select(STATE_NAME, STATE_FIPS, geometry)
            
            # subset random sample points
            rpt_sub <- st_intersection(rpt, state_shp) 
            
            # extract to points
            ap <- st_as_sf(raster::extract(sras_stack, rpt_sub, sp = T, cellnumbers = T, method = "simple"))
            unlink(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i]), recursive = TRUE, force = TRUE) 
            return(ap)
          }
}

extract_from_state <- function(cores, all_states){
  
   foreach(i=1:length(all_states), 
           
          .packages= c("raster","foreach", "sf", "tidyverse"),
          #.combine = rbind.data.frame, 
          .export = c("all_states", "states")) %dopar% {

              year <- as.numeric(substr(all_states[i], 42, 45)) # check subsetting
              state_name <- gsub(".grd", "", sub(paste0(".*", year,"_"), "", all_states[i]))
              ras <- raster(all_states[i])  # could take brick here
              
              dir.create(file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",year,state_name))) #creates unique filepath for temp directory
              rasterOptions(tmpdir=file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",year,state_name)))  #sets tmpdir
              
              s <- states %>% filter(STATE_NAME == state_name)
              s <- as(s, "Spatial")
              
              r <- raster::crop(ras, s)
              r <- raster::mask(r, s)
              
              vals <- raster::getValues(r)
              prop <- as.data.frame(table(vals))
              colnames(prop) <- c("CROP_VALUE", "FREQ")
              prop$N <- length(vals[!is.na(vals)])
              prop$PERC <- (prop$FREQ/prop$N)*100
              prop$STATE_NAME <- state_name
              prop$YEAR <- year

              unlink(paste0("C:/Users/eburchf/Desktop/TRASH/",year,state_name), recursive = TRUE, force = TRUE) # delete temporary files
              prop
              
          }
}

clean_AP_data <- function(crop_name) {
  
  ap <- readRDS(paste0("./out/", crop_name, "_AP.RDS")) 
  ap <- st_as_sf(ap)
  
  ap <- ap %>%
    dplyr::select(STATE_NAME, STATE_FIPS, !!as.symbol(crop_name), TYPE, geometry)
  
  colnames(ap) <- c("STATE_NAME", "STATE_FIPS", "FREQ", "AP", "geometry")
  return(ap)
}

run_historical_models <- function(mtry_dir, out_ext, data_list) {
  
  # runs historical RF models on full panel
  
  mtry <- readRDS(mtry_dir)
  
  for (i in 1:length(crop_names)) {
    crop_name <- crop_names[[i]]
    m <- mtry %>% 
      filter(CROP == crop_name) %>% 
      pull(MTRYTRAIN)
    m <- ifelse(crop_name == "HAY", m[1], m)
    
    train <- data_list[[crop_name]]
    
    rf_ap <- randomForest(AP ~ ., 
                          data = train,
                          mtry = m, 
                          ntree = 500)
    
    saveRDS(rf_ap, paste0("./out/RF_historical/", crop_name, out_ext))
    
  }
  
}

rf_performance <- function(test_list, out_ext) {
  
  perf_final <- list()
  
  for (i in 1:length(crop_names)) {
    
    crop_name <- crop_names[[i]]
    test <- test_list[[crop_name]]
    
    rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext))
    pred <- predict(rf, test, type = "class")
    t <- with(test, table(pred, test$AP))
    rf_perf <- (t[1,1] + t[2,2])/nrow(test)
    out <- cbind.data.frame(PERF = rf_perf, CROP = crop_name)
    perf_final[[i]] <- out
    
  }
  
  perf_final <- do.call("rbind.data.frame", perf_final)
  
  return(perf_final)
  
  
}

build_imp <- function(out_ext) {
  
  imp_final <- list()
  
  for (i in 1:length(crop_names)) {
    
    crop_name <- crop_names[[i]]
    
    rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext))
    imp <- as.data.frame(randomForest::importance(rf))
    imp$varnames <- rownames(imp) # row names to column
    imp$CROP <- crop_name
    
    imp_final[[i]] <- imp
    
  }
  
  imp_final <- do.call("rbind.data.frame", imp_final)
  return(imp_final)
}

build_historical_grids <- function(out_ext, vars) {
  
  hist_preds <- brick("./out/attr_grids/5k_hist_attr_WC.grd") # remove WC, go with daymet or hi_res for viz
  hist_preds[["T_SAND"]] <- hist_preds[["TS_SAND"]]
  
  for (i in 1:length(crop_names)) {
    
    crop_name <- crop_names[[i]]
    print(crop_name)
    rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext))
    full <- full_list[[crop_name]]
    apsub <- subset(hist_preds, vars)
    prob <- predict(apsub, rf, type = "prob", index = 2) # prob of value 1, col 2
    writeRaster(prob, paste0("./out/hist_grids/prob_", crop_name, substr(out_ext, 1, nchar(out_ext)-4), "_WC.tif"), overwrite=T)
    
  }
}

build_future_stack <- function(dir) {
  
  names <- c("Alfalfa", "Corn", "Cotton", "Hay", "Soy", "Spring wheat", "Winter wheat")
  m1 <- stack(Sys.glob(dir))
  names(m1) <- names
  m1 <- dropLayer(m1, "Spring.wheat")
  return(m1)
  
}

run_human <- function(crop_name, nint){{
  
  crop_perf <- list()
  crop_imp <- list()
  crop_name <- crop_names[[c]]
  
  for (i in 1:length(varlist)) {
    
    lname <- names(varlist)[[i]]

    census_sub <- full_list[[crop_name]] %>% select(AP, varlist[[i]])
    
    perf_list <- list()
    imp_list <- list()
    
    for (s in 1:nint) {
      
      set.seed(s) # different random subsets of data

      random_rn <- sample(nrow(census_sub), ceiling(nrow(census_sub)*.25))
      train <- census_sub[-random_rn,]
      test <- census_sub[random_rn,] 
      rf_ap <- randomForest(AP ~ ., data = train, ntree = 500)  # change to other alrogithms
      
      # percentage correctly classified, PCC
      pred <- predict(rf_ap, test, type = "class")
      t <- with(test, table(pred, test$AP))
      pcc <- (t[1,1] + t[2,2])/nrow(test)

      # sensitivity (percent of presences correctly classified)
      cm <- caret::confusionMatrix(pred, test$AP)
      sens <- cm$byClass["Sensitivity"]
      
      # specificity (the percentage of absenses correctly classified)
      spec <- cm$byClass["Specificity"]
      
      # kappa (a measure of the agreement between predicted presences and absenses corrected for agreement that might be due to chane alone)
      kappa <- cm$overall["Kappa"]
      
      # auc (area under ROC)
      pred2 <- prediction(as.numeric(pred), test$AP)
      auc <- performance(pred2, measure = "auc")
      auc <- auc@y.values
      

      perf <- cbind.data.frame(pcc, sens, spec, kappa, auc)
      rownames(perf) <- NULL
      colnames(perf) <- c("PCC", "SENSITIVITY", "SPECIFICITY", "KAPPA", "AUC")
      perf_list[[s]] <- perf
      
      # importance
      imp <- as.data.frame(randomForest::importance(rf_ap))
      imp$varnames <- rownames(imp) # row names to column
      imp$varnames <- as.factor(imp$varnames)
      rownames(imp) <- NULL
      imp$ITERATION <- s
      imp_list[[s]] <- imp
      
    }
    
    perf <- do.call(rbind.data.frame, perf_list)
    # colnames(perf) <- c("PCC", "SENSITIVITY", "SPECIFICITY", "KAPPA", "AUC")
    
    imp <- do.call(rbind.data.frame, imp_list)
    colnames(imp) <- c("MEANDECREASEGINI", "VARIABLES", "ITERATION")
    
    crop_perf[[lname]] <- perf
    crop_imp[[lname]] <- imp
    
  }
  
  perf_df <- do.call(rbind.data.frame, crop_perf)
  perf_df$VARSET <- rownames(perf_df)
  rownames(perf_df) <- NULL
  perf_df <- perf_df %>% separate(VARSET, c("VARSET", "ITERATION"), sep = "[.]")
  perf_df$CROP <- crop_name
  # all_perf[[crop_name]] <- perf_df
  saveRDS(perf_df, paste0("./out/human/", crop_name, "_perf_final_wclim.RDS"))
  
  imp_df <- do.call(rbind.data.frame, crop_imp)
  imp_df$VARSET <- rownames(imp_df)
  rownames(imp_df) <- NULL
  imp_df <- imp_df %>% separate(VARSET, c("VARSET", "remove"), sep = "[.]") %>% select(-c(remove))
  imp_df$CROP <- crop_name
  # all_imp[[crop_name]] <- imp_df
  saveRDS(imp_df, paste0("./out/human/", crop_name, "_imp_final_wclim.RDS"))
  
}}

run_those_models <- function(crop_names, nint) {
  foreach(c = 1:length(crop_names),
          
          .export = c("varlist", "nint", "full_list", "run_human"),
          .packages = c("randomForest", "tidyverse", "ROCR")) %dopar% {
            
            crop_name = crop_names[[c]]
            run_human(crop_name, nint)
            
          }
  
}

build_hist_vis <- function(out_ext, yl, ds = "_DayMet") {
  
  all_ras <- list()
  library(viridis)
  coul2 <- viridis::inferno(25)
  
  r <- stack(Sys.glob(paste0("./out/hist_grids/prob_*", substr(out_ext, 1, nchar(out_ext)-4), ds, ".tif")))
  r <- dropLayer(r, paste0("prob_SWHEAT", substr(out_ext, 1, nchar(out_ext)-4), ds))

  names(r) <- c("Alfalfa", "Corn", "Cotton", "Hay", "Soy", "Wheat")

  state_sp <- as(st_transform(states, r@crs), "Spatial")

  hist_clim <- rasterVis::levelplot(r, 
                         margin = F,
                         xlab = "",
                         ylab = yl,
                         main = "",
                         # xlab = list(names(all_ras), space = "top"),  # winter wheat
                         # ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2000-2020"), rot=0), 
                         scales = list(draw=F), 
                         layout = c(3, 2), 
                         names.attr = names(r), #rep("", dim(allvis)[3]),
                         col.regions = rev(terrain.colors(20)),
                         colorkey=list(space="bottom"),
                         at = seq(0, 1, length.out = 20)) + # inferno viridis
    latticeExtra::layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  #pdf(paste0("./figs/hist_", substr(out_ext, 1, nchar(out_ext)-4), ".pdf"))
  #print(hist_clim)
  #dev.off()
  out <- list(hist_clim, r)
  
  return(out)
  
  
}

build_error_rasters <- function(out_ext, source = "", fact = 3, thdf = bio_th) {
  
  # load necessary data
  actual <- readRDS("./out/cdl_extract.RDS")
  ras_list <- Sys.glob(paste0("./out/hist_grids/prob_*", out_ext, source, ".tif"))
  hs <- stack(ras_list) # reference raster 
  hs <- aggregate(hs, fact = fact, fun = mean, na.rm=T)
  actual <- st_transform(actual, hs@crs)
  
  # tha <- find_threshold(out_ext)
  
  crop_names_sub <- c("CORN", "COTTON", "SOY", "WWHEAT", "ALFALFA", "HAY")
  
  for (i in 1:length(crop_names_sub)) {
    
    crop <- crop_names_sub[[i]]
    print(crop)
    
    actual_ras <- rasterize(actual[,crop], hs[[1]], fun = sum, background = 0)
    writeRaster(actual_ras[[2]], paste0("./out/hist_error/extent_", crop, ".tif"), overwrite=T)
    
    pred_ras <- raster(paste0("./out/hist_grids/prob_", crop, "_", out_ext, source, ".tif"))
    pred_ras <- aggregate(pred_ras, fact = fact, fun = mean, na.rm=T)
    
    crop_error <- data.frame(PREDV = getValues(pred_ras), ACTV = getValues(actual_ras[[2]]))
    crop_error$ACTV[is.na(crop_error$PREDV)] <- NA
    
    th <- thdf %>% filter(CROP == crop) %>% pull(TH)
 
    crop_error$AS <- ifelse(crop_error$PREDV >= th  & crop_error$ACTV > 0, 1, 0) # correct classification of suitability
    crop_error$FS <- ifelse(crop_error$PREDV >= th  & crop_error$ACTV == 0, 1, 0) # falsely classified suitability
    crop_error$MS <- ifelse(crop_error$PREDV <= th & crop_error$ACTV > 0, 1, 0) # falsely classified unsuitability
    crop_error$AUS <- ifelse(crop_error$PREDV <= th & crop_error$ACTV == 0, 1, 0) # correctly classified unsuitability
    
    crop_error$CATS <- ifelse(crop_error$AS == 1, 1, 0) # green
    crop_error$CATS <- ifelse(crop_error$FS == 1, 2, crop_error$CATS) # red
    crop_error$CATS <- ifelse(crop_error$MS == 1, 3, crop_error$CATS) # blue
    crop_error$CATS <- ifelse(crop_error$AUS == 1, 4, crop_error$CATS) # gray
    
    out_ras <- setValues(actual_ras[[2]], crop_error$CATS)
    writeRaster(out_ras, paste0("./out/hist_error/cat_", crop, "_", out_ext, "_", fact, "_CSTH.tif"), overwrite=T)
    
  }
  
}

build_future_grids <- function(out_ext = "_BIO.RDS", full_list) {
  
  attr <- brick("./out/attr_grids/future_attr.grd") # just resampled to new proj
  attr[["T_SAND"]] <- attr[["TS_SAND"]]
  fut <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_grids/*.grd")) # cropped, masked and named
  
  for (j in 1:length(fut)) {
    
    fr <- brick(fut[[j]])
    fn <- strsplit(fut[[j]], "/")[[1]][9]
    fn <- substr(fn, 1, nchar(fn)-4)
    all_attr <- stack(fr, attr)
    print(fn)
    
    for (i in 1:length(crop_names)) {
      
      crop_name <- crop_names[[i]]
      print(crop_name)
      rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext))
      apsub <- subset(all_attr, names(full_list[[crop_name]]))
      prob <- predict(apsub, rf, type = "prob", index = 2) # prob of value 1, col 2
      writeRaster(prob, paste0("./out/future_grids/prob_", fn, "_", crop_name, substr(out_ext, 1, nchar(out_ext)-4), ".tif"), overwrite=T)
      
    }
    
  }
  
}

build_ensemble_figure <- function(ssp) {
  
  names <- c("Alfalfa 2021-2040", "Alfalfa 2041-2060", "Alfalfa 2061-2080", "Alfalfa 2081-2100", 
             "Corn 2021-2040", "Corn 2041-2060", "Corn 2061-2080", "Corn 2081-2100", 
             "Cotton 2021-2040", "Cotton 2041-2060", "Cotton 2061-2080", "Cotton 2081-2100", 
             "Hay 2021-2040", "Hay 2041-2060", "Hay 2061-2080", "Hay 2081-2100", 
             "Soy 2021-2040", "Soy 2041-2060", "Soy 2061-2080", "Soy 2081-2100", 
             "Wheat 2021-2040", "Wheat 2041-2060", "Wheat 2061-2080", "Wheat 2081-2100") 
  
  all_mean <- stack(Sys.glob(paste0("./out/future_ens/mean_stack_", ssp, "*.tif")))
  names(all_mean) <- names
  all_sd <- stack(Sys.glob(paste0("./out/future_ens/sd_stack_", ssp, "*.tif")))
  names(all_sd) <- names
  
  state_sp <- as(st_transform(states, all_mean@crs), "Spatial")
  
  future_clim <- levelplot(all_mean, 
                           margin = F,
                           xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                           ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                           scales = list(draw=F), 
                           layout = c(4, 6), 
                           names.attr = rep("", dim(all_mean)[3]),
                           col.regions = rev(terrain.colors(20)),
                           at = seq(0, 1, length.out = 20)) + 
    layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  pdf(paste0("./figs/future_probs_ens_mean_", ssp, ".pdf"))
  print(future_clim)
  dev.off()
  
  coul2 <- colorRampPalette(brewer.pal(11, "BuPu"))
  
  future_clim <- levelplot(all_sd, 
                           margin = F,
                           xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                           ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                           scales = list(draw=F), 
                           layout = c(4, 6), 
                           names.attr = rep("", dim(all_mean)[3]),
                           col.regions = coul2) +
    #  at = seq(0, 1, length.out = 20)) + # inferno viridis
    layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  pdf(paste0("./figs/future_probs_ens_sd_", ssp, ".pdf"))
  print(future_clim)
  dev.off()
  
}

build_ensemble_delta_figure <- function(ssp, out_ext = "_BIO") {
  
  names <- c("Alfalfa 2021-2040", "Alfalfa 2041-2060", "Alfalfa 2061-2080", "Alfalfa 2081-2100", 
             "Corn 2021-2040", "Corn 2041-2060", "Corn 2061-2080", "Corn 2081-2100", 
             "Cotton 2021-2040", "Cotton 2041-2060", "Cotton 2061-2080", "Cotton 2081-2100", 
             "Hay 2021-2040", "Hay 2041-2060", "Hay 2061-2080", "Hay 2081-2100", 
             "Soy 2021-2040", "Soy 2041-2060", "Soy 2061-2080", "Soy 2081-2100", 
             "Wheat 2021-2040", "Wheat 2041-2060", "Wheat 2061-2080", "Wheat 2081-2100") 
  
  all_mean <- stack(Sys.glob(paste0("./out/future_ens/mean_stack_", ssp, "*.tif")))
  names(all_mean) <- names
  
  hs <- stack(Sys.glob("./out/hist_grids/prob_*_FULL_WC.tif")) # compare with FULL since best performing model
  hs <- dropLayer(hs, "prob_SWHEAT_FULL_WC")
  
  hsf <- stack(stack(replicate(4, hs[[1]])), stack(replicate(4, hs[[2]])), stack(replicate(4, hs[[3]])), stack(replicate(4, hs[[4]])), 
               stack(replicate(4, hs[[5]])), stack(replicate(4, hs[[6]])))
  names(hsf) <- names
  
  state_sp <- as(st_transform(states, all_mean@crs), "Spatial")
  
  allvis <- all_mean - hsf
  
  coul2 <- colorRampPalette(brewer.pal(11, "BrBG"))
  
  delta_clim <- levelplot(allvis, 
                          margin = F,
                          xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                          ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                          scales = list(draw=F), 
                          layout = c(4, 6), 
                          names.attr = rep("", dim(allvis)[3]),
                          col.regions = coul2,
                          at = seq(-0.6,  0.6, length.out = 11)) +
    layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  g <- grid.arrange(delta_clim, ncol=1, nrow=1)
  ggsave(paste0("./figs/future_diff_", ssp, out_ext, "_WC.png"), g, width=10, height=10)
  
}

find_threshold <- function(out_ext, index = "Accuracy") {

  th_final <- list()

  for (i in 1:length(crop_names)) {

    crop = crop_names[[i]]

    rf <- readRDS(paste0("./out/RF_historical/", crop, out_ext))
    pred <- cbind.data.frame(predict(rf, test_list[[i]], type = "prob"), test_list[[i]]$AP)
    colnames(pred) <- c("ABSENSE", "PRESENCE", "ACTUAL")
    pred$ACTUAL <- as.factor(pred$ACTUAL)
    levels(pred$ACTUAL) <- c(0, 1)

    th_list <- seq(0.05, 1.0, by = 0.05)
    th_out <- list()

    for (t in 1:(length(th_list))) {

      p <- pred
      p$PRED <- as.factor(ifelse(p$PRESENCE > th_list[[t]], 1, 0))
      levels(p$PRED) <- c(0,1)

      out <- confusionMatrix(p$PRED, p$ACTUAL)
      out <- as.data.frame(out$overall)
      out$INDEX <- rownames(out)
      out$TH <- th_list[[t]]

      th_out[[t]] <- out
    }

    ops_out <- do.call("rbind.data.frame", th_out)
    colnames(ops_out) <- c("VALUE", "INDEX", "TH")
    rownames(ops_out) <- NULL

    max_th <- ops_out %>% filter(INDEX == rlang::sym(index)) %>% filter(VALUE == max(VALUE, na.rm=T))
    max_th$CROP <- crop
    max_th$CCROP <- clean_crop_names[[i]]
    if (nrow(max_th) > 1) {
      max_th <- max_th %>% filter(TH == min(TH))

    } else {
      max_th <- max_th
    }
    th_final[[i]] <- max_th
  }


  th <- do.call("rbind.data.frame", th_final)
  
  th <- th %>% filter(CROP != "SWHEAT")

  return(th)

}

build_pdf_trimmed <- function(varname, out_ext = "_BIO") {
  
  pdf <- readRDS(paste0("./out/pdp/", varname, out_ext, ".RDS"))
  sd <- apmc %>%
    group_by(CROP) %>%
    summarize(hisd = (mean(!!rlang::sym(varname), na.rm=T) + 2*sd(!!rlang::sym(varname), na.rm=T)),
              losd = (mean(!!rlang::sym(varname), na.rm=T) - 2*sd(!!rlang::sym(varname), na.rm=T)))
  st_geometry(sd) <- NULL
  sd$CCROP <- NA
  
  for (i in 1:length(crop_names)) {
    
    sd$CCROP[sd$CROP == crop_names[[i]]] <- clean_crop_names[[i]]
    
  }
  
  pdf <- merge(pdf, sd, by.x ="CROP", by.y = "CCROP")

  hicotton <- apmc %>%
    group_by(CROP) %>%
    summarize(hisd = (mean(!!rlang::sym(varname), na.rm=T) + 1*sd(!!rlang::sym(varname), na.rm=T))) %>%
    filter(CROP == "COTTON") %>% pull(hisd)
  
  locotton <- apmc %>%
    group_by(CROP) %>%
    summarize(losd = (mean(!!rlang::sym(varname), na.rm=T) - 1*sd(!!rlang::sym(varname), na.rm=T))) %>%
    filter(CROP == "COTTON") %>% pull(losd)
  
  pdf$hisd[pdf$CROP == "COTTON"] <- hicotton
  pdf$losd[pdf$CROP == "COTTON"] <- locotton
  
  pdf$y[pdf$x > pdf$hisd] <- NA
  pdf$y[pdf$x < pdf$losd] <- NA
  
  pdf <- pdf %>% filter(CROP.y != "SWHEAT")
  
  return(pdf)
  
}

future_cat_viz <- function(ssp = "585", th = bio_th_k) {
  
  names <- c("Alfalfa 2021-2040", "Alfalfa 2041-2060", "Alfalfa 2061-2080", "Alfalfa 2081-2100", 
             "Corn 2021-2040", "Corn 2041-2060", "Corn 2061-2080", "Corn 2081-2100", 
             "Cotton 2021-2040", "Cotton 2041-2060", "Cotton 2061-2080", "Cotton 2081-2100", 
             "Hay 2021-2040", "Hay 2041-2060", "Hay 2061-2080", "Hay 2081-2100", 
             "Soy 2021-2040", "Soy 2041-2060", "Soy 2061-2080", "Soy 2081-2100", 
             "Wheat 2021-2040", "Wheat 2041-2060", "Wheat 2061-2080", "Wheat 2081-2100") 
  
  all_mean <- stack(Sys.glob(paste0("./out/future_ens/mean_stack_", ssp, "*.tif")))
  names(all_mean) <- names
  
  hs <- stack(Sys.glob("./out/hist_grids/prob_*_FULL_DayMet.tif")) 
  hs <- dropLayer(hs, "prob_SWHEAT_FULL_DayMet") 
  
  hsf <- stack(stack(replicate(4, hs[[1]])), stack(replicate(4, hs[[2]])), stack(replicate(4, hs[[3]])), stack(replicate(4, hs[[4]])), 
               stack(replicate(4, hs[[5]])), stack(replicate(4, hs[[6]])))
  names(hsf) <- names

  sras <- stack()
  
  for (i in 1:dim(all_mean)[[3]]) {
    
    print(names[[i]])
    
    # pull crop specific th
    cn <- strsplit(names[[i]], " ")[[1]][1]
    ths <- th %>% filter(CCROP == cn) %>% pull(TH)
    
    hist <- hsf[[i]]
    fut <- all_mean[[i]]
    
    hist[hist >= ths] <- 1
    hist[hist < ths] <- 0
    
    fut[fut >= ths] <- 11
    fut[fut < ths] <- 0
    
    # 11-1 = 10 remains suitable
    # 11-0 = 11 becomes suitable
    # 0-1 = -1 becomes unsuitable
    # 0-0 = 0 remains unsuitable

    r <- fut - hist
    r <- ratify(r)
    # rat <- levels(r)[[1]]
    rat <- data.frame(ID = c(-1, 0, 10, 11))
    
    rat$TYPE <- NA
    rat$TYPE[rat$ID == -1] <- "Becomes unsuitable" #red
    rat$TYPE[rat$ID == 0] <- "Remains unsuitable" # gray
    rat$TYPE[rat$ID == 10] <- "Remains suitable" # green
    rat$TYPE[rat$ID == 11] <- "Becomes suitable" # blue
    
    levels(r) <- rat
    sras <- stack(sras, r)
    
  }
  
  names(sras) <- names
  cols <- c("orange", "light gray", "#228B22", "purple")
  
  state_sp <- as(st_transform(states, sras@crs), "Spatial")
  
  fcat_viz <- levelplot(sras, 
                        margin = F,
                        xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                        ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                        scales = list(draw=F), 
                        layout = c(4, 6), 
                        names.attr = rep("", dim(sras)[3]),
                        colorkey = T,
                        col.regions = cols) +
    layer(sp.polygons(state_sp, lwd = 0.7, col = "black"))
  
  pdf(paste0("./figs/cat_fut_proj_", ssp, "_CTH.pdf"))
  print(fcat_viz)
  dev.off()
  
  return(sras)
  
  
}

future_cat_viz_ppl <- function(ssp = "585", th = full_th_k, diff_ras) {
  
  # add human impact to climate projections (difference in historical)
  # and compare to the suitability projected by the full historical models
  
  names <- c("Alfalfa 2021-2040", "Alfalfa 2041-2060", "Alfalfa 2061-2080", "Alfalfa 2081-2100", 
             "Corn 2021-2040", "Corn 2041-2060", "Corn 2061-2080", "Corn 2081-2100", 
             "Cotton 2021-2040", "Cotton 2041-2060", "Cotton 2061-2080", "Cotton 2081-2100", 
             "Hay 2021-2040", "Hay 2041-2060", "Hay 2061-2080", "Hay 2081-2100", 
             "Soy 2021-2040", "Soy 2041-2060", "Soy 2061-2080", "Soy 2081-2100", 
             "Wheat 2021-2040", "Wheat 2041-2060", "Wheat 2061-2080", "Wheat 2081-2100") 
  
  all_mean <- stack(Sys.glob(paste0("./out/future_ens/mean_stack_", ssp, "*.tif")))
  names(all_mean) <- names
  
  hs <- diff_ras
  all_diff <- stack(stack(replicate(4, hs[[1]])), stack(replicate(4, hs[[2]])), stack(replicate(4, hs[[3]])), stack(replicate(4, hs[[4]])), 
                    stack(replicate(4, hs[[5]])), stack(replicate(4, hs[[6]])))
  names(all_diff) <- names
  
  all_mean <- all_mean + all_diff # this adds the human influence surface for each crop
  
  hs <- stack(Sys.glob("./out/hist_grids/prob_*_FULL_DayMet.tif")) # try WC and DayMet
  hs <- dropLayer(hs, "prob_SWHEAT_FULL_DayMet")
  
  hsf <- stack(stack(replicate(4, hs[[1]])), stack(replicate(4, hs[[2]])), stack(replicate(4, hs[[3]])), stack(replicate(4, hs[[4]])), 
               stack(replicate(4, hs[[5]])), stack(replicate(4, hs[[6]])))
  names(hsf) <- names
  
  sras <- stack()
  
  for (i in 1:dim(all_mean)[[3]]) {
    
    print(names[[i]])
    
    # pull crop specific th
    cn <- strsplit(names[[i]], " ")[[1]][1]
    ths <- th %>% filter(CCROP == cn) %>% pull(TH)
    
    hist <- hsf[[i]]
    fut <- all_mean[[i]]
    
    hist[hist >= ths] <- 1
    hist[hist < ths] <- 0
    
    fut[fut >= ths] <- 11
    fut[fut < ths] <- 0
    
    # 11-1 = 10 remains suitable
    # 11-0 = 11 becomes suitable
    # 0-1 = -1 becomes unsuitable
    # 0-0 = 0 remains unsuitable
    
    r <- fut - hist
    r <- ratify(r)
    # rat <- levels(r)[[1]]
    rat <- data.frame(ID = c(-1, 0, 10, 11))
    
    rat$TYPE <- NA
    rat$TYPE[rat$ID == -1] <- "Becomes unsuitable" #red
    rat$TYPE[rat$ID == 0] <- "Remains unsuitable" # gray
    rat$TYPE[rat$ID == 10] <- "Remains suitable" # green
    rat$TYPE[rat$ID == 11] <- "Becomes suitable" # blue
    
    levels(r) <- rat
    sras <- stack(sras, r)
    
  }
  
  names(sras) <- names
  
  cols <- c("orange", "light gray", "#228B22", "purple")
  
  state_sp <- as(st_transform(states, sras@crs), "Spatial")
  
  fcat_viz <- levelplot(sras, 
                        margin = F,
                        xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                        ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                        scales = list(draw=F), 
                        layout = c(4, 6), 
                        names.attr = rep("", dim(sras)[3]),
                        colorkey = T,
                        col.regions = cols) +
    layer(sp.polygons(state_sp, lwd = 0.7, col = "black"))
  
  pdf(paste0("./figs/cat_fut_proj_", ssp, "_FULL_CSTH.pdf"))
  print(fcat_viz)
  dev.off()
  
  return(sras)
  
  
}

future_bin_viz <- function(ssp = "585", thb = thb, thf=thf, diff_ras) {
  
  # it's only in comparison to historical baseline that makes sense... 
  
  names <- c("Alfalfa 2021-2040", "Alfalfa 2041-2060", "Alfalfa 2061-2080", "Alfalfa 2081-2100", 
             "Corn 2021-2040", "Corn 2041-2060", "Corn 2061-2080", "Corn 2081-2100", 
             "Cotton 2021-2040", "Cotton 2041-2060", "Cotton 2061-2080", "Cotton 2081-2100", 
             "Hay 2021-2040", "Hay 2041-2060", "Hay 2061-2080", "Hay 2081-2100", 
             "Soy 2021-2040", "Soy 2041-2060", "Soy 2061-2080", "Soy 2081-2100", 
             "Wheat 2021-2040", "Wheat 2041-2060", "Wheat 2061-2080", "Wheat 2081-2100") 
  
  all_mean <- stack(Sys.glob(paste0("./out/future_ens/mean_stack_", ssp, "*.tif")))
  names(all_mean) <- names
  
  hs <- diff_ras
  all_diff <- stack(stack(replicate(4, hs[[1]])), stack(replicate(4, hs[[2]])), stack(replicate(4, hs[[3]])), stack(replicate(4, hs[[4]])), 
                    stack(replicate(4, hs[[5]])), stack(replicate(4, hs[[6]])))
  names(all_diff) <- names
  
  all_mean_full <- all_mean + all_diff # this adds the human influence surface for each crop
  
  all_mean_stack <- stack()
  all_mean_full_stack <- stack()
  
  for (i in 1:dim(all_mean_full)[[3]]) {
    
    print(names[[i]])
    
    # pull crop specific th
    cn <- strsplit(names[[i]], " ")[[1]][1]
    thfull <- thf %>% filter(CCROP == cn) %>% pull(TH)
    thbio <- thb %>% filter(CCROP == cn) %>% pull(TH)
    
    all_mean[[i]][all_mean[[i]] >= thbio] <- 1
    all_mean[[i]][all_mean[[i]] < thbio] <- 0
    
    all_mean_full[[i]][all_mean_full[[i]] >= thfull] <- 1
    all_mean_full[[i]][all_mean_full[[i]] <- thfull] <- 0
    
    all_mean_stack <- stack(all_mean_stack, all_mean[[i]])
    all_mean_full_stack <- stack(all_mean_full_stack, all_mean_full[[i]])
  }
  
  names(all_mean_stack) <- names
  names(all_mean_full_stack) <- names
    
  cols <- c("white", "#77815C")
  
  state_sp <- spTransform(state_sp, all_mean@crs)

  fcat_viz <- levelplot(all_mean_stack, 
                        margin = F,
                        xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                        ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                        scales = list(draw=F), 
                        layout = c(4, 6), 
                        names.attr = rep("", dim(all_mean)[3]),
                        colorkey = F,
                        col.regions = cols) +
    layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  pdf(paste0("./figs/bin_fut_proj_", ssp, "_WC.pdf"))
  print(fcat_viz)
  dev.off()
  
  fcat_viz <- levelplot(all_mean_full_stack, 
                        margin = F,
                        xlab = list(c("2021-2040", "2041-2060", "2061-2080", "2081-2100"), space = "top"),  # winter wheat
                        ylab=list(c("Wheat", "Soy", "Hay", "Cotton", "Corn", "Alfalfa"), rot=0), 
                        scales = list(draw=F), 
                        layout = c(4, 6), 
                        names.attr = rep("", dim(all_mean_full)[3]),
                        colorkey = F,
                        col.regions = cols) +
    layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
  
  pdf(paste0("./figs/bin_fut_proj_", ssp, "_WC_FULL.pdf"))
  print(fcat_viz)
  dev.off()
  
  
}

extract_suit <- function(ras_stack, ssp) {
  
  r2e <- extract(ras_stack, frr_prj, df=T)
  colnames(r2e) <- c("ID", 
                     "Alfalfa_2000-2020", "Corn_2000-2020", "Cotton_2000-2020", "Hay_2000-2020", "Soy_2000-2020", "Winter.wheat_2000-2020",
                     "Alfalfa_2021-2040", "Corn_2021-2040", "Cotton_2021-2040", "Hay_2021-2040", "Soy_2021-2040", "Winter.wheat_2021-2040",
                     "Alfalfa_2041-2060", "Corn_2041-2060", "Cotton_2041-2060", "Hay_2041-2060", "Soy_2041-2060", "Winter.wheat_2041-2060",
                     "Alfalfa_2061-2080", "Corn_2061-2080", "Cotton_2061-2080", "Hay_2061-2080", "Soy_2061-2080", "Winter.wheat_2061-2080",
                     "Alfalfa_2081-2100", "Corn_2081-2100", "Cotton_2081-2100", "Hay_2081-2100", "Soy_2081-2100", "Winter.wheat_2081-2100")
  
  out <- pivot_longer(r2e, cols = 2:ncol(r2e), names_to = "FULL", values_to = "SUITABILITY") 
  out <- out %>%
    separate(FULL, into = c("CROP", "YEAR"), sep = "_")
  
  out$SSP <- ssp
  
  saveRDS(out, paste0("./out/future_suit/ssp_", ssp, ".RDS"))
  
}











########################################################################################################################################
# Archive
########################################################################################################################################

# build_state_bricks <- function(cores, state_names) {
#   # this function constructs the clean state-bricks from the extent clips in the /States folder
#   foreach(i=1:length(state_names), 
#           
#           .packages= c("raster","foreach", "sf", "tidyverse"),
#           .export = c("state_names", "states")) %dopar% {
#             
#             dir.create(file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i]))) 
#             rasterOptions(tmpdir=file.path(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i])))  
#             
#             sras_files <- ifelse(state_names[i] == "Virginia", 
#                                  paste0("C:/Users/eburchf/Desktop/Data/CDL/States/*_Virginia.grd"),
#                                  paste0("C:/Users/eburchf/Desktop/Data/CDL/States/*", state_names[i], ".grd"))
#             sras_files <- Sys.glob(sras_files)
#             sras_stack <- stack(sras_files) # stack for multiple files, brick for single file
#             state_shp <- states %>% filter(STATE_NAME == state_names[i])
#             
#             # crop rasters
#             cropped <- raster::crop(sras_stack, state_shp) # clip to extent
#             masked <- raster::mask(cropped, state_shp) # include only pixels in shapefile boundaries
#             writeRaster(masked, filename = paste0("C:/Users/eburchf/Desktop/Data/CDL/State_stacks/", state_names[i], ".gri"),
#                         overwrite = T)
#             remove(cropped, masked)
#             unlink(paste0("C:/Users/eburchf/Desktop/TRASH/",state_names[i]), recursive = TRUE, force = TRUE) # delete temporary files
#           }
# }


# build_error_rasters <- function(crop_error, model) {
#   
#   error_rasters <- list()
#   error_shp <- list()
#   
#   #crop_names <- crop_names[crop_names != "SWHEAT"]
#   
#   for (i in 1:length(crop_names)) {
#     
#     crop <- crop_names[[i]]
#     print(crop)
#     
#     varname1 = paste0(crop, "_PRED")
#     varname2 = paste0(crop, "_ACT")
#     
#     crop_error$AS <- ifelse(crop_error[[varname1]] == crop & crop_error[[varname2]] == crop, 1, 0) # correct classification of suitability
#     crop_error$FS <- ifelse(crop_error[[varname1]] == crop & crop_error[[varname2]] != crop, 1, 0) # falsely classified suitability
#     crop_error$MS <- ifelse(crop_error[[varname1]] != crop & crop_error[[varname2]] == crop, 1, 0) # falsely classified unsuitability
#     crop_error$AUS <- ifelse(crop_error[[varname1]] != crop & crop_error[[varname2]] != crop, 1, 0) # correctly classified unsuitability
#     
#     crop_error$CATS <- ifelse(crop_error$AS == 1, 1, 0) # green
#     crop_error$CATS <- ifelse(crop_error$FS == 1, 2, crop_error$CATS) # red
#     crop_error$CATS <- ifelse(crop_error$MS == 1, 3, crop_error$CATS) # blue
#     crop_error$CATS <- ifelse(crop_error$AUS == 1, 4, crop_error$CATS) # gray
#     
#     error_ras <- rasterize(crop_error[,"CATS"], hs[[1]], fun = modal, background = 0)
#     error_rasters[[crop]] <- error_ras[[2]]
#     writeRaster(error_ras[[2]], paste0("./out/hist_error/", crop, model, ".tif"), overwrite=T)
#     
#     error_shp[[crop]] <- crop_error
#     
#     # total <- rasterize(st_coordinates(crop_error), hs[[1]], fun = "count", background = 0) # total points in each grid cell
#     # total[total == 0] <- NA
#     # total_correct <- rasterize(st_coordinates(crop_error %>% filter(AS == 1 | AUS == 1)), hs[[1]], fun = "count", background = 0)
#     # total_rasters[[crop]] <- total_correct/total
#     # writeRaster(perc_correct, paste0("./out/hist_error/pc_", crop, model, ".tif"), overwrite=T)
#     
#   }
#   
#   return(error_shp)
#   
# }
# 
# 

# build_future_grids <- function(out_ext = "_CLIM.RDS", full_list) {
#   
#   attr <- brick("./out/attr_grids/future_attr.grd") # just resampled to new proj
#   fut <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_grids/*.grd")) # cropped, masked and named
#   
#   for (j in 1:length(fut)) {
#     
#     fr <- brick(fut[[j]])
#     fn <- strsplit(fut[[j]], "/")[[1]][9]
#     fn <- substr(fn, 1, nchar(fn)-4)
#     all_attr <- stack(fr, attr)
#     print(fn)
#     
#     for (i in 1:length(crop_names)) {
#       
#       crop_name <- crop_names[[i]]
#       print(crop_name)
#       rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext)) 
#       apsub <- subset(all_attr, names(full_list[[crop_name]]))
#       prob <- predict(apsub, rf, type = "prob", index = 2) # prob of value 1, col 2
#       writeRaster(prob, paste0("./out/future_grids/prob_", fn, "_", crop_name, substr(out_ext, 1, nchar(out_ext)-4), ".tif"), overwrite=T)
#       
#     }
#     
#   }
#   
# }

# build_projection_figure <- function(scenario, model) {
#   
#   m1 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2041-2060_*BIO.tif"))
#   m2 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2061-2080_*CLIM.tif"))
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2081-2100_*CLIM.tif"))
#   
#   allvis <- stack(hs, m1, m2, m3) # hs made outside function
#   # allvis[allvis == 0] <- NA # for final figure, takes a second
#   coul2 <- inferno(25)
#   
#   future_clim <- levelplot(allvis, 
#                            margin = F,
#                            xlab = list(names(hs), space = "top"),  # winter wheat
#                            ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2000-2020"), rot=0), 
#                            scales = list(draw=F), 
#                            layout = c(6, 4), 
#                            names.attr = rep("", dim(allvis)[3]),
#                            col.regions = coul2) + # inferno viridis
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   pdf(paste0("./figs/", model, "_", scenario, ".pdf"))
#   print(future_clim)
#   dev.off()
#   
#   return(future_clim)
#   
# }
# 
# 
# 
# build_hist_vis_roi <- function(out_ext, state_abbr) {
#   
#   all_ras <- list()
#   library(viridis)
#   coul2 <- viridis::inferno(25)
#   
#   for (i in 1:length(crop_names)) {
#     r <- raster(paste0("./out/hist_grids/prob_", crop_names[[i]], substr(out_ext, 1, nchar(out_ext)-4), ".tif"))
#     all_ras[[i]] <- r
#   }
#   
#   all_ras <- stack(all_ras[[1]], all_ras[[2]], all_ras[[3]], all_ras[[5]], all_ras[[6]], all_ras[[7]])
#   names(all_ras) <- clean_crop_names[!(clean_crop_names %in% "Spring wheat")]
#   
#   state_sp <- state_sp %>% filter(STATE_ABBR == state_abbr)
#   all_ras <- raster::crop(all_ras, state_sp)
#   
#   hist_clim <- rasterVis::levelplot(all_ras, 
#                                     margin = F,
#                                     xlab = "",
#                                     ylab = "",
#                                     # xlab = list(names(all_ras), space = "top"),  # winter wheat
#                                     # ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2000-2020"), rot=0), 
#                                     scales = list(draw=F), 
#                                     layout = c(3, 2), 
#                                     names.attr = names(all_ras), #rep("", dim(allvis)[3]),
#                                     col.regions = coul2) + # inferno viridis
#     latticeExtra::layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   pdf(paste0("./figs/hist_", substr(out_ext, 1, nchar(out_ext)-4), ".pdf"))
#   print(hist_clim)
#   dev.off()
#   out <- list(hist_clim, all_ras)
#   
#   return(out)
#   
#   
# }
# 
# build_projection_figure_roi <- function(scenario, model, state_abbr) {
#   
#   m1 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2041-2060_*BIOCLIM.tif"))
#   m2 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2061-2080_*BIOCLIM.tif"))
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/prob_wc2.1_2.5m_bioc_", model, "_ssp", scenario, "_2081-2100_*BIOCLIM.tif"))
#   
#   allvis <- stack(hs, m1, m2, m3) # hs made outside function
#   # allvis[allvis == 0] <- NA # for final figure, takes a second
#   coul2 <- inferno(25)
#   
#   state_sp <- state_sp %>% filter(STATE_ABBR == state_abbr)
#   allvis <- raster::crop(allvis, state_sp)
#   
#   future_clim <- levelplot(allvis, 
#                            margin = F,
#                            xlab = list(names(hs), space = "top"),  # winter wheat
#                            ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2000-2020"), rot=0), 
#                            scales = list(draw=F), 
#                            layout = c(6, 4), 
#                            names.attr = rep("", dim(allvis)[3]),
#                            col.regions = coul2) + # inferno viridis
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   #pdf(paste0("./figs/", model, "_", scenario, ".pdf"))
#   #print(future_clim)
#   #dev.off()
#   
#   return(future_clim)
#   
# }

# build_ensembles <- function() {
#   
#   # ended up not using this
#   
#   ext_list <- c("ssp245_2021-2040", "ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
#                 "ssp370_2021-2040", "ssp370_2041-2060", "ssp370_2061-2080", "ssp370_2081-2100",
#                 "ssp585_2021-2040", "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")
#   
#   for (i in 1:length(ext_list)) {
#     
#     ext <- ext_list[[i]] # will need list of these
#     print(ext)
#     
#     rlist <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_grids/*", ext, ".grd"))
#     
#     n <- length(rlist)
#     print(n)
#     
#     outras <- stack(rlist[[1]])
#     for (i in 2:length(n)) {
#       
#       r <- stack(rlist[[n]])
#       outras <- outras + r
#       
#     }
#     
#     outras <- outras/n 
#     
#     writeRaster(outras, paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_ens/", ext, ".grd"), overwrite=T)
#     
#     
#   }
# }
# 
# build_future_grids_ens <- function(out_ext = "_CLIM.RDS", full_list) {
#   
#   attr <- brick("./out/attr_grids/5k_hist_attr_WC.grd") # just resampled to new proj
#   attr[["T_SAND"]] <- attr[["TS_SAND"]]
#   attr <- subset(attr, c(varlist[["Agricultural"]], topo, soil_sub))
#   fut <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/Data/WorldClim/future/CS_ens/*.grd")) # cropped, masked and named
#   
#   for (j in 1:length(fut)) {
#     
#     fr <- brick(fut[[j]])
#     names(fr) <- climate_full
#     fn <- strsplit(fut[[j]], "/")[[1]][9]
#     fn <- substr(fn, 1, nchar(fn)-4)
#     all_attr <- stack(fr, attr)
#     print(fn)
#     
#     for (i in 1:length(crop_names)) {
#       
#       crop_name <- crop_names[[i]]
#       print(crop_name)
#       rf <- readRDS(paste0("./out/RF_historical/", crop_name, out_ext)) 
#       apsub <- subset(all_attr, varlist[["Biophysical"]])
#       prob <- predict(apsub, rf, type = "prob", index = 2) # prob of value 1, col 2
#       writeRaster(prob, paste0("./out/future_grids/ens_prob_", fn, "_", crop_name, substr(out_ext, 1, nchar(out_ext)-4), "_WC.tif"), overwrite=T)
#       
#     }
#     
#   }
#   
#   
# }
# 
# build_ensemble_figure <- function(ssp) {
#   
#   hs <- stack("./out/hist_grids/prob_BIO_DayMet_RS.tif")
#   
#   state_sp <- as(states, "Spatial")
#   state_sp <- spTransform(state_sp, hs@crs)
#   
#   m0 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2021-2040_*BIO_WC.tif"))
#   m1 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2041-2060_*BIO_WC.tif"))
#   m2 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2061-2080_*BIO_WC.tif"))
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2081-2100_*BIO_WC.tif"))
#   
#   allvis <- stack(m0, m1, m2, m3) # hs made outside function
#   coul2 <- inferno(20)
#   
#   future_clim <- levelplot(allvis, 
#                            margin = F,
#                            xlab = list(names(m0), space = "top"),  # winter wheat
#                            ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2021-2040"), rot=0), 
#                            scales = list(draw=F), 
#                            layout = c(6, 4), 
#                            names.attr = rep("", dim(allvis)[3]),
#                            # col.regions = coul2,
#                            col.regions = terrain.colors(25),
#                            at = seq(0, 1, length.out = 20)) + # inferno viridis
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   pdf(paste0("./figs/future_probs_", ssp, ".pdf"))
#   print(future_clim)
#   dev.off()
#   
#   out <- list(future_clim, allvis)
#   return(out)
#   
#   
# }
# make_mf <- function(vname) {
#   
#   pv <- colnames(census)[vname]
#   
#   fm <- as.formula(paste("AP ~ ", pv))
#   mod <- glm(fm, data = census, family = binomial)
#   
#   r2 <- with(summary(mod), 1-deviance/null.deviance)
#   coef <- round(mod$coefficients[2], 2)
#   
#   q <- quantile(census[,pv], c(0.02, 0.98))
#   by = (q[[2]] - q[[1]])/1000
#   pred <- data.frame(out = seq(from = q[[1]], to = q[[2]], 
#                                by = by))
#   colnames(pred) <- pv
#   pred$AP <- predict(mod, newdata = pred, type = "response")
#   
#   ggplot(pred, aes(x = pred[,pv], y = AP)) +
#     geom_point(alpha = 0.2) +
#     #geom_line() +
#     #geom_smooth(se = F, method.args = (list(family = binomial))) +
#     project_theme +
#     xlab(pv) +
#     ggtitle(paste0(coef, ", R2: ", round(r2, 3))) +
#     ylim(c(0,1)) +
#     ylab("") +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank())
#   
# }
# 
# 

# build_ensemble_delta_figure <- function(ssp, out_ext = "_BIO") {
#   
#   hs <- stack(Sys.glob("./out/hist_grids/prob_*_BIO_WC.tif")) # probably need to go back to WC
#   
#   
#   if (out_ext == "_FULL") {
#     
#     hs <- stack(Sys.glob("./out/hist_grids/prob_*_FULL_WC.tif"))
#   }
#   
#   m0 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2021-2040_*BIO_WC.tif"))
#   m1 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2041-2060_*BIO_WC.tif"))
#   m2 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2061-2080_*BIO_WC.tif"))
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, "_2081-2100_*BIO_WC.tif"))
#   
#   allvis <- stack(m0-hs, m1-hs, m2-hs, m3-hs)
#   
#   state_sp <- as(states, "Spatial")
#   state_sp <- spTransform(state_sp, allvis@crs)
#   
#   coul2 <- colorRampPalette(brewer.pal(11, "BrBG"))
#   
#   future_clim <- levelplot(allvis, 
#                            margin = F,
#                            xlab = list(names(m0), space = "top"),  # winter wheat
#                            ylab=list(c("2081-2100", "2061-2080", "2041-2060", "2021-2040"), rot=0), # check order here...
#                            scales = list(draw=F), 
#                            layout = c(6, 4), 
#                            names.attr = rep("", dim(allvis)[3]),
#                            col.regions = coul2,
#                            at = seq(-0.9,  0.9, length.out = 11)) +
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   
#   g <- grid.arrange(future_clim, ncol=1, nrow=1)
#   ggsave(paste0("./figs/future_diff_", ssp, out_ext, ".png"), g, width=10, height=10)
# }


# future_cat_viz <- function(year_ext = "_2021-2040_", ssp = "585", th = 0.5) {
#   
#   # hs <- stack("./out/hist_grids/prob_BIO_DayMet_RS.tif")
#   hs <- stack(Sys.glob("./out/hist_grids/prob_*_BIO.tif"))
#   
#   state_sp <- as(states, "Spatial")
#   state_sp <- spTransform(state_sp, hs@crs)
#   
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, year_ext, "*BIO_WC.tif"))
#   m3 <- projectRaster(m3, hs, method = "bilinear")
#   
#   hs[hs >= th] <- 1 # suitable
#   hs[hs < th] <- 0 # unsuitable
#   
#   m3[m3 >=th] <- 11 # suitable
#   m3[m3 < th] <- 0 # unsuitable
#   
#   # 11-1 = 10 remains suitable
#   # 11-0 = 11 becomes suitable
#   # 0-1 = -1 becomes unsuitable
#   # 0-0 = 0 remains unsuitable
#   
#   sras <- stack()
#   
#   for (i in 1:dim(m3)[[3]]) {
#     
#     r <- m3[[i]] - hs[[i]]
#     r <- ratify(r)
#     # rat <- levels(r)[[1]]
#     rat <- data.frame(ID = c(-1, 0, 10, 11))
#     
#     rat$TYPE <- NA
#     rat$TYPE[rat$ID == -1] <- "Becomes unsuitable" #red
#     rat$TYPE[rat$ID == 0] <- "Remains unsuitable" # gray
#     rat$TYPE[rat$ID == 10] <- "Remains suitable" # green
#     rat$TYPE[rat$ID == 11] <- "Becomes suitable" # blue
#     
#     levels(r) <- rat
#     sras <- stack(sras, r)
#     
#     
#   }
#   names(sras) <- names(m3)
#   
#   cols <- c("red", "light gray", "#228B22", "blue")
#   state_sp <- spTransform(state_sp, sras@crs)
#   
#   fcat_viz <- levelplot(sras, 
#                         margin = F,
#                         xlab = "",
#                         ylab = "",
#                         names.attr = names(sras),
#                         colorkey=T,
#                         #  ylab=list(c("SSP585", "SSP245"), space = "top"),  # winter wheat
#                         # xlab=list(c("2000-2020", "2041-2060", "2061-2080", "2081-2100"), rot=0),
#                         scales = list(draw=F), 
#                         #  layout = c(4, 2), 
#                         # names.attr = rep("", dim(allvis)[3]),
#                         col.regions = cols) + # inferno viridis
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   
#   pdf(paste0("./figs/cat_fut_proj", year_ext, ssp, ".pdf"))
#   print(fcat_viz)
#   dev.off()
#   
#   return(sras)
#   
#   
# }
# 
# future_cat_viz_ppl <- function(year_ext = "_2021-2040_", ssp = "585", th = 0.5, diff_ras = all_ras) {
#   
#   hs <- stack("./out/hist_grids/prob_FULL_DayMet_RS.tif")
#   
#   state_sp <- as(states, "Spatial")
#   state_sp <- spTransform(state_sp, hs@crs)
#   
#   m3 <- build_future_stack(dir = paste0("./out/future_grids/ens_prob_ssp", ssp, year_ext, "*BIO_WC.tif"))
#   m3 <- m3 + all_ras # add people filter
#   
#   hs[hs >= th] <- 1
#   hs[hs < th] <- 0
#   
#   m3[m3 >=th] <- 11
#   m3[m3 < th] <- 0
#   
#   sras <- stack()
#   
#   for (i in 1:dim(m3)[[3]]) {
#     
#     r <- m3[[i]] - hs[[i]]
#     r <- ratify(r)
#     # rat <- levels(r)[[1]]
#     rat <- data.frame(ID = c(-1, 0, 10, 11))
#     
#     rat$TYPE <- NA
#     rat$TYPE[rat$ID == -1] <- "Becomes unsuitable" #red
#     rat$TYPE[rat$ID == 0] <- "Remains unsuitable" # gray
#     rat$TYPE[rat$ID == 10] <- "Remains suitable" # green
#     rat$TYPE[rat$ID == 11] <- "Becomes suitable" # blue
#     
#     levels(r) <- rat
#     sras <- stack(sras, r)
#     
#     
#   }
#   names(sras) <- names(hs)
#   
#   cols <- c("red", "light gray", "#228B22", "blue")
#   state_sp <- spTransform(state_sp, sras@crs)
#   
#   fcat_viz <- levelplot(sras, 
#                         margin = F,
#                         xlab = "",
#                         ylab = "",
#                         names.attr = names(sras),
#                         colorkey=F,
#                         #  ylab=list(c("SSP585", "SSP245"), space = "top"),  # winter wheat
#                         # xlab=list(c("2000-2020", "2041-2060", "2061-2080", "2081-2100"), rot=0),
#                         scales = list(draw=F), 
#                         #  layout = c(4, 2), 
#                         # names.attr = rep("", dim(allvis)[3]),
#                         col.regions = cols) + # inferno viridis
#     layer(sp.polygons(state_sp, lwd = 0.7, col = "#696969"))
#   
#   
#   pdf(paste0("./figs/cat_fut_proj", year_ext, ssp, "_HUM.pdf"))
#   print(fcat_viz)
#   dev.off()
#   
#   return(sras)
#   
#   
# }
