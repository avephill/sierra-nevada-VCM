# This file contains the functions used in sdm_main. This is essentially the
# static workhorse of the sdm. All changes for different runs should be set in sdm_main


## First we'll build an SDM on historic climate and weislander dist and then 
# project it forward

# For changing variable names between bio1 and MAT
bioclim.dic <- read.csv("Data/Dictionaries/bioclim_dictionary.csv")


# Here's a function I'm going to specify that preps occurrence data for feeding
# sdm. Let's get this out of the way here:
# This function takes occurrence lonlat dataframe: data.frame(x,y) where x = lon
# and y = lat and a stack of the predictor rasters and spits out the sdm stack 
# and the results from the vifstep function which weeds out collinear variables
sdmdataPrep <- function(pres_lonlat.df,
                        abs_lonlat.df,
                        predictor.stack, 
                        predictor_names){
  # These are the predictor values at locations of species presence
  presvals <- raster::extract(predictor.stack, pres_lonlat.df)
  
  # predictor values at random locations
  absvals <- raster::extract(predictor.stack, abs_lonlat.df)
  
  # We know that probability of presence is 1 for areas where the species
  # was found, and we assume it's 0 for the random background points
  Y <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
  sdmdata.df <- data.frame(cbind(Y, rbind(presvals, absvals))) %>% 
    dplyr::select(Y, all_of(predictor_names))

  sdmdata <- sdmData(Y ~ ., train = sdmdata.df)
  return(list(sdmdata = sdmdata))
}



checkCollinearity <- function(predictor.stack){
  
  predictor_studyarea.df <- values(predictor.stack) %>% 
    data.frame() %>% 
    na.omit()
  
  # Instead of looking for collinearity across the occurrences we're going to weed
  # out collinearity with regards to the whole study area
  # vif.out <- vifstep(dplyr::select(sdmdata.df, -Y), th = 10)
  vif.out <- vifstep(predictor_studyarea.df, th = 10)
  varnames <- names(predictor_studyarea.df)
  vif.list <- list()
  for(i in 1:50){
    vif <- vifstep(predictor_studyarea.df, th = 10)
    vif.list[[i]] <- vif@results
  }
  
  # This tells us which variables switch back and forth, i.e. collinear with each other
  vif.namecount <- table(unlist(lapply(vif.list, function(x)x[,'Variables'])))
  
  # This tells us which variables are the least collinear on avg
  vif.df <- data.frame(Variables=NA, MeanVIF=NA)
  for (i in 1:length(varnames)){
    env.name <- varnames[i]
    vif_avg <- mean(unlist(lapply(vif.list,function(x)x[x$Variables==env.name,'VIF'])))
    vif.df <- rbind(vif.df,c(env.name,vif_avg))
  }
  
  vif.tab <- vif.df[order(vif.df$MeanVIF),] %>% na.omit %>% filter(MeanVIF != "NaN")
  names(vif.tab) <- c('Climate Variable',"Mean VIF")
  
  collinearity_results <- list(vif_avg.df = vif.tab, vif_namecount.tab = vif.namecount)
  return(collinearity_results)
}



# Make SDM from historical data (weislander) with historical climate (prism)
# and project them to current time using current climate (prism)
# Then produce plots that show the difference between present prediction and past
# prediction
buildSDM <- function(tree.grid,
                     repMethod = NA, train_N = "all",
                     run.name, n_replicates = 0,
                     sdm_method,
                     predictor_names,
                     predictor.stack,
                     curr_predict.stack,
                     removeCollinearity = T){
  
  # First let's make the directory where all our results will go
  results.dir <- paste0("Results/SDM/",run.name,
                        "_repMethod", repMethod,
                        "_trn",train_N,"_boot",n_replicates,"/")
  dir.create(results.dir)
  
  
  
  weis_tree.pts <- data.frame(rasterToPoints(tree.grid, 
                                  fun=function(x){x == 1})[,1:2])
  
  weis_not_tree.pts <- data.frame(rasterToPoints(tree.grid, 
                                      fun=function(x){x == 0})[,1:2])
  
  # Use all points unless a smaller N is set
  if(train_N == "all"){
    # These points will be used to train the model
    allPosTrain.N <- nrow(weis_tree.pts) #* (1 - test_train.ratio)
    pres_weisTree.df <- data.frame(
      weis_tree.pts[sample(nrow(weis_tree.pts), allPosTrain.N), ])
    
    allAbsTrain.N <- nrow(weis_not_tree.pts) #* (1 - test_train.ratio)
    abs_weisTree.df <- data.frame(
      weis_not_tree.pts[sample(nrow(weis_not_tree.pts), allAbsTrain.N), ])
    
  } else {
    pres_weisTree.df <- data.frame(
      weis_tree.pts[sample(nrow(weis_tree.pts), train_N), ])
    
    abs_weisTree.df <- data.frame(
      weis_not_tree.pts[sample(nrow(weis_not_tree.pts), train_N), ])
    print("if this throws an error then decrease N")
    
  }
  
  print(paste0("training presence points: ", nrow(pres_weisTree.df)))
  print(paste0("training absence points: ", nrow(abs_weisTree.df)))
  
  
  weisTree_dataprep.results <- sdmdataPrep(pres_lonlat.df = pres_weisTree.df, 
                                           abs_lonlat.df = abs_weisTree.df,
                                           predictor.stack = predictor.stack,
                                           predictor_names = predictor_names)
  
  weisTree_sdmdata <- weisTree_dataprep.results$sdmdata
  
  
  
  # ENSEMBLE let's make an ensemble model of all of them
  # a number of methods can be used, see documentation, but they include 
  # weighted mean, unweighted mean, median, entropy, etc
  
  # We can't run sdm with x and y in sdmdata so we need to remove that from the formula
  predictor.names <- weisTree_sdmdata@features.name
  weisTree_sdmFormula = reformulate(predictor.names, response="Y")
  print("Training SDM now...")
  # If no n_replicates is specified let's not do any bootstraps
  if(n_replicates > 1){
    if(repMethod == "cv_spblock"){
      # Spatial Blocking of CV
      
      pre.sf <- rbind(pres_weisTree.df %>% mutate(Y = 1), 
                      abs_weisTree.df %>% mutate(Y = 0))
      block_sp.sf <- st_as_sf(pre.sf, coords = c("x", "y"), crs = crs(predictor.stack))
      
      # It's important that pre.sf and block_sp.sf are in same order
      sp_extract.df <- raster::extract(predictor.stack, pre.sf %>% select(-Y), df = T) %>% 
        mutate(Y = pre.sf$Y)
      
      # rangeExplorer(rasterLayer = predictor.stack,
      #               species = "Y",
      #               speciesData = block_sp.sf,
      #               rangeTable = NULL,
      #               minRange = 30000, # limit the search domain
      #               maxRange = 100000)
      
      # sac <- spatialAutoRange(rasterLayer = predictor.stack[[pred_names]],
      #                         sampleNumber = 5000,
      #                         doParallel = TRUE,
      #                         showPlots = TRUE)
      # Says 148km block size
      
      # spatial blocking by rows and columns with checkerboard assignment
      sb <- spatialBlock(speciesData = block_sp.sf,
                         species = "Y",
                         rasterLayer = predictor.stack[[pred_names]],
                         theRange = 70000, # size of the blocks
                         k = n_replicates,
                         selection = "random",
                         iteration = 100, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, # shift the blocks horizontally
                         yOffset = 0)
      
      folds <- sb$foldID
      
      # create an empty vector to store the AUC of each fold
      models = list()
      for(k in seq_len(5)){
        # extracting the training and testing indices
        # this way only works with foldID
        trainSet <- which(folds != k) # training set indices
        testSet <- which(folds == k) # testing set indices
  
        # Leaving ID for later sanity check
        sdm.data <-  sdmData(Y ~ ., 
                             train = sp_extract.df[trainSet,] %>% select(-ID),
                             test = sp_extract.df[testSet,] %>% select(-ID))
        
        print("about to actually train spatial block cv")
        KweisTree_sdm <- sdm::sdm(weisTree_sdmFormula, 
                                  data = sdm.data, 
                                  methods = sdm_method#, 
                                  #parallelSettings = list(ncore = 6)
                                  )
        
        
        models <- c(models, KweisTree_sdm)
        

        
      }

      
      weisTree_sdm <- models
      
      
    } else {
      # Regular cv
      weisTree_sdm <- sdm::sdm(weisTree_sdmFormula, data = weisTree_sdmdata, 
                               methods=sdm_method, 
                               replication=c(repMethod), cv.folds = n_replicates,
                               parallelSettings = list(ncore = 10))
    }
    
    
  } else {
    # No replicates
    weisTree_sdm <- sdm::sdm(weisTree_sdmFormula, data = weisTree_sdmdata, 
                             methods=sdm_method,
                             parallelSettings = list(ncore = 10))
  sd}
  
  
  #### NOW WE ENSEMBLE
  
  # This ensemble is a prediction, it is no longer an SDM anymore
  ensembleFile <- paste0(results.dir,"weis", run.name, "_historic.tif")
  if(file.exists(ensembleFile)) file.remove(ensembleFile)


  # If there are bootstraps or multiple algorithms used, make an ensemble of 
  # the predictions, otherwise just use the single prediction
  if(n_replicates < 1){
    
    weisTree_ensemble <- predict(weisTree_sdm, predictor.stack)
    crs(weisTree_ensemble) <- crs(predictor.stack)
    writeRaster(weisTree_ensemble, filename = ensembleFile)
    
  } else if(repMethod == "cv_spblock"){
    
    weisTree_ensemble <- diyEnsemble(sdm_list = weisTree_sdm, 
                                     predictor.stack = predictor.stack, 
                                     write.fp = ensembleFile)
    
  } else {

    weisTree_ensemble <- sdm::ensemble(weisTree_sdm, predictor.stack, 
                                       setting=list(method="weighted", stat="AUC"),
                                       filename = ensembleFile)
  }
  
  
  
  
  ## PROJECT TO PRESENT DAY ###
  
  
  # We're not overlaying current distribution quite yet
  # calveg.grid <- raster("Data/EVT/Rasters/CalvegEVT.grd")
  # calveg_tree.grid <- calveg.grid[calveg.grid %in% whr_tree.ints]
  # calveg_tree.grid[!(calveg_tree.grid %in% whr_tree.ints)] <- NA
  
  # This ensemble is a prediction, it is no longer an SDM anymore
  ensembleFile <- paste0(results.dir,"weis",run.name,"_to_Present.tif")
  if(file.exists(ensembleFile)) file.remove(ensembleFile)
  # If there are bootstraps or multiple algorithms used, make an ensemble of 
  # the predictions, otherwise just use the single prediction
  if(n_replicates < 1){
    weisTree_sdm_current.prediction <- predict(weisTree_sdm, curr_predict.stack)
    crs(weisTree_sdm_current.prediction) <- crs(curr_predict.stack)
    writeRaster(weisTree_sdm_current.prediction, filename = ensembleFile)
    
  } else if(repMethod == "cv_spblock"){
    
    weisTree_sdm_current.prediction <- diyEnsemble(sdm_list = weisTree_sdm, 
                                     predictor.stack = curr_predict.stack, 
                                     write.fp = ensembleFile)
    
  } else {
    weisTree_sdm_current.prediction <- sdm::ensemble(weisTree_sdm, curr_predict.stack, 
                                                     setting=list(method="weighted", stat="AUC"),
                                                     filename = ensembleFile)
  }
  
  
  
  # Let's plot the SDM ensemble prediction
  hist.p <- gplot(weisTree_ensemble) + geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'blue', high = 'green', limits=c(0, 1), breaks=seq(0,1,by=0.2)) +
    coord_equal()
  
  curr.p <- gplot(weisTree_sdm_current.prediction) + geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'blue', high = 'green', limits=c(0, 1), breaks=seq(0,1,by=0.2)) +
    coord_equal()
  
  hist_pres_pred.plot <- ggarrange(hist.p, curr.p, 
                                   labels = c("Historic Prediction", "Current Prediction"),
                                   ncol = 2, nrow = 1, common.legend = T)
  
  
  par(mfrow=c(1,1))
  weisTree_dif.grid <- weisTree_sdm_current.prediction - weisTree_ensemble
  weisTree_dif.fn <- paste0(results.dir, "weis",run.name,"_habDif.tif")
  writeRaster(weisTree_dif.grid, filename = weisTree_dif.fn, overwrite = T)
  
  
  # Data prep for geo_raster() and geom_sf()
  weisTree_dif.df <- weisTree_dif.grid %>% 
    as("SpatialPixelsDataFrame") %>%
    as.data.frame() %>% rename(value = 'layer')
  alt.contour <- st_as_sf(
    rasterToContour(alt, levels = seq(0, 5000, 1000)))
  ca.sf <- st_as_sf(ca.shp)
  
  delta_suit.plot <- ggplot() + 
    geom_raster(data = weisTree_dif.df, aes(x,y, fill = value, alpha = 1)) +
    scale_fill_gradient2(name = paste("Δ ",run.name, " Suitability"), 
                         low = 'red', mid = 'white', high = 'blue', midpoint = 0) +
    coord_equal() +
    geom_sf(data = alt.contour, aes(color = level) , size = .75) +
    scale_colour_brewer(name = "Alt Band", type = "seq",
                        palette = "Blues", direction = -1) +
    geom_sf(data = ca.sf, aes(), size = 1, fill = NA) +
    xlim(c(min(weisTree_dif.df$x),max(weisTree_dif.df$x))) +
    ylim(c(min(weisTree_dif.df$y),max(weisTree_dif.df$y))) +
    guides(alpha=FALSE) + theme_minimal() + 
    ggtitle(paste0(run.name, " Δ Habitat Suitability (1930-2015)"))
  
  sdm_results <- list(delta_suit.plot = delta_suit.plot,
                      hist_pres_pred.plot = hist_pres_pred.plot,
                      weis.sdm = weisTree_sdm,
                      sdm_dataprep.results = weisTree_dataprep.results)
  save(sdm_results, file=paste0(results.dir,"sdm_results.RData"))
  return(sdm_results)
}


# This function ensembles lists of sdms because i can't replicate an SDM object
# after doing spatial blocking
diyEnsemble <- function(sdm_list, predictor.stack, write.fp, SDM_DIR = sdm_dir){
  # We're going to do a weighted mean ourself. 
  # after inspecting the ensembling source code
  # getMethod("ensemble", signature = c(x="sdmModels", newdata="Raster")) %>% View()
  # it looks like they just multiply the predicted raster by the AUC or whatever metric,
  # then average
  k_performance.df = data.frame()
  sdmstats <- c('AUC', 'COR', 'Deviance')
  
  hist_hsm.stack <- raster::stack()
  for (i in 1:length(sdm_list)){
    k.sdm <- sdm_list[[i]]
    
    # If there are more than 1 method in each fold, then we'll ensemble the 
    # methods within the fold and then weighted average the folds
    if(length(getModelId(k.sdm)) > 1){
      weisTree_pred_i <- sdm::ensemble(k.sdm, predictor.stack, 
                                       setting=list(method="weighted", stat="AUC"),
                                       filename = paste0(SDM_DIR, "junk"), overwrite=T)
      
      # Extract model performance from each sdm item, going to average across 
      # sdm methods
      model_eval <- getEvaluation(k.sdm, stat = sdmstats) %>%
        select(-modelID) %>% 
        summarise(across(everything(), mean))
      
    } else {
      weisTree_pred_i <- predict(k.sdm, predictor.stack, nc = 10)
      
      # Extract model performance from each sdm item
      model_eval <- getEvaluation(k.sdm, stat = sdmstats) %>%
        select(-modelID)
      
    }

    crs(weisTree_pred_i) <- crs(predictor.stack)
    hist_hsm.stack <- addLayer(hist_hsm.stack, weisTree_pred_i)
    
    k_performance.df <- rbind(k_performance.df, model_eval)
    
  }
  
  
  
  # shift AUCs up such that the highest AUC is 1 and lower AUCs have weights below 1
  k_weights <- k_performance.df %>% 
    # mutate(weights = AUC + (1 - max(AUC))) %>% 
    # mutate(weights = AUC * 10) %>% 
    # pull(weights) 
    pull(AUC)
  
  sdm_ensemble <- weighted.mean(hist_hsm.stack, w = k_weights)
  # weisTree_ensemble <- mean(hist_hsm.stack)
  writeRaster(sdm_ensemble, filename = write.fp)
  return(sdm_ensemble)
}

# extract and format evaluation metrics from ragtag SDM list
evalSDMlist <- function(sdm_list, sdmstats){
  sdm_eval.df <- data.frame()
  for(k in 1:length(sdm_list)){
    
    k.sdm <- sdm_list[[k]]
    k_sdm_eval.df <- getEvaluation(k.sdm, stat = sdmstats) %>%
      tibble() %>% 
      rename(AUCtest = AUC) %>% 
      mutate(modelBlockID = k,
             AUCtrain = getEvaluation(k.sdm, stat = "AUC", wtest = "training") %>% 
               pull("AUC")) %>% 
      mutate(`AUCtest - AUCtrain` = AUCtest - AUCtrain)
    
    sdm_eval.df  <- rbind(sdm_eval.df, k_sdm_eval.df)
  }
  
  return(sdm_eval.df)
}

climImportSDMlist <- function(sdm_list){
  sdm_varimp.df <- data.frame()
  
  method_vec <- sdm_list[[1]]@run.info$method %>% as.character()
  
  method.dic <- tibble(method = method_vec ,
                      slotNum = 1:length(method_vec ))
  # For each item in sdm list ()
  for(k in 1:length(sdm_list)){
    
    k.sdm <- sdm_list[[k]]
    k.varimp <- getVarImp(k.sdm)
    k.varimp.list <- k.varimp@varImportanceList
    slotNums <- names(k.varimp.list)
    k.df.list <- lapply(k.varimp.list, FUN = function(x){
      k.method.varimp <- x@varImportance %>%
        data.frame() %>%
        tibble()
    }) 
    
    k.varimp.df <-  do.call("rbind", k.df.list)
    
    # Here we assume method order is preserved and add column specifying the method
    methodVec <- c()
    for(j in slotNums){
      methodName <- method.dic %>% 
        filter(slotNum == j) %>% 
        pull(method)
      methodVec <- c(methodVec, rep(methodName, 
                                    nrow(k.varimp.df %>% 
                                           distinct(variables))))
    } 
                                                        
    k.varimp_fin.df <-  k.varimp.df %>% 
      mutate(method = methodVec,
             modelBlockID = k)
    
    
    sdm_varimp.df  <- rbind(sdm_varimp.df , k.varimp_fin.df)
  }
  
  return(sdm_varimp.df)
}



PCtoConifOCC <- function(con.r, notcon.r,
                         abs_percentcover_threshold = .1, 
                         pres_percentcover_threshold = .01,
                         majorityRules = F){
  # Need to produce raster with 1s for presence, 0s for absence, and NA for everything else
  # There are no zeros in the rasters, only values close to zeros. All else is NA
  
  # How much of a cell needs to be sampled in order for it to be an absence?
  # with 4km the area is 16 km^2, 10% is 400m^2 abs right?
  # abs_area_percent_cutoff <- .1 # 10%
  
  # If a cell has enough to be considered a presence or an absence, it is a 
  # presence
  
  # makes background layer of 0s and NAs
  background.r <- con.r + notcon.r
  # If absence (background - con.r) percent cover is less than the absence
  # threshold, then make it NA
  # BUT ONLY IF conifer presence is also less than the presence threshold because
  # some 0s will turn to 1s in the next step and if they are NAs then they won't 
  # be updated with the presence layer
  background.r[(background.r - con.r) < abs_percentcover_threshold &
                 con.r < pres_percentcover_threshold] <- NA
  # If absence percent cover is greater than the absence threshold than it is
  # an absence, also make conifer PC above threshold 0 as well because it will
  # be updated in next step
  background.r[(background.r - con.r) >= abs_percentcover_threshold |
                 con.r >= pres_percentcover_threshold] <- 0
  
  # Now add conifer PC and make it binary
  conifer.grid <- background.r + con.r
  
  # Decide if we want to do majority rules (absence vs presence) or if we 
  # want to be sensitive to presences (as in the HSM)
  if(majorityRules == T){
    # Only if conifer percent cover is greater than not-conifer do we consider
    # it a presence
    conifer.grid[conifer.grid >= notcon.grid] <- 1
    conifer.grid[conifer.grid < notcon.grid] <- 0
    
  } else {
    
    conifer.grid[conifer.grid >= pres_percentcover_threshold] <- 1
    conifer.grid[conifer.grid < pres_percentcover_threshold] <- 0
  }
  
  
  return(conifer.grid)
}





multiMethodEvaluation <- function(sdm_object){
  sdmstats <- c('AUC', 'COR', 'Deviance', 'obs.prevalence', 'threshold', 
                'sensitivity', 'specificity', 'TSS', 'Kappa', 'NMI', 
                'phi', 'ppv', 'npv', 'ccr', 'prevalence')
  
  method_vec <- sdm_object[[1]]@run.info$method %>% as.character()
  
  shotgun_eval.df <- evalSDMlist(sdm_object, sdmstats) %>% 
    left_join(tibble(modelID = 1:length(method_vec),
                     sdmMethod = method_vec))
  
  shotgun_eval.means <- shotgun_eval.df %>% 
    group_by(sdmMethod) %>% 
    summarise_at(vars(-modelID, -modelBlockID), mean) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), round, 3))
  
  shotgun_eval.sd_error <- shotgun_eval.df %>% 
    group_by(sdmMethod) %>% 
    summarise_at(vars(-modelID, -modelBlockID), sd) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), round, 3))
  
  sdm_eval.final <- paste(as.matrix(shotgun_eval.means), 
                          "±", 
                          as.matrix(shotgun_eval.sd_error)) %>% 
    matrix(nrow=nrow(shotgun_eval.means), 
           dimnames = dimnames(shotgun_eval.means)) %>% 
    data.frame() %>% 
    tibble() %>% 
    mutate(sdmMethod = shotgun_eval.means$sdmMethod)
  
  return(sdm_eval.final)
}
