# Here we build an ensemble SDM based on weislander trees and old climate data.
# using the sdm package

# This script sources sdm_functions.R because the majority of the sdm work is done 
# in that script, and the configuration is primarily done here


library(sdm)
library(rgdal)
library(dismo)
library(usdm)
library(rasterVis)
library(ggpubr)
library(PresenceAbsence)
library(smoothr)
library(patchwork)
library(rJava)
library(sf)
library(units)
library(stars)
library(blockCV)
library(parallel)
library(raster)
library(tidyverse)
library(gt)
library(cowplot)


if(!exists("project_directory")) project_directory <- 
    "" # Enter location of project directory within quotation marks
setwd(project_directory)

# Turn spherical geometry off
sf::sf_use_s2(FALSE)


## Let's load in some data that will likely be used in all SDM runs:

# This is a CRS where the coordinates are in longitude and latitude
wgs84.crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

## New Study Area
eco_sn.sp <- st_read("Data/Masks/EcoSierra/EcoSierraMask.shp")
  
# This is a dictionary for WHR codes that should work for both CALVEG and Weislander
whr_dic <- read.csv("Data/Dictionaries/whr_dictionary.csv", stringsAsFactors = F)
source("Scripts/0-Data Prep/WHR_function.R")

# historical PRISM Data
hist_prism.stack <- raster::brick("")
names(hist_prism.stack) <- paste0("bio", 1:19)
hist_prism.stack <- hist_prism.stack %>% 
  mask(eco_sn.sp) %>% 
  crop(extent(eco_sn.sp))

# current PRISM Data
curr_prism.stack <- raster::stack("")
names(curr_prism.stack) <- paste0("bio", 1:19)
curr_prism.stack <- curr_prism.stack %>% 
  mask(eco_sn.sp) %>% 
  crop(extent(eco_sn.sp))


## Let's load in the sdm workhorse functions
source("Scripts/1-Analysis/sdm_functions.R")

# Get all variables with VIF < 10 and run vif.step 50 times and make a decision
# about which ones to keep
coll.results <- checkCollinearity(hist_prism.stack)
coll.results
pred_names <- names(coll.results$vif_namecount.tab)


## Source gymnosperm occurrence data

# Percent cover grid of conifer / non conifer at resolution of predictors
con.grid <- raster("Data/Wieslander/conifer_mhc2_percentCover.tif")
notcon.grid <- raster("Data/Wieslander/not_conifer_mhc2_percentCover.tif")

# This converts percent cover estimates to presence/absence based on thresholds
sdm_occ_Input.grid <- PCtoConifOCC(con.r = con.grid, notcon.r = notcon.grid, 
                                   abs_percentcover_threshold = .05,
                                   pres_percentcover_threshold = .05)
# plot(sdm_occ_Input.grid)

print("# Points available")
print(values(sdm_occ_Input.grid) %>% table())
# getmethodNames()

# Wrapper function so that parameters are constant across different sdm method combos
buildSDM_w_CVBlock <- function(sdmMethods, runName){
  gymno_results <- buildSDM(tree.grid = sdm_occ_Input.grid,
                            repMethod = "cv_spblock", 
                            train_N = "all",#
                            n_replicates = 5,
                            run.name = runName,
                            predictor_names = pred_names,
                            predictor.stack = hist_prism.stack[[pred_names]],
                            curr_predict.stack = curr_prism.stack[[pred_names]],
                            sdm_method = sdmMethods)
  return(gymno_results)
}




### TRAIN GAM
gymno_results <- buildSDM_w_CVBlock(sdmMethods = c("gam"),
                                    runName = "GAM")



### Evaluate GAM

# SDM results directory
SDM_dir <- "Results/..."
load(paste0(SDM_dir, "sdm_results.RData"))
gam_gymno.sdm<- sdm_results$weis.sdm
sdmstats <- c('AUC', 'COR', 'Deviance', 'obs.prevalence', 'threshold', 
              'sensitivity', 'specificity', 'TSS', 'Kappa', 'NMI', 
              'phi', 'ppv', 'npv', 'ccr', 'prevalence')


sdm_eval.means <- colMeans(sdm_eval.df %>% select(-modelBlockID, -modelID))
sdm_eval.sd_error <- sapply(sdm_eval.df %>% dplyr::select(-modelBlockID, -modelID), sd)
sdm_eval.final <- paste(round(sdm_eval.means, 3), "+/-", round(sdm_eval.sd_error, 3))
names(sdm_eval.final) <- names(sdm_eval.means)
sdm_eval.final
write_csv(sdm_eval.final %>% data.frame() %>% t() %>% data.frame(), 
          file = paste0(SDM_dir, "evaluationMetrics.csv"))




## TRAIN MULTIPLE METHODS (For Supplementary Materials)

shotgunMethods <- c("gam","rf", "brt","glm","mars","bio")
gymno_results <- buildSDM_w_CVBlock(sdmMethods = shotgunMethods,
                                    runName = "SHOTGUN")








