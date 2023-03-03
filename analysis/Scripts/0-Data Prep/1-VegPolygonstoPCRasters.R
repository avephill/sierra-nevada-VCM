# This Script takes CALVEG and Wieslander vegetation vector data and converts them 
# to rasters where the cell values are percent cover of angiosperm dominated or gymnosperm
# dominated vegetation

library(sf)
library(raster)
library(stars)
library(tidyverse)
library(rgdal)


if(!exists("project_directory")) project_directory <- 
    "" # Enter location of project directory within quotation marks
setwd(project_directory)

source("Scripts/0-Data Prep/WHR_function.R")
source("Scripts/0-Data Prep/prep_functions.R")

# This is a CRS where the coordinates are in longitude and latitude
wgs84.crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

## New Study Area
eco_sn.sp <- st_read("Data/Masks/EcoSierra/EcoSierraMask.shp")

source("Scripts/0-Data Prep/WHR_function.R")

# Turn spherical geometry off
sf::sf_use_s2(FALSE)

# Get one climate layer for which to base the output rasters on
clim.grid <- raster::brick("")[[1]]
values(clim.grid) <- 0 # Blank it out
bl.grid <- clim.grid %>% mask(eco_sn.sp) %>% crop(eco_sn.sp)

# Polygonize the blank raster into square polygons
bl.sf <- rasterToPolygons(bl.grid) %>% as("sf") %>% 
  st_transform(wgs84.crs) %>% 
  mutate(CELL_AREA = st_area(.))

# Read in Statewide Wieslander Data
wies.sf <- read_sf("")


wies_study.sf <- wies.sf %>% st_transform(st_crs(eco_sn.sp)) %>% 
  st_intersection(eco_sn.sp)





# Here we distinguish Conifers vs. Not Conifers
# We pull MHC (montane hardwood-conifer) out so we can split it between the two
con.sf <- wies_study.sf %>% filter(WHR1 %in% (getWHR("gymno") %>% 
                                                pull("CWHR.Type.Code") %>% 
                                                str_subset("MHC", negate = T)))
notcon.sf <- wies_study.sf %>% filter(WHR1 %in% (getWHR("angio") %>% 
                                                   pull("CWHR.Type.Code") %>% 
                                                   str_subset("MHC", negate = T)))
mixed.sf <- wies_study.sf %>% filter(WHR1 == "MHC")
print(paste0("Montane HW-CON: ", mixed.sf %>% pull(AREA_M2) %>% sum() / 1000000, " km2"))



# So we take the sf and then we st_union them into one feature, and
# then we calculate the percent cover of the sf within each cell of bl.grid
con_simp.shp <- con.sf %>% 
  st_union() %>% 
  as("Spatial")
con.r <- rasterize(con_simp.shp, bl.grid, getCover = T, background = 0, progress = "text")

notcon_simp.shp <- notcon.sf %>% 
  st_union() %>% 
  as("Spatial")
notcon.r <- rasterize(notcon_simp.shp, bl.grid, getCover = T, background = 0, progress = "text")

mixed_simp.shp <- mixed.sf %>% 
  st_union() %>% 
  as("Spatial")
mixed.r <- rasterize(mixed_simp.shp, bl.grid, getCover = T, background = 0, progress = "text")


con_mhc.r <- con.r + mixed.r/2
notcon_mhc.r <- notcon.r + mixed.r/2

# We're writing the mixed to disk so that we can separate it out later if we want to from con.r and notcon.r
writeRaster(mixed.r, "Data/Wieslander/mixed_percentCover.tif", overwrite = T)
writeRaster(con_mhc.r, "Data/Wieslander/conifer_mhc2_percentCover.tif", overwrite = T)
writeRaster(notcon_mhc.r, "Data/Wieslander/not_conifer_mhc2_percentCover.tif", overwrite = T)




#### NOW CALVEG

# Produced this from calveg .gdb files, wherein I loaded them into QGIS,
# Vector > Data Management Tools > Merge Vector Layers
# Then saved as shapefile



# So we take the sf and then we st_union them into one feature, and
# then we calculate the percent cover of the sf within each cell of bl.grid
gymnos <- getWHR("gymno") %>% 
  pull("CWHR.Type.Code") %>% 
  str_subset("MHC", negate = T)



angios <- getWHR("angio") %>% 
  pull("CWHR.Type.Code") %>% 
  str_subset("MHC", negate = T)

calcon.r <- getCALVEGvegSF(gymnos)
calnotcon.r <- getCALVEGvegSF(angios)
calmixed.r <- getCALVEGvegSF("MHC")

calcon_mhc.r <- calcon.r + calmixed.r/2
calnotcon_mhc.r <- calnotcon.r + calmixed.r/2

# We're writing the calmixed to disk so that we can separate it out later if we want to from calcon.r and calnotcon.r
writeRaster(calmixed.r, "Data/CALVEG/mixed_percentCover.tif", overwrite = T)
writeRaster(calcon_mhc.r, "Data/CALVEG/conifer_mhc2_percentCover.tif", overwrite = T)
writeRaster(calnotcon_mhc.r, "Data/CALVEG/not_conifer_mhc2_percentCover.tif", overwrite = T)










