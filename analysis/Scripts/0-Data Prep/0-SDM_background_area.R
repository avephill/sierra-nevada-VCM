# We define the study area based on EPA Ecotype "Northwestern Forested Mountains"
# and put a buffer around it 

library(sf)

# This is a CRS where the coordinates are in longitude and latitude
wgs84.crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
# Turn spherical geometry off
sf::sf_use_s2(FALSE)

sierra_full.shp <- st_read("Data/Masks/SierraFire_Mask/") %>% 
  st_transform(wgs84.crs)
sierra_boundingbox <- extent(sierra_full.shp)
sierra_boundingbox@ymax <- 40.9
eco_sn.sp <- st_read("../Commons/Boundaries/Natural/NA_CEC_Eco_Level1/NA_CEC_Eco_Level1.shp") %>% 
  st_transform(wgs84.crs) %>% 
  st_crop(sierra_boundingbox) %>% 
  filter(NA_L1NAME == "NORTHWESTERN FORESTED MOUNTAINS") %>% 
  st_buffer(45/60)

st_write(eco_sn.sp[1], "Data/Masks/EcoSierra/EcoSierraMask40_9.shp", delete_layer = T)

# check
plot(readOGR("Data/Masks/EcoSierra/EcoSierraMask40_9.shp"))

# So this is all well and good, however, the weislander survey only has bits
# and spurts above 40 degrees latitude in the sierras. this is very close to mt
# Lassen. Furthermore, if we limit it to south of 40 degrees then we don't have
# to do any spatial extrapolation in the predictions. Going to cut the model
# to below 40 degrees



buffered.bb <- extent(eco_sn.sp)
buffered.bb@ymax <- 40.0
eco_sn40.sp <- eco_sn.sp %>% 
  st_crop(buffered.bb)

st_write(eco_sn40.sp[1], "Data/Masks/EcoSierra/EcoSierraMask.shp", delete_layer = T)

# check
plot(readOGR("Data/Masks/EcoSierra/EcoSierraMask.shp"))
