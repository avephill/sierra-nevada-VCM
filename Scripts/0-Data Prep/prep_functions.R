# Various functions shared by scripts in this folder


# This function stitches together the 3 CALVEG vegetation files that makeup the entirety of the study area

calvegSource.dir <- ""
getCALVEGvegSF <- function(WHRs, type = "raster"){
  print("loading north interior...")
  norinterior.sf <- st_read(paste(calvegSource.dir, "S_USA.EVMid_R05_NorthInterior.gdb/")) %>% 
    select(CWHR_TYPE) %>% 
    filter(CWHR_TYPE %in% WHRs) %>%  
    st_union() %>% 
    st_transform(st_crs(eco_sn.sp)) %>% 
    st_intersection(eco_sn.sp) %>% 
    st_sf()
  
  print("loading north sierra...")
  norsierra.sf <- st_read(paste(calvegSource.dir,"S_USA.EVMid_R05_NorthSierra.gdb/")) %>% 
    select(CWHR_TYPE) %>% 
    filter(CWHR_TYPE %in% WHRs) %>% 
    st_union() %>% 
    st_transform(st_crs(eco_sn.sp)) %>% 
    st_intersection(eco_sn.sp) %>% 
    st_sf()
  
  print("loading south sierra...")
  sousierra.sf <- st_read(paste(calvegSource.dir,"S_USA.EVMid_R05_SouthSierra.gdb/")) %>% 
    select(CWHR_TYPE) %>% 
    filter(CWHR_TYPE %in% WHRs) %>% 
    st_union() %>% 
    st_transform(st_crs(eco_sn.sp)) %>% 
    st_intersection(eco_sn.sp) %>% 
    st_sf()
  
    
  
  print("Merging...")
  if(type == "raster"){
    norinterior.r <- rasterize(norinterior.sf %>% as("Spatial"), bl.grid, 
                               getCover = T, background = 0, progress = "text")
    norsierra.r <- rasterize(norsierra.sf %>% as("Spatial"), bl.grid, 
                             getCover = T, background = 0, progress = "text")
    sousierra.r <- rasterize(sousierra.sf %>% as("Spatial"), bl.grid, 
                             getCover = T, background = 0, progress = "text")
    
    tot <- norinterior.r + norsierra.r + sousierra.r
  } else if (type == "polygon"){
    tot <- rbind(norinterior.sf, norsierra.sf, sousierra.sf) %>% 
      st_union()
  }
  
  return(tot)
}
