# Function that's called throughout the project that returns the list of WHR types
# that are either gymnosperm or angiosperm dominated
getWHR <- function(vegGroup){
  # This is a dictionary for WHR codes that should work for both CALVEG and Weislander
  whr_dic <- read.csv("Data/Dictionaries/whr_dictionary.csv",stringsAsFactors = F)
  
  # Exclude Palm Oasis, Joshua Tree, Desert Wash, Montane and Foothill Riparia, Sagebrush, Juniper
  # Urban, Riparian, and types that are also listed as trees e.g. (ASP, JUN)
  unwantedWHR <- c("ADS","POS", "JST", "DSW", "VRI", "MRI", "SGB", "JUN", "ASP", "DRI","URB")
  
  if(vegGroup == "gymno"){
    # Separate Lifeform codes into those that are trees vs shrubs
    gymnoWHR_LF = c("WHR_CON") 
    
    # Separate WHR into trees and shrubs using corresponding Lifeform codes
    gymnoWHR = unique(whr_dic[whr_dic$CWHR.Lifeform.Code %in% gymnoWHR_LF,"CWHR.Type.Code"])
    # Also exclude Blue-Oak Foothill Pine
    gymnoWHR <- gymnoWHR[!(gymnoWHR %in% c(unwantedWHR, "BOP"))] %>% c("MHC")
    
    veg.df <- whr_dic %>% filter(CWHR.Type.Code %in% gymnoWHR) %>% 
      dplyr::select(CWHR.Type.Description, CWHR.Type.Code, Raster_VALUE) %>% 
      na.omit %>% distinct()
    
  } else if(vegGroup == "angio"){
    
    angioWHR_LF = c("WHR_SHB", "WHR_HDW", "WHR_MIX")
    
    angioWHR = unique(whr_dic[whr_dic$CWHR.Lifeform.Code %in% angioWHR_LF, "CWHR.Type.Code"])
    angioWHR <- angioWHR[!(angioWHR %in% unwantedWHR)]
    
    veg.df <- whr_dic %>% filter(CWHR.Type.Code %in% angioWHR) %>% 
      dplyr::select(CWHR.Type.Description, CWHR.Type.Code, Raster_VALUE) %>% 
      na.omit %>% distinct()
    
  }
  
  return(veg.df)
}
