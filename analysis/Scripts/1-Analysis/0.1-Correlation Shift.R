# Calculates correlation matrices and correlation shift for the HSM variables

library(tidyverse)
library(rgdal)
library(sp)
library(sf)
library(raster)

if(!exists("project_directory")) project_directory <- 
    "" # Enter location of project directory within quotation marks
setwd(project_directory)

## Study Area
eco_sn.sp <- st_read("Data/Masks/EcoSierra/EcoSierraMask.shp")

sf::sf_use_s2(FALSE)

# This is a CRS where the coordinates are in longitude and latitude
wgs84.crs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# Historical Climate data source
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

varNames <- hist_prism.stack %>% names()

hist.df <- values(hist_prism.stack) %>% 
  data.frame() %>% 
  na.omit()
hist_cor.df <- cor(hist.df)
  
curr.df <- values(curr_prism.stack) %>% 
  data.frame() %>% 
  na.omit()
curr_cor.df <- cor(curr.df)

# Calculate correlation shift from historic to present
cor_shift.df <- curr_cor.df - hist_cor.df

# Plotting function for each correlation matrix
corPlot <- function(cor.df, var_names = varNames){
  cor.df[upper.tri(cor.df)] <- NA
  cor.p <- cor.df %>% 
    data.frame() %>% 
    mutate(var1 = row.names(cor.df)) %>% 
    tibble() %>%
    pivot_longer(cols = !var1, names_to = "var2", values_to = "Value") %>% 
    mutate(var1 = factor(var1, levels = var_names), 
           var2 = factor(var2, levels = var_names)) %>% 
    ggplot(aes(x=var1, y=var2, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation",
                         na.value = "white") +
    theme_minimal()+ 
    coord_fixed() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.6),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  return(cor.p)
}

hist_cor.p <- corPlot(hist_cor.df)
curr_cor.p <- corPlot(curr_cor.df)
cor_shift.p <- corPlot(cor_shift.df)

hist_cor.df %>% data.frame() %>% 
  mutate(VAR1 = row.names(hist_cor.df)) %>% 
  pivot_longer(!VAR1, names_to = "VAR2", values_to = "COR") %>% 
  filter(COR != 1) %>% 
  mutate(MIN_COR = min(abs(COR)), 
         MAX_COR = max(abs(COR))) %>% 
  filter(COR %in% c(MIN_COR, MAX_COR, MIN_COR*-1, MAX_COR*-1))

curr_cor.df %>% data.frame() %>% 
  mutate(VAR1 = row.names(curr_cor.df)) %>% 
  pivot_longer(!VAR1, names_to = "VAR2", values_to = "COR") %>% 
  filter(COR != 1) %>% 
  mutate(MIN_COR = min(abs(COR)), 
         MAX_COR = max(abs(COR))) %>% 
  filter(COR %in% c(MIN_COR, MAX_COR, MIN_COR*-1, MAX_COR*-1))


cor_shift.df %>% data.frame() %>% 
  mutate(VAR1 = row.names(cor_shift.df)) %>% 
  pivot_longer(!VAR1, names_to = "VAR2", values_to = "COR") %>% 
  filter(COR != 0) %>% 
  mutate(MIN_COR = min(abs(COR)), 
         MAX_COR = max(abs(COR))) %>% 
  filter(COR %in% c(MIN_COR, MAX_COR, MIN_COR*-1, MAX_COR*-1))



hist_cor.df[floor(hist_max.i/19), hist_max.i%%19]
hist_cor.df %>% arrange()


hist_cor.p
curr_cor.p
cor_shift.p
ggsave("Results/Collinearity_historic.png", hist_cor.p, width = 7, height = 7)
write_csv(hist_cor.df %>% data.frame, file = "Results/Collinearity_historic.csv")
ggsave("Results/Collinearity_current.png", curr_cor.p, width = 7, height = 7)
write_csv(curr_cor.df %>% data.frame, file = "Results/Collinearity_current.csv")
ggsave("Results/Collinearity_shift.png", cor_shift.p, width = 7, height = 7)
write_csv(cor_shift.df %>% data.frame, file = "Results/Collinearity_shift.csv")

cor_shift.p_adj <- cor_shift.p +
  theme(legend.position = c(.5, .7))

hist_cor.p_adj <- hist_cor.p + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6)) +
  guides(fill = "none")

curr_cor.p_adj <- curr_cor.p + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6)) +
  guides(fill = "none")

cor_composite.p <- 
  
  (cor_shift.p_adj) +
  (hist_cor.p_adj / curr_cor.p_adj) +
  plot_annotation(tag_levels = "a") 


ggsave("Results/Collinearity_Composite.png", cor_composite.p, width = 10, height = 8)

