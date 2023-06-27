## Load packages
  library(akgfmaps)
  library(sf)
  library(ggmap)
  library(rgdal)
  library(tidyverse)
  library(RColorBrewer)
  library(patchwork)
  library(ggpubr)

## Load map layers
  map_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")


## Trim survey grid to survey area
  map_layers$survey.grid %>%
    st_transform(crs = st_crs(map_layers$survey.area)) %>%
    st_intersection(map_layers$survey.area) -> survey.grid
  
## Specify plot boundary, transform to map crs
  plot.boundary <- akgfmaps::transform_data_frame_crs(data.frame(x = c(-180, -150), 
                                                                 y = c(54.5, 67)), 
                                                      out.crs = map_layers$crs) 

## Plot 
  ggplot() +
    geom_sf(data = map_layers$bathymetry, color=alpha("grey70")) +
    geom_sf(data = map_layers$survey.grid, fill=NA, color=alpha("grey70"), linewidth = 1)+
    geom_sf(data = map_layers$survey.area, fill = NA, linewidth = 1) +
    geom_sf(data = map_layers$akland, fill = "grey80", size=0.1) +
    geom_sf(data = map_layers$survey.area, fill = NA) +
    
    scale_x_continuous(breaks = c(-180, -175, -170, -165, -160, -155, -150), labels = paste0(c(180, 175, 170, 165, 160, 155, 150), "°W")) + 
    #scale_y_continuous(breaks = c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), labels = paste0(c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), "°N")) +
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    theme_bw()+
    theme(panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_rect(fill = NA, color = "black"),
          legend.key = element_rect(fill = NA, color = "grey70"),
          legend.key.size = unit(0.65,'cm'),
          legend.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10), 
          legend.title = element_text(size = 10),
          plot.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank()) -> study_map
  
## Save map
  ggsave(plot = study_map, "./Figures/study_map.png", height = 7, width = 8, units = "in")
  
