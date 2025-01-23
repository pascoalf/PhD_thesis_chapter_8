## Map of sampling stations

library(dplyr)
library(ggplot2)
library(ggOceanMaps)
library(ggspatial)
#
reds_vector <- RColorBrewer::brewer.pal(9, "Reds")

# labels data frame
map_labels <- ASVs_diversity %>% 
  ungroup() %>% 
  select(Station, PelagicLayer, year, Latitude, Longitude) %>% 
  distinct()

#
basemap(limits = c(3, 13, 78, 79.5), 
        shapefiles = "Arctic",
        glaciers = TRUE, 
        bathymetry = TRUE,
        #        bathy.style = "rbg",
        land.col = "#eeeac4", 
        gla.col = "cadetblue", 
        land.border.col = NA, 
        gla.border.col = NA,
        grid.size = 0.05) +
  ggspatial::geom_spatial_point(data = map_labels, 
                                aes(x = Longitude, 
                                    y = Latitude),
                                size = 2) + 
  ggspatial::geom_spatial_label_repel(data = map_labels, 
                                aes(x = Longitude, 
                                    y = Latitude,
                                    label = Station), max.overlaps = 1000,
                                size = 2)+
  facet_grid(rows = c("year", "PelagicLayer")) + 
  theme(#legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 15)) + 
  scale_color_gradient(low = reds_vector[3], high = reds_vector[9])+
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tl", 
                                    which_north = "true",
                                    height = unit(0.5, "cm"),
                                    width = unit(0.5, "cm")) + 
  labs(x = "Longitude",
       y = "Latitude")
