# This script creates a map of the Antrim Shale outcrop area in Michigan, USA, showing sampling sites and geological features.

#### Load libraries ####
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(mapview)
library(mapedit)
library(ggimage)
#setwd("/Users/quinnmoon/Downloads/Antrim_Microbiome")

#### Load geographic data ####


# U.S. states and Canadian provinces
us_states <- ne_states(country = "United States of America", returnclass = "sf")
canada_provinces <- ne_states(country = "Canada", returnclass = "sf")

# World lakes (includes Great Lakes and others)
lakes <- ne_download(scale = 10, type = "lakes", category = "physical", returnclass = "sf")

# Filter Great Lakes plus Lake St. Clair
great_lakes <- lakes %>%
  filter(name %in% c("Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie", "Lake Ontario", "Lake Saint Clair"))

basin_url <- "https://gisp.mcgi.state.mi.us/arcgis/rest/services/DEQ/Data_Miner/MapServer/26/query?where=1=1&outFields=*&outSR=4326&f=geojson"
sedimentary_basins <- st_read(basin_url)
antrim_basin <- sedimentary_basins %>%
  filter(grepl("Antrim", LABEL, ignore.case = TRUE))


#### Manually define label positions for zoomed extent ####

state_labels <- data.frame(
  name = c("Michigan", "Ohio", "Indiana", "Illinois", "Wisconsin"),
  lon = c(-84.8, -84.0, -86, -88.1, -88.0),
  lat = c(43.5, 41.6, 41.6, 41.6, 44.4)
)

province_labels <- data.frame(
  name = "Ontario",
  lon = -83,
  lat = 46.5
)

lake_labels <- data.frame(
  name = c("Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie", "Lake Ontario"),  # no St. Clair label
  lon = c(-86.9, -87.0, -82.5, -82.4, -77.8),
  lat = c(46.65, 43.0, 44.8, 41.65, 44.2)
)

site_labels <- data.frame(
  names = c("Sitting Bull D4-4 (Bull)", "Conant D1-10 (Conant)", "Gregory C3-12 (Greg)", "Little Big Horn B4-23 (Horn)", "Mortensen A4-24 (Mortensen)", "Mancelona West D2-31 (West)"),
  lon = c(-85.1551, -85.2819, -85.2181, -85.1144, -85.2139, -85.083),
  lat = c(44.9344, 45.0938, 45.0105, 44.8983, 44.9879, 44.8598)
)

zoom_box <- st_as_sf(
  st_sfc(
    st_polygon(list(rbind(
      c(-85.42, 44.75),
      c(-84.70, 44.75),
      c(-84.70, 45.33),
      c(-85.42, 45.33),
      c(-85.42, 44.75)  # Close the polygon
    ))),
    crs = 4326
  )
)



# This will open an interactive map to draw polygons. I used Wuchter et al (2013) as reference
# to draw gas producing regions
my_shape <- editMap()
plot(my_shape)
st_write(my_shape, "my_drawn_shape.geojson", append = FALSE)


# Prepare legend_shapes with consistent feature column
antrim_basin2 <- antrim_basin %>%
  mutate(feature = "Antrim Basin") %>%
  select(feature, geometry)

my_shape2 <- my_shape %>%
  mutate(feature = "Custom Shape") %>%
  select(feature, geometry)

legend_shapes <- rbind(antrim_basin2, my_shape2)

# Add feature column to site_labels
site_labels$feature <- "Site"

#### Main Map Plot ####
main_plot <- ggplot() +
  geom_sf(data = canada_provinces, fill = "gray85", color = "black", size = 0.7) +
  geom_sf(data = us_states, fill = "gray90", color = "black", size = 0.4) +
  geom_sf(data = great_lakes, fill = "lightblue", color = "black", alpha = 0.7) +
  
  # Polygon legend items (Antrim & Biogenic)
  geom_sf(data = legend_shapes, aes(fill = feature), color = NA, alpha = 0.6) +
  
  # Zoom box
  geom_sf(data = zoom_box, fill = NA, color = "red", linewidth = 0.6, linetype = "solid") +
  
  # Labels
  geom_text(data = state_labels, aes(x = lon, y = lat, label = name),
            size = 4, fontface = "bold", color = "black") +
  geom_text(data = province_labels, aes(x = lon, y = lat, label = name),
            size = 4, fontface = "bold", color = "black") +
  geom_text(data = lake_labels, aes(x = lon, y = lat, label = name),
            size = 3, fontface = "italic", color = "blue4") +
  
  # Sampling well points (with fill-based legend)
  geom_point(data = site_labels, aes(x = lon, y = lat, fill = feature),
             shape = 17, color = "red2", size = 1, show.legend = TRUE) +
  
  coord_sf(xlim = c(-88.5, -82.0), ylim = c(41.5, 46.7), expand = FALSE) +
  
  # Unified fill legend
  scale_fill_manual(
    name = NULL,
    values = c(
      "Antrim Basin" = "black",
      "Custom Shape" = "goldenrod",
      "Site" = "red2"
    ),
    labels = c(
      "Antrim Basin" = "Antrim Shale Outcrop",
      "Custom Shape" = "Biogenic Gas",
      "Site" = "Sampling Well"
    )
  ) +
  
  # Unified guide for fill aesthetic (polygons and triangle)
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = c(NA, NA, 17),              # Only Sampling Well gets shape
        color = c(NA, NA, "red2"),
        fill = c("black", "goldenrod", NA),
        size = c(NA, NA, 3)                 # Make triangle larger in legend
      ),
      order = 1
    )
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = c(1.05, 0),                # Right, bottom aligned
    legend.justification = c("left", "bottom"),
    legend.box.just = "bottom",
    legend.direction = "vertical",
    plot.margin = margin(5, 60, 5, 5)
  )

#main_plot

# Save the main plot
#ggsave("Antrim_Sampling_Map.pdf", plot = main_plot, width = 10, height = 8, dpi = 300)

#### Zoom-in plot (inset) ####
zoom_plot <- ggplot() +
  geom_sf(data = canada_provinces, fill = "gray85", color = "black", size = 0.7) +
  geom_sf(data = us_states, fill = "gray90", color = "black", size = 0.4, alpha = 0.7) +
  geom_sf(data = great_lakes, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_sf(data = antrim_basin, fill = "black", alpha = 0.6) +
  geom_sf(data = my_shape, fill = "goldenrod", color = "goldenrod", alpha = 0.4) +
  geom_point(
    data = site_labels,
    aes(x = lon, y = lat),
    color = "red2",
    size = 3,
    shape = 17  # Star shape
  ) +
  geom_label_repel(
    data = site_labels,
    aes(x = lon, y = lat, label = names),
    size = 2,
    fontface = "bold",
    color = "red2",
    fill = "white",         # White background behind text
    nudge_x = 0.05,
    direction = "y",
    hjust = 0,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.5
  ) +
  coord_sf(xlim = c(-85.42, -84.7), ylim = c(44.75, 45.33), expand = FALSE) +
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
    panel.border = element_rect(color = "red", fill = NA, linewidth = 1),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

#zoom_plot

#ggsave("Antrim_Sampling_Map_zoom.pdf", plot = zoom_plot_with_border, width = 10, height = 8, dpi = 300)

#### Combined Map Plot ####
combined_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(zoom_plot, 
            x = 0.51,   # slight shift left to balance increased width
            y = 0.30,   # slight shift down if necessary
            width = 0.45,  # increased size
            height = 0.45) +
  draw_line(x = c(0.50, 0.6653), y = c(0.70, 0.7394), color = "red", size = 0.5) +  # top-left
  draw_line(x = c(0.5023, 0.6653), y = c(0.647, 0.315), color = "red", size = 0.5)   # bottom-left

#combined_plot

ggsave("Antrim_Sampling_Map_combined.pdf", plot = combined_plot, width = 21, height = 8)


