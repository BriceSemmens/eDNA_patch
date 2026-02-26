# CalCOFI Station Map

# Load the necessary libraries
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(ggspatial)
library(ggOceanMaps)
library(grid)
library(cowplot)
library(ggnewscale)
library(stringr)

# Load the datasets
station_map <- read.csv("/Users/nastassiapatin/OneDrive - UC San Diego/CalCOFI eDNA Data Management/CalCOFIStationOrder.csv")
data <- read.csv("/Users/nastassiapatin/OneDrive - UC San Diego/Manuscripts/Delphinid detection probability/intercal_ALL_metadata_011526.csv")

save_png = "/Users/nastassiapatin/OneDrive - UC San Diego/Manuscripts/Delphinid detection probability/intercal_station_samples.png"
save_svg = "/Users/nastassiapatin/OneDrive - UC San Diego/Manuscripts/Delphinid detection probability/intercal_station_samples.svg"

# Define the map limits for the California coast
limits <- c(-124, -116.5, 32.5, 37)  

# Create the basemap with ggOceanMaps
basemap_ocean <- basemap(limits = limits, bathy.style = "rcb", 
                         bathymetry = TRUE, grid.col = NA)

### Map of stations and filters taken

# Get stations for each bio rep
station_bioreps <- data %>% 
  select(unique_biorep_numeric, station) %>%
  distinct(unique_biorep_numeric, station) %>%
  group_by(station) %>%
  summarize(num_bioreps = n())

# Merge filtered eDNA data with station coordinates

# Add "station" column to station sheet
station_map <- station_map %>%
  dplyr::rename(station = Station_key) 

# Fix station decimal to match station map
# Remove underscore to match eDNA data
station_map$station <- station_map$station %>%
  str_replace_all("70", "70.0") %>%
  str_replace_all("_", " ")

# Merge station map with sample data
station_map_ecalcofi <- station_map %>%
  inner_join(station_bioreps, by = "station") %>%
  select(station, Lat, Lon, num_bioreps) %>%
  distinct(station, Lat, Lon, num_bioreps)

# Make the map!

map_ecalcofi <- basemap_ocean +
  scale_fill_gradientn(colours = c("#08306B", "#08519C", "#2171B5",
                                   "#4292C6", "#6BAED6", "#9ECAE1",
                                   "#C6DBEF"),
                       name = "Bottom Depth (m)",
                       trans = "reverse") +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = station_map_ecalcofi,
    aes(x = Lon, y = Lat, fill = num_bioreps),
    shape = 21, color = "black",
    stroke = 1, size = 3
  ) +
  scale_fill_gradient(
    low = "lightgray", high = "darkred",
    limits = c(6, 19),
    breaks = seq(6, 19, by = 4),
    name = "No. filters (all depths)") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank()
  ) #+
  #ggtitle("Number of filters collected at each station")

ggsave(filename = save_png, plot = map_ecalcofi,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")

ggsave(filename = save_svg, plot = map_ecalcofi,
       width = 12, height = 8,dpi = 300, units = "in")

### Map delphinid detections by depth

# Get stations and depths for each bio rep
station_depth_biorep <- data %>% 
  select(unique_biorep_numeric, station, Depth_m) %>%
  distinct(unique_biorep_numeric, station, Depth_m)

# Sum detections for the technical replicates of every filter
data_filt <- data %>%
  select(delphinus_all, unique_biorep_numeric) %>%
  group_by(unique_biorep_numeric) %>%
  summarize(delphinid_techreps = sum(delphinus_all)) 

# Classify detection presence/absence of every filter
data_filt <- data_filt %>%
  mutate(delph_biorep_any = if_else(delphinid_techreps > 0, 1, 0))

# Aggregate by station and depth
data_map <- data_filt %>%
  left_join(station_depth_biorep, by = "unique_biorep_numeric") %>%
  group_by(station, Depth_m) %>%
  summarize(delph_sum = sum(delph_biorep_any))

# Add in coordinates for stations
data_map <- data_map %>%
  left_join(station_map_ecalcofi, by = "station")

# Group depths
data_map_depthgroups <- data_map %>%
  mutate(depth_group = case_when(
    Depth_m < 11 ~ "1-10",
    Depth_m >= 11 & Depth_m <= 101 ~ "15-100",
    TRUE ~ "101-515m")
  )

# Make maps for each depth group

# 1-10 m
data_map_1 <- data_map_depthgroups %>%
  filter(depth_group == "1-10")

map_edna_d1 <- basemap_ocean +
  scale_fill_gradientn(colours = c("#08306B", "#08519C", "#2171B5",
                                   "#4292C6", "#6BAED6", "#9ECAE1",
                                   "#C6DBEF"),
                       name = "Bottom Depth (m)",
                       trans = "reverse") +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = data_map_1,
    aes(x = Lon, y = Lat, fill = delph_sum),
    shape = 21, color = "black", stroke = 1, size = 3
  ) +
  scale_fill_gradientn(
    colors = c("white", "salmon", "maroon", "darkred"),
    limits = c(0, 5),
    breaks = c(0, 1, 2, 3, 4, 5),
    name = "No. Biological Replicates") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank()
  ) +
  ggtitle("Delphinid detections, sample depth group: 1-10 m")

ggsave(filename = "delph_detection_map_1.png", plot = map_edna_d1,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")

ggsave(filename = "delph_detection_map_1.svg", plot = map_edna_d1,
       width =8, height = 5, dpi = 300, units = "in")

# 11-100 m

data_map_2 <- data_map_depthgroups %>%
  filter(depth_group == "15-100")

map_edna_d2 <- basemap_ocean +
  scale_fill_gradientn(colours = c("#08306B", "#08519C", "#2171B5",
                                   "#4292C6", "#6BAED6", "#9ECAE1",
                                   "#C6DBEF"),
                       name = "Bottom Depth (m)",
                       trans = "reverse") +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = data_map_2,
    aes(x = Lon, y = Lat, fill = delph_sum),
    shape = 21, color = "black", stroke = 1, size = 3
  ) +
  scale_fill_gradientn(
    colors = c("white", "salmon", "maroon", "darkred"),
    limits = c(0, 5),
    breaks = c(0, 1, 2, 3, 4, 5),
    name = "No. Biological Replicates") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank()
  ) +
  ggtitle("Delphinid detections, sample depth group: 15-100 m")

ggsave(filename = "delph_detection_map_2.png", plot = map_edna_d2,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")

ggsave(filename = "delph_detection_map_2.svg", plot = map_edna_d2,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")

# 101-515 m

data_map_3 <- data_map_depthgroups %>%
  filter(depth_group == "101-515m")

map_edna_d3 <- basemap_ocean +
  scale_fill_gradientn(colours = c("#08306B", "#08519C", "#2171B5",
                                   "#4292C6", "#6BAED6", "#9ECAE1",
                                   "#C6DBEF"),
                       name = "Bottom Depth (m)",
                       trans = "reverse") +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = data_map_3,
    aes(x = Lon, y = Lat, fill = delph_sum),
    shape = 21, color = "black", stroke = 1, size = 3
  ) +
  scale_fill_gradientn(
    colors = c("white", "salmon", "maroon", "darkred"),
    limits = c(0, 5),
    breaks = c(0, 1, 2, 3, 4, 5),
    name = "No. Biological Replicates") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_blank()
  ) +
  ggtitle("Delphinid detections, sample depth group: 101-515 m")

ggsave(filename = "delph_detection_map_3.png", plot = map_edna_d3,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")

ggsave(filename = "delph_detection_map_3.svg", plot = map_edna_d3,  
       width = 8, height = 5, dpi = 300, units = "in",
       bg = "white")
