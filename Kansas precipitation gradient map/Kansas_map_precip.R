# Kansas Precipitation Map
# Nichole Ginnan (nginn001@ucr.edu)
# July 18, 2025

library(tigris)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)

options(tigris_use_cache = TRUE)

# Load precipitation data
precip_data <- read_csv("/Users/nicholeginnan/Documents/KU/Luteibacter/scripts/PrecipTotal-2010-2020.csv")
# calculate averages
county_avg_precip <- precip_data %>%
  group_by(County) %>%
  summarize(mean_precip = mean(precip, na.rm = TRUE)) %>%
  ungroup()

#Load Kansas county shapefile
ks_counties <- counties(state = "KS", cb = TRUE, class = "sf") %>%
  mutate(NAME = tolower(NAME))  # Normalize names for joining

# Make sure County names are also lowercase to match
county_avg_precip <- county_avg_precip %>%
  mutate(County = tolower(County))

#Join the precipitation data to the map
ks_map <- left_join(ks_counties, county_avg_precip, by = c("NAME" = "County"))
# Define site coordinates
site_points <- data.frame(
  Site = c("SVR", "HAY", "KNZ", "TLI"),
  Lat = c(38.87416, 38.83564, 39.10183, 38.76780),
  Long = c(-100.98318, -99.30332, -96.59624, -97.5664))
# Convert to sf, using WGS84 (EPSG:4326), then transform to match your map
site_sf <- st_as_sf(site_points, coords = c("Long", "Lat"), crs = 4326)
# Match CRS with map shapefile
site_sf <- st_transform(site_sf, crs = st_crs(ks_map))
# Plot the map
ggplot(ks_map) +
  geom_sf(aes(fill = mean_precip), color = "grey", size = 0.2) +
  scale_fill_gradientn(
    colors = c("#7d3049","#aa4465", "#ffa69e", "#3d5a80","#2b415d","#1d2e43"),
    name = "Mean Annual Precip 2010-2020 (mm)",
    na.value = "gray90") +
  geom_sf(data = site_sf, color = "white", size = 3) +
  geom_sf_text(data = site_sf, aes(label = Site), color = "white", size = 7, fontface = "bold", nudge_y = -0.3)+
  theme_minimal() +
  labs(fill = "Precip (mm)") +
 theme(plot.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom")
