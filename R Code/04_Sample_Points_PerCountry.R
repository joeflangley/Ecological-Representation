# Load libraries
library(units)
library(giscoR)
library(tidyverse)
library(sf)
library(units)

# Read in the country names vector and then the country boundaries
countries <- st_read("data/raw_data/country_boundaries/countries.shp")

# Calculate the number of sampled points proportional to each country's area
samp_points <- vector("numeric")

# Loop through each country
for (i in 1:nrow(countries)) {
  # Calculate the area of the country in km^2 and drop the units
  country_area <- st_area(countries[i, ]) %>% set_units(km2) %>% drop_units()
  
  # Calculate the proportional value for the number of sample points, by dividing the country's
  # area by 25 km^2 and rounding to the nearest integer
  prop <- round(country_area / 5^2, 0)
  
  # Ensure that the number of sample points is at least 1000
  if (prop < 1000) {
    prop <- 1000
    # Ensure that the number of sample points does not exceed 50000
  } else if (prop > 50000) {
    prop <- 50000
  }
  samp_points[i] <- prop
}

# Add the number of sampled points to a new column in countries
countries$samp_points <- samp_points

# Create a function that takes a stratified sample of points
sample_func <- function(x, n) {
  sp::spsample(as_Spatial(x), n, type = "stratified")
}

set.seed(9199) # set seed

spatial_points <- list()

for (i in 1:nrow(countries)) {
  # Take a sample of points from each country, equivalent to the proportional value we calculated
  dat <- sample_func(x = countries[i, ], n = countries$samp_points[i]) %>% 
    as.data.frame(.) %>% 
    
    # Add the ISO3 code, country name and number of sampled points to the dataframe
    cbind(., 
          ISO3_CODE = countries$ISO3_CODE[i],
          country_name = countries$NAME_ENGL[i],
          IUCN_name = countries$IUCN_name[i],
          samp_points = countries$samp_points[i]) %>% 
    
    # Give coordinate columns x and y names
    rename(x = x1, y = x2)
  
  spatial_points[[i]] <- dat
}

# Combine into one dataframe
spatial_points <- do.call(rbind, spatial_points)

# Save as a csv
write.csv(spatial_points, "data/raw_data/spatial_points.csv", row.names = F)
