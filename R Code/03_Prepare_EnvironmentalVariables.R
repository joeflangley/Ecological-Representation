# Load libaries
library(tidyverse)
library(geodata)
library(terra)
library(giscoR)
library(DescTools)
library(car)
library(sf)

######################### ---- Environmental Data Download ---- ########################## 
# Download the worldclim data
worldclim_global(var = "bio", res = 0.5, path = "data/raw_data/") # bio1:bio19

# Download the soil data
soil_world(var = "sand", depth = 15, stat = "mean", path = "data/raw_data/") # sand content
soil_world(var = "silt", depth = 15, stat = "mean", path = "data/raw_data/") # silt content
soil_world(var = "clay", depth = 15, stat = "mean", path = "data/raw_data/") # clay content
soil_world(var = "phh2o", depth = 15, stat = "mean", path = "data/raw_data/") # water pH
soil_world(var = "nitrogen", depth = 15, stat = "mean", path = "data/raw_data/") # total nitrogen
soil_world(var = "soc", depth = 15, stat = "mean", path = "data/raw_data/") # soil organic carbon
soil_world(var = "ocd", depth = 30, stat = "mean", path = "data/raw_data/") # organic carbon density

# Download the elevation data
elevation_global(res = 0.5, path = "data/raw_data/")

# Calculate terrain variables from the elevation model
elevation <- rast("data/raw_data/wc2.1_30s/wc2.1_30s_elev.tif")
elevation[is.na(elevation)] <- 0 # replace NA values with 0

terrain(x = elevation, v = "slope", neighbors = 8, unit = "degrees", filename = "data/raw_data/wc2.1_30s/slope.tif", overwrite = T) # slope
terrain(x = elevation, v = "aspect", neighbors = 8, unit = "degrees", filename = "data/raw_data/wc2.1_30s/aspect.tif", overwrite = T) # aspect
terrain(x = elevation, v = "TPI", filename = "data/raw_data/wc2.1_30s/TPI.tif", overwrite = T) # topographic position index
terrain(x = elevation, v = "TRI", filename = "data/raw_data/wc2.1_30s/TRI.tif", overwrite = T) # terrain ruggedness index

############################## ---- Cropping Rasters ---- ################################ 
# Read in the country names vector and then the country boundaries
countries <- read.csv("data/output_data/country_names.csv")[,1]
countries <- gisco_get_countries(country = countries, year = "2020", resolution = "01") %>% # country boundaries
  dplyr::select(ISO3_CODE, NAME_ENGL) 

# Read in the soil, bioclimatic and topographic spatial rasters
soil_files <- list.files("data/raw_data/soil_world", pattern = ".tif", full.names = TRUE)
soil_rasters <- lapply(soil_files, rast)
names(soil_rasters) <- str_remove_all(basename(soil_files), ".tif")

bio_topo_files <- list.files("data/raw_data/wc2.1_30s/", pattern = ".tif", full.names = TRUE)
bio_topo_rasters <- lapply(bio_topo_files, rast)
names(bio_topo_rasters) <- str_remove_all(basename(bio_topo_files), ".tif")

# Put all the variables in one list
env_dat <- c(soil_rasters, bio_topo_rasters)

# Convert country boundaries to a spatial vector
v <- vect(countries)

# Crop each environmental raster to the extent of the country vector
cropped_rasters <- vector("list", length = length(env_dat))
for (i in 1:length(env_dat)) {
  cropped_rast <- crop(env_dat[[i]], ext(v))
  cropped_rasters[[i]] <- cropped_rast
}

# Mask the cropped rasters by the country vector to remove data outside the studied countries
masked_rasters <- vector("list", length = length(env_dat))
for (i in 1:length(env_dat)) {
  masked_rast <- mask(cropped_rasters[[i]], v)
  masked_rasters[[i]] <- masked_rast
}
names(masked_rasters) <- names(env_dat)

# Project and save these masked rasters
for (i in 1:length(masked_rasters)) {
  r <- project(masked_rasters[[i]], "+proj=moll")
  file_name <- paste0("data/output_data/masked_env_data/", names(masked_rasters)[i], ".tif")
  writeRaster(r, filename = file_name, overwrite = TRUE)
}

cropped_rasters <- list()
for (i in 1:4) {
  cropped_rast <- crop(bio_topo_rasters[[i]], ext(v))
  cropped_rasters[[i]] <- cropped_rast
}

masked_rasters <- list()
for (i in 1:4) {
  masked_rast <- mask(cropped_rasters[[i]], v)
  masked_rasters[[i]] <- masked_rast
}
names(masked_rasters) <- names(bio_topo_rasters[1:4])

# Project and save these masked rasters
for (i in 1:length(masked_rasters)) {
  r <- project(masked_rasters[[i]], "+proj=moll")
  file_name <- paste0("data/output_data/masked_env_data/", names(masked_rasters)[i], ".tif")
  writeRaster(r, filename = file_name, overwrite = TRUE)
}

############################ ---- Environmental Sample ---- ############################## 
countries <- countries %>% 
  st_transform(., crs = "+proj=moll")

spatial_sample <- sp::spsample(as_Spatial(countries), 50000, type = "stratified") # stratified sample of points
sample_vect <- vect(spatial_sample)

masked_rasters <- list.files("data/output_data/masked_env_data/", full.names = T)
masked_rasters <- lapply(masked_rasters, rast) # read in the masked environmental rasters

env_dat <- rast(masked_rasters) # combine into one spatial raster
env_dat <- extract(env_dat, sample_vect) # extract the environmental data from the sampled points

env_dat <- env_dat %>% 
  drop_na() %>% # remove NAs
  dplyr::select(-ID) # remove ID column

######################## ---- Tests for Multicollinearity ---- ########################### 
# Create a correlation matrix
CorMat <- cor(env_dat, method='spearman') 

# Remove highly correlated variables 
colRemove <- FindCorr(CorMat, cutoff = 0.7, verbose = FALSE)
env_slim <- env_dat[,-colRemove]

# Test using Variance Inflation Factor
env_slim <- env_slim %>% mutate(dummy = c(1)) # add a dummy response variable for the linear model
model <- lm(dummy ~ ., data = env_slim)
vif(model) # model shows <5 VIF for all variables, which is good

# # Alternative method for removing collinearity
# library(virtualspecies)
# a <- removeCollinearity(env_dat, select.variables = TRUE, sample.points = TRUE, nb.points = 20000)

######################## ---- Save Retained Variable Names ---- ########################## 
# Get the column names from the reduced environmental dataset
retained_variables <- colnames(env_slim)[-14] # remove dummy variable

# Save them as a csv file
write.csv(retained_variables, "data/output_data/retained_variable_names.csv", row.names = FALSE)

# Read the retained variables into a vector
retained_variables <- read.csv("data/output_data/retained_variable_names.csv")[,1]
pattern <- paste(retained_variables, collapse = "|")
env_dat <- list.files(path = "data/output_data/masked_env_data/", pattern = pattern, full.names = TRUE)

env_dat <- lapply(env_dat, rast)
env_dat <- rast(env_dat)
