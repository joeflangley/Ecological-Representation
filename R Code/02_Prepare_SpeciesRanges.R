# Load libraries
library(sf)
library(tidyverse)
library(rCAT)
sf_use_s2(FALSE)

############################## ---- Threatened Species ---- ##############################
# Read in the threatened species data and add columns characterising the taxon and threat status
t_birds <- st_read("all_data/species_data/threatened/birds/data_0.shp") %>% 
  mutate(taxon = "Birds", status = "Threatened") 
  
t_mammals <- st_read("all_data/species_data/threatened/mammals/data_0.shp") %>% 
  mutate(taxon = "Mammals", status = "Threatened")

t_herps <- st_read("all_data/species_data/threatened/herptiles/data_0.shp") %>% 
  mutate(taxon = "Herptiles", status = "Threatened")
  
t_plants <- st_read("all_data/species_data/threatened/plants/data_0.shp") %>% 
  mutate(taxon = "Plants", status = "Threatened") 

# Combine into one dataframe and filter PRESENCE and ORIGIN based on p. 204 of the KBA 
# guidelines
threatened_species <- rbind(t_birds, t_mammals, t_herps, t_plants) %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, taxon, status) %>% # retain relevant columns
  st_transform(., crs = "+proj=moll") # reproject to Mollweide projection

# Identify any invalid geometries and fix them
invalid_geom <- !st_is_valid(threatened_species) # identify invalid geometries
table(invalid_geom)
threatened_species$geometry[invalid_geom] <- st_make_valid(threatened_species$geometry[invalid_geom]) # fix only the invalid geometries

# Group by species and merge geometries while retaining the values of non-geometry columns
threatened_species <- threatened_species %>%
  group_by(SCI_NAME) %>%
  summarise(
    across(-geometry, first), # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the taxon and status columns should be the same for each record of the same species
    geometry = st_union(geometry), # unionise the geometries for each species
  ) %>% 
  relocate(ID_NO) # put ID column back in the first position

# Read in the species summary data
threatened_summary <- rbind(read.csv("all_data/species_data/threatened/birds/simple_summary.csv"),
                            read.csv("all_data/species_data/threatened/mammals/simple_summary.csv"),
                            read.csv("all_data/species_data/threatened/herptiles/simple_summary.csv"),
                            read.csv("all_data/species_data/threatened/plants/simple_summary.csv")) %>% 
  rename(SCI_NAME = scientificName, ID_NO = internalTaxonId) %>% # rename columns to align with polygon column names
  dplyr::select(ID_NO, redlistCategory) # retain relevant columns

# Merge the geometry data with the summary data to get each species' Red List category
threatened_species <- left_join(threatened_species, threatened_summary, by = "ID_NO") %>% 
  relocate(geometry, .after = redlistCategory)
write_sf(threatened_species, "data/output_data/prepared_species_ranges/threatened_species.gpkg")

############################ ---- Non Threatened Species ---- ############################
# Here we will do taxonomic groups separately due to the large number of species in each

# Function to make geometries valid, dissolve by SCI_NAME and add summary data
prepare_species <- function(species_data) {
  # Identify any invalid geometries and fix them
  invalid_geom <- !st_is_valid(species_data) # identify invalid geometries
  species_data$geometry[invalid_geom] <- st_make_valid(species_data$geometry[invalid_geom]) # fix only the invalid geometries
  
  # Group by species and merge geometries while retaining the values of non-geometry columns
  species_data <- species_data %>%
    group_by(SCI_NAME) %>%
    summarise(
      across(-geometry, first),        # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the taxon and status columns should be the same for each record of the same species
      geometry = st_union(geometry),   # union the geometries for each species
    ) %>% 
    relocate(ID_NO) # put ID column back in the first position
  
  # Merge the geometry data with the summary data to get each species' Red List category
  species_data <- left_join(species_data, non_threatened_summary, by = "ID_NO") %>% 
    relocate(geometry, .after = redlistCategory)
}

# Read in the summary data for all non-threatened species
non_threatened_summary <- rbind(read.csv("all_data/species_data/non_threatened/plants/part1/simple_summary.csv"),
                                read.csv("all_data/species_data/non_threatened/plants/part2/simple_summary.csv"),
                                read.csv("all_data/species_data/non_threatened/plants/part3/simple_summary.csv"),
                                read.csv("all_data/species_data/non_threatened/birds/simple_summary.csv"),
                                read.csv("all_data/species_data/non_threatened/herptiles/simple_summary.csv"),
                                read.csv("all_data/species_data/non_threatened/mammals/simple_summary.csv")) %>% 
  rename(SCI_NAME = scientificName, ID_NO = internalTaxonId) %>% # rename columns to align with polygon column names
  dplyr::select(ID_NO, redlistCategory) %>%  # retain relevant columns
  distinct()

# Read in the species data, add columns characterising the taxon and threat status, and
# filter as we did for the threatened data
nt_birds <- rbind(st_read("all_data/species_data/non_threatened/birds/data_0.shp"),
                  st_read("all_data/species_data/non_threatened/birds/data_1.shp"),
                  st_read("all_data/species_data/non_threatened/birds/data_2.shp")) %>%
  st_drop_geometry()
  mutate(taxon = "Birds", status = "Not Threatened") %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, taxon, status) %>% # retain relevant columns
  st_transform(., crs = "+proj=moll") # reproject to Mollweide projection

nt_mammals <- st_read("all_data/species_data/non_threatened/mammals/data_0.shp") %>% 
  st_drop_geometry()
  mutate(taxon = "Mammals", status = "Not Threatened") %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, taxon, status) %>% # retain relevant columns
  st_transform(., crs = "+proj=moll")

nt_herps <- rbind(st_read("all_data/species_data/non_threatened/herptiles/data_0.shp"),
                  st_read("all_data/species_data/non_threatened/herptiles/data_1.shp")) %>% 
  st_drop_geometry()
  mutate(taxon = "Herptiles", status = "Not Threatened") %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, taxon, status) %>% # retain relevant columns
  st_transform(., crs = "+proj=moll") 

nt_plants <- rbind(st_read("all_data/species_data/non_threatened/plants/part1/data_0.shp"),
                   st_read("all_data/species_data/non_threatened/plants/part1/data_1.shp"),
                   st_read("all_data/species_data/non_threatened/plants/part2/data_0.shp"),
                   st_read("all_data/species_data/non_threatened/plants/part3/data_0.shp")) %>% 
  st_drop_geometry()
  mutate(taxon = "Plants", status = "Not Threatened") %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, taxon, status) %>% # retain relevant columns
  st_transform(., crs = "+proj=moll") 

# Apply the function to each taxonomic group and save the output
nt_birds <- prepare_species(nt_birds)
write_sf(nt_birds, "data/output_data/prepared_species_ranges/nt_birds.gpkg")

nt_mammals <- prepare_species(nt_mammals)
write_sf(nt_mammals, "data/output_data/prepared_species_ranges/nt_mammals.gpkg")

nt_herps <- prepare_species(nt_herps)
write_sf(nt_herps, "data/output_data/prepared_species_ranges/nt_herps.gpkg")

nt_plants <- prepare_species(nt_plants)
write_sf(nt_plants, "data/output_data/prepared_species_ranges/nt_plants.gpkg")

################################ ---- subLocRapoport ---- ################################
# Plant species are not comprehensively assessed on the IUCN Red List, and many assessed
# species do not have published ranges. We will create range estimates for these species

# 1) THREATENED PLANTS
# Read in the prepared species ranges and filter to keep plants only
threatened_plants <- st_read("data/output_data/prepared_species_ranges/threatened_species.gpkg") %>% 
  filter(taxon == "Plants")

# Read in the plant summary data
plants_summary <- read.csv("all_data/species_data/threatened/plants/simple_summary.csv") %>% 
  rename(SCI_NAME = scientificName, ID_NO = internalTaxonId) %>% # rename columns to align with polygon column names
  dplyr::select(ID_NO, SCI_NAME, redlistCategory) # retain relevant columns

# Find the species with summary data but without Red List ranges
species_without_ranges <- setdiff(plants_summary$SCI_NAME, threatened_plants$SCI_NAM)

# Read in the species occurrence point data
plants_points <- read.csv("all_data/species_data/threatened/plants/points_data.csv") %>% 
  rename(SCI_NAME = sci_name, ID_NO = id_no) %>% 
  mutate(subspecies = str_squish(subspecies),
         subpop = str_squish(subpop)) %>% 
  mutate(subspecies = if_else(subspecies %in% c("0", "", "<NULL>"), NA_character_, subspecies),
         subpop = if_else(subpop %in% c("0", "", "<NULL>"), NA_character_, subpop)) %>% 
  filter(presence %in% c(1, 2), # match KBA guidelines on presence
         origin %in% c(1, 2, 6), # match KBA guidelines on origin 
         SCI_NAME %in% species_without_ranges) %>% # retain species that have summary data
  dplyr::select(ID_NO, SCI_NAME, latitude, longitude) # retain relevant columns

species_with_point_data <- unique(plants_points$SCI_NAME) # vector of species names that DO NOT have ranges but DO have occurrence data

# Using the subLocRapoport function for species with 3 or more occurrences, and st_buffer 
# for species with 1/2 occurrences, create our own range polygons for species in the 
# species_with_point_data vector

# library(emstreeR)
rap_list <- list()
for (i in seq_along(species_with_point_data)) {
  # Select rows in turn for each species 
  species <- plants_points %>% 
    filter(SCI_NAME == species_with_point_data[i])
  
  # Create a dataframe with the latitude and longitude of each distinct occurrence
  ll <- data.frame(lat = species$latitude,
                   long = species$longitude) %>% 
    distinct()
  
  # If the species has 3 or more distinct occurrences, generate Rapoport ranges
  if (nrow(ll) >= 3) {
    thepoints <- simProjWiz(ll, returnV = "S")
    # barrier_dis <- ComputeMST(thepoints) # for manually assigning a barrier distance in the rapoport
    # barrier_dis <- max(barrier_dis$distance)*0.1
    sfs <- subLocRapoport(thepoints, returnV = "SF") # default buffer and barrier distance
    sfs <- sfs$buffers %>% 
      st_transform(., crs = "+proj=moll") %>% 
      st_as_sf() %>% 
      mutate(ID_NO = species$ID_NO[1],
             SCI_NAME = species$SCI_NAME[1]) %>% 
      rename(geometry = x)
  }
  
  # If the species has only 1/2 distinct occurrences, add a buffer around each point
  else {
    sfs <- species %>% 
      st_as_sf(., coords = c("longitude", "latitude"), crs = "EPSG:4326") %>% # add spatial features to the dataframe 
      st_transform(., crs = "+proj=moll") # reproject to Mollweide
    sfs <- st_buffer(sfs, dist = 10000) # add 10 km buffer
    
    # If there are 2 buffered points for the same species, dissolve their geometries into one
    if (nrow(sfs) > 1) {
      sfs <- sfs %>% 
        group_by(SCI_NAME) %>%
        summarise(
          across(-geometry, first),        # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the taxon and status columns should be the same for each record of the same species
          geometry = st_union(geometry),   # union the geometries for each species
        ) %>% 
        relocate(ID_NO) %>%  # put ID column back in the first position
        st_as_sf()
    }
  }

  rap_list[[i]] <- sfs
}

df <- do.call(rbind, rap_list) %>% 
  mutate(taxon = "Plants", 
         status = "Threatened") %>% 
  left_join(., plants_summary[, -2], by = "ID_NO") %>% 
  relocate(geometry, .after = redlistCategory)

write_sf(df, "data/output_data/prepared_species_ranges/threatened_plants_rapoport.gpkg")

# 2) NON-THREATENED PLANTS
# Run the same analyses but for non-threatened plants
non_threatened_plants <- st_read("data/output_data/prepared_species_ranges/nt_plants.gpkg")

plants_summary <- rbind(read.csv("all_data/species_data/non_threatened/plants/part1/simple_summary.csv"),
                        read.csv("all_data/species_data/non_threatened/plants/part2/simple_summary.csv"),
                        read.csv("all_data/species_data/non_threatened/plants/part3/simple_summary.csv")) %>% 
  rename(SCI_NAME = scientificName, ID_NO = internalTaxonId) %>% # rename columns to align with polygon column names
  dplyr::select(ID_NO, SCI_NAME, redlistCategory) %>%  # retain relevant columns
  distinct()

species_without_ranges <- setdiff(plants_summary$SCI_NAME, non_threatened_plants$SCI_NAM) 

plants_points <- read.csv("all_data/species_data/non_threatened/plants/points_data.csv") %>% 
  rename(SCI_NAME = sci_name, ID_NO = id_no) %>% 
  mutate(subspecies = str_squish(subspecies),
         subpop = str_squish(subpop)) %>% 
  mutate(subspecies = if_else(subspecies %in% c("0", "", "<NULL>"), NA_character_, subspecies),
         subpop = if_else(subpop %in% c("0", "", "<NULL>"), NA_character_, subpop)) %>% 
  filter(presence %in% c(1, 2), # match KBA guidelines on presence
         origin %in% c(1, 2, 6), # match KBA guidelines on origin 
         SCI_NAME %in% species_without_ranges) %>% # retain species that have summary data
  dplyr::select(ID_NO, SCI_NAME, latitude, longitude) # retain relevant columns

species_with_point_data <- unique(plants_points$SCI_NAME) # vector of species names that DO NOT have ranges but DO have occurrence data

rap_list <- list()
for (i in seq_along(species_with_point_data)) {
  # Select rows in turn for each species 
  species <- plants_points %>% 
    filter(SCI_NAME == species_with_point_data[i])
  
  # Create a dataframe with the latitude and longitude of each distinct occurrence
  ll <- data.frame(lat = species$latitude,
                   long = species$longitude) %>% 
    distinct()
  
  # If the species has 3 or more distinct occurrences, generate Rapoport ranges
  if (nrow(ll) >= 3) {
    thepoints <- simProjWiz(ll, returnV = "S")
    # barrier_dis <- ComputeMST(thepoints) # for manually assigning a barrier distance in the rapoport
    # barrier_dis <- max(barrier_dis$distance)*0.1
    sfs <- subLocRapoport(thepoints, returnV = "SF") # default buffer and barrier distance
    sfs <- sfs$buffers %>% 
      st_transform(., crs = "+proj=moll") %>% 
      st_as_sf() %>% 
      mutate(ID_NO = species$ID_NO[1],
             SCI_NAME = species$SCI_NAME[1]) %>% 
      rename(geometry = x)
  }
  
  # If the species has only 1/2 distinct occurrences, add a buffer around each point
  else {
    sfs <- species %>% 
      st_as_sf(., coords = c("longitude", "latitude"), crs = "EPSG:4326") %>% # add spatial features to the dataframe 
      st_transform(., crs = "+proj=moll") # reproject to Mollweide
    sfs <- st_buffer(sfs, dist = 10000) # add 10 km buffer
    
    # If there are 2 buffered points for the same species, dissolve their geometries into one
    if (nrow(sfs) > 1) {
      sfs <- sfs %>% 
        group_by(SCI_NAME) %>%
        summarise(
          across(-geometry, first),        # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the taxon and status columns should be the same for each record of the same species
          geometry = st_union(geometry),   # union the geometries for each species
        ) %>% 
        relocate(ID_NO) %>%  # put ID column back in the first position
        st_as_sf()
    }
  }
  # file_name <- paste0("data/output_data/prepared_species_ranges/", species$SCI_NAME[1], ".shp")
  # write_sf(sfs, file_name)
  
  rap_list[[i]] <- sfs
}

df <- do.call(rbind, rap_list) %>% 
  mutate(taxon = "Plants", 
         status = "Not Threatened") %>% 
  left_join(., plants_summary[, -2], by = "ID_NO") %>% 
  relocate(geometry, .after = redlistCategory)

write_sf(df, "data/output_data/prepared_species_ranges/not_threatened_plants_rapoport.gpkg")
