# Load libraries
library(sf)
library(tidyverse)
library(giscoR)
library(vegan)
library(pairwiseAdonis)
library(adehabitatHR)
library(terra)
library(writexl)
library(readxl)
library(units)
library(gstat)
library(ggpubr)
library(ggpattern)
library(lsmeans)
library(ggh4x)
library(rstatix)
sf_use_s2(FALSE)

############################# ---- Data Preparation ---- ################################# 
# Read in the country names vector and then the country boundaries
countries <- st_read("data/raw_data/country_boundaries/countries.shp")

# Read in the prepared PCA, KBA, IBA and IPA data
pa_data <- st_read("data/output_data/prepared_areas/wdpa_wdoecm_prepared.shp")
osm_pa_data <- st_read("data/output_data/prepared_areas/IND_CHN_osm_prepared.shp") %>% 
  rename("NAME" = "name")
pa_data <- bind_rows(pa_data, osm_pa_data) # combine the PCA data into one dataframe

ipa_data <- st_read("data/output_data/prepared_areas/important_plant_areas_prepared.shp") %>% 
  rename(IUCN_name = IUCN_nm)
kba_data <- st_read("data/output_data/prepared_areas/key_biodiversity_areas_prepared.shp")
iba_data <- st_read("data/output_data/prepared_areas/important_bird_areas_prepared.shp")

# Break the PCA, KBA, IBA and IPA data down into lists, where each element is a country
pa_list <- split(pa_data, pa_data$IUCN_name)
kba_list <- split(kba_data, kba_data$IUCN_name)
iba_list <- split(iba_data, iba_data$IUCN_name)
ipa_list <- split(ipa_data, ipa_data$IUCN_name)

# Read in the sampled points and their coordinates
sampled_points <- read.csv("data/raw_data/spatial_points.csv")
sampled_points <- split(sampled_points, sampled_points$IUCN_name) # split into a list with an element for each country

# Convert the sampled points into sf objects
for (i in 1:nrow(countries)) {
  sampled_points[[i]] <- st_as_sf(sampled_points[[i]], 
                                  coords = c("x", "y"), 
                                  crs = st_crs(countries))
}

# Read in the environmental data
retained_variables <- read.csv("data/output_data/retained_variable_names.csv")[,1]
pattern <- paste(retained_variables, collapse = "|")
env_dat <- list.files(path = "data/output_data/masked_env_data/", pattern = pattern, full.names = TRUE)
env_dat <- lapply(env_dat, rast)
env_dat <- rast(env_dat)

# Plot and view the environmental data
plot(env_dat)

# Extract the environmental data from the sampled points
env_extractions <- list()
for (i in 1:nrow(countries)) {
  dat <- terra::extract(env_dat, sampled_points[[i]], xy = TRUE) %>%
    drop_na()
  env_extractions[[i]] <- dat
}

# Save the environmental extractions in excel
sheet_names <- countries$IUCN_name
named_env_extractions <- setNames(env_extractions, sheet_names)

write_xlsx(named_env_extractions, path = "data/output_data/env_extractions.xlsx",
           col_names = TRUE, format = TRUE)

# Read the environmental extractions back in
env_extractions <- lapply(excel_sheets("data/output_data/env_extractions.xlsx"), 
                          function(sheet) {
                            read_excel("data/output_data/env_extractions.xlsx", sheet = sheet)
                          })

# Standardise the environmental data
standardised_extractions <- list()
for (i in 1:nrow(countries)) {
  # Standardise the data
  dat <- decostand(env_extractions[[i]][, -c(1,15,16)], method = "standardize")
  
  # Add the other columns back in
  dat <- cbind(dat,
               ID = env_extractions[[i]][,1],
               x = env_extractions[[i]][, 15],
               y = env_extractions[[i]][, 16])
  
  standardised_extractions[[i]] <- dat
}

# Save the standardised environmental extractions in excel
sheet_names <- countries$IUCN_name
named_env_extractions <- setNames(standardised_extractions, sheet_names)

write_xlsx(named_env_extractions, path = "data/output_data/standardised_env_extractions.xlsx",
           col_names = TRUE, format = TRUE)

# Read the standardised environmental extractions back in
standardised_extractions <- lapply(excel_sheets("data/output_data/standardised_env_extractions.xlsx"), 
                                   function(sheet) {
                                     read_excel("data/output_data/standardised_env_extractions.xlsx", sheet = sheet)
                                   })

# Function to perform intersection check for all points with each polygon dataset
check_intersections <- function(points_sf, pa_polygons, kba_polygons, iba_polygons, ipa_polygons) {
  # Find intersections using st_intersects (returns logical matrices)
  pa_intersections <- st_intersects(points_sf, pa_polygons, sparse = FALSE)
  kba_intersections <- st_intersects(points_sf, kba_polygons, sparse = FALSE)
  iba_intersections <- st_intersects(points_sf, iba_polygons, sparse = FALSE)
  ipa_intersections <- st_intersects(points_sf, ipa_polygons, sparse = FALSE)
  
  # Convert logical matrices to dfs and add area names
  pa_intersect_points <- points_sf[apply(pa_intersections, 1, any), ] %>%
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "PA") # Add area name
  
  kba_intersect_points <- points_sf[apply(kba_intersections, 1, any), ] %>%
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "KBA") # Add area name
  
  iba_intersect_points <- points_sf[apply(iba_intersections, 1, any), ] %>%
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "IBA") # Add area name
  
  ipa_intersect_points <- points_sf[apply(ipa_intersections, 1, any), ] %>%
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "IPA") # Add area name
  
  # Combine all intersecting points into a single dataframe
  combined_df <- rbind(pa_intersect_points, kba_intersect_points, iba_intersect_points, ipa_intersect_points)
  
  # Return the combined dataframe
  return(combined_df)
}

env_intersections <- list()
for (i in 1:nrow(countries)) {
  # Get environmental extractions df
  points_data <- standardised_extractions[[i]] 
  
  # Get respective country's PA, KBA, IBA and IPA data
  pa_dat <- pa_list[[i]]
  kba_dat <- kba_list[[i]]
  iba_dat <- iba_list[[i]]
  ipa_dat <- ipa_list[[i]]
  
  # Covert the environmental extractions df to an sf object
  points_sf <- st_as_sf(points_data, coords = c("x", "y"), crs = "+proj=moll") 
  
  # Use the function above to find sampled points inside PAs, KBAs, IBAs and IPAs
  intersect_results <- check_intersections(points_sf, pa_dat, kba_dat, iba_dat, ipa_dat)
  
  # Prepare the background environmental extractions df to match the format of intersect_results
  background_df <- points_data %>% 
    mutate(Area = "Background") %>%
    dplyr::select(-x, -y)
  
  # Combine the background_df and intersect_results
  intersect_results <- rbind(background_df, intersect_results)
  
  env_intersections[[i]] <- intersect_results
}

# Save the standardised environmental intersections in excel
sheet_names <- countries$IUCN_name
named_env_intersections <- setNames(env_intersections, sheet_names)

write_xlsx(named_env_intersections, path = "data/output_data/standardised_env_intersections.xlsx",
           col_names = TRUE, format = TRUE)

# Read the standardised environmental intersections back in
env_intersections <- lapply(excel_sheets("data/output_data/standardised_env_intersections.xlsx"), 
                            function(sheet) {
                              read_excel("data/output_data/standardised_env_intersections.xlsx", sheet = sheet)
                            })

######################## ---- Principal Component Analysis ---- ########################## 
# Conduct principal component analysis (PCA) on each country's environmental data
pca_list <- list()
for (i in 1:nrow(countries)) {
  # Conduct PCA on the background points
  background_pca <- subset(env_intersections[[i]], Area == "Background")[,-c(14, 15)] %>% 
    prcomp()
  
  # Get summary and model data from the background_pca
  pca_model <- background_pca
  pca_summary <- summary(background_pca)
  
  # Filter the standardised data to exclude the background points
  others_data <- subset(env_intersections[[i]], Area != "Background")[,-c(14, 15)]
  
  # Conduct PCA on the PCA, KBA, IBA and IPA points
  others_pca <- data.frame(predict(background_pca, newdata = others_data))
  
  # Add the Area column to the others_pca
  others_pca <- cbind(Area = subset(env_intersections[[i]], Area != "Background")[,15],
                      others_pca)
  
  # Add the Area column to the background_pca
  background_pca <- cbind(Area = subset(env_intersections[[i]], Area == "Background")[,15],
                          data.frame(background_pca$x))
  
  # Combine the background and others PCA dataframes into one
  dat <- rbind(background_pca, others_pca)
  
  pca_list[[i]] <- list(PCA_df = dat, PCA_model = pca_model, PCA_summary = pca_summary)
}

# Create a list with the variance explained by PC1 and PC2 for each country
explained_variance_list <- list()
for (i in 1:length(pca_list)) {
  # Find explained variance for PC1 as a percentage
  explained_variance_PC1 <- round(pca_list[[i]]$PCA_summary$importance[2, 1]*100, 1)
  
  # Find explained variance for PC2 as a percentage
  explained_variance_PC2 <- round(pca_list[[i]]$PCA_summary$importance[2, 2]*100, 1)
  
  explained_variance_list[[i]] <- list(PC1 = explained_variance_PC1, PC2 = explained_variance_PC2)
}

# Calculate 99% minimum convex polygons for each country and its conservation areas on the 
# first two principal components
calculate_minimum_convex_polygons <- function(pca_data, area, pct) {
  subset(pca_data, area == Area) %>% 
    dplyr::select(PC1, PC2) %>% 
    SpatialPointsDataFrame(coords = ., 
                           data = .,
                           proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% 
    mcp(., percent = pct)
}

mcp_list <- list()
for (i in 1:nrow(countries)) {
  # Calculate minimum convex polygons for the background environmental space
  background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "Background", pct = 100)
  
  # Calculate minimum convex polygons for each area (PA, KBA, IBA and IPA)
  pa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "PA", pct = 99)
  kba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "KBA", pct = 99)
  iba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "IBA", pct = 99)
  ipa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "IPA", pct = 99)

  # Combine minimum convex polygons into a single dataframe, labelled by area
  dat <- rbind(as.data.frame(background_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                 mutate(Area = "Background"),
               as.data.frame(pa_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                 mutate(Area = "PCA"),
               as.data.frame(kba_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                 mutate(Area = "KBA"),
               as.data.frame(iba_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                 mutate(Area = "IBA"),
               as.data.frame(ipa_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                 mutate(Area = "IPA")) %>% 
    relocate(Area)
  
  mcp_list[[i]] <- dat
}

################################# ---- PCA Plots ---- #################################### 
# Combine the PCA and MCP data into one dataframe per country
ggplot_friendly_list <- list()
for (i in 1:nrow(countries)) {
  # Merge the PCA and MCP data
  dat <- rbind(pca_list[[i]]$PCA_df[,1:3] %>% # retain only the Area, PC1 and PC2 columns
                 mutate(Type = "point"), # assign Type value of "point" to denote PCA points
               mcp_list[[i]] %>% 
                 mutate(Type = "polygon")) # assign Type value of "polygon" to denote MCPs
  
  # Add a country column
  dat <- dat %>% 
    mutate(Country = countries$IUCN_name[i])
  
  ggplot_friendly_list[[i]] <- dat
}

names(ggplot_friendly_list) <- countries$IUCN_name

# Create plots from the PCA and save each country as a separate file for viewing
for (i in 1:nrow(countries)) {
  plot <- ggplot() +
    geom_point(data = subset(ggplot_friendly_list[[i]], Type == "point"),
               mapping = aes(x = PC1, y = PC2, color = Area),
               alpha = 0.1, size = 0.8, shape = 1)+
    geom_polygon(data = subset(ggplot_friendly_list[[i]], Type == "polygon"),
                 mapping = aes(x = PC1, y = PC2, color = Area),
                 alpha = 0, linewidth = 0.7)+
    scale_color_manual(values = c("Background" = "grey", # color palette should be colorblind safe, based on Wong (2011)
                                  "PCA" = "#0072B2",
                                  "KBA" = "#E69F00",
                                  "IPA" = "#009E73",
                                  "IBA" = "#CC79A7")) +
    labs(x = "PC1", y = "PC2", color = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black", size = 10),
          legend.position = "bottom",
          legend.text = element_text(size = 10))
  
  file_name <- paste0("outputs/Environmental_Overlaps/Individual_Countries/", countries$IUCN_name[i], ".png")
  ggsave(filename = file_name, plot = plot, width = 180, height = 180, units = "mm")
    
}

# Create two facet plots for the supplementary materials (including IBAs)
ggplot_friendly_df_a <- do.call(rbind, ggplot_friendly_list[1:15])
ggplot_friendly_df_b <- do.call(rbind, ggplot_friendly_list[16:28])

supp_a <- ggplot() +
  geom_point(data = subset(ggplot_friendly_df_a, Type == "point"),
             mapping = aes(x = PC1, y = PC2, color = Area),
             alpha = 0.1, size = 0.8, shape = 1)+
  geom_polygon(data = subset(ggplot_friendly_df_a, Type == "polygon"),
               mapping = aes(x = PC1, y = PC2, color = Area),
               alpha = 0, linewidth = 0.7)+
  scale_color_manual(values = c("Background" = "grey", # color palette should be colorblind safe, based on Wong (2011)
                                "PCA" = "#0072B2",
                                "KBA" = "#E69F00",
                                "IPA" = "#009E73",
                                "IBA" = "#CC79A7")) +
  labs(x = "PC1", y = "PC2", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3, scales = "free")

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/supp_pca_a.png",
       plot = supp_a,
       width = 180, height = 297, units = "mm",
       dpi = 600)

supp_b <- ggplot() +
  geom_point(data = subset(ggplot_friendly_df_b, Type == "point"),
             mapping = aes(x = PC1, y = PC2, color = Area),
             alpha = 0.1, size = 0.8, shape = 1)+
  geom_polygon(data = subset(ggplot_friendly_df_b, Type == "polygon"),
               mapping = aes(x = PC1, y = PC2, color = Area),
               alpha = 0, linewidth = 0.7)+
  scale_color_manual(values = c("Background" = "grey", # color palette should be colorblind safe, based on Wong (2011)
                                "PCA" = "#0072B2",
                                "KBA" = "#E69F00",
                                "IPA" = "#009E73",
                                "IBA" = "#CC79A7")) +
  labs(x = "PC1", y = "PC2", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3, scales = "free")

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/supp_pca_b.png",
       plot = supp_b,
       width = 180, height = 240, units = "mm",
       dpi = 600)

# Create part B of the main figure plot (PCAs of four countries without IBAs)
ggplot_friendly_df_c <- do.call(rbind, ggplot_friendly_list) %>% 
  filter(Country %in% c("United Kingdom", "Algeria", "Libya", "Türkiye"))

names(explained_variance_list) <- countries$IUCN_name

country_facet <- ggplot()+
  geom_point(data = subset(ggplot_friendly_df_c, Type == "point" & Area != "IBA"),
             mapping = aes(x = PC1, y = PC2, color = Area),
             alpha = 0.1, size = 0.8, shape = 1)+
  geom_polygon(data = subset(ggplot_friendly_df_c, Type == "polygon" & Area != "IBA"),
               mapping = aes(x = PC1, y = PC2, color = Area),
               alpha = 0, linewidth = 0.7)+
  scale_color_manual(values = c("Background" = "grey", # color palette should be colorblind safe, based on Wong (2011)
                                "PCA" = "#0072B2",
                                "KBA" = "#E69F00",
                                "IPA" = "#009E73")) +
  labs(x = "PC1", y = "PC2", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~ Country, scales = "free", ncol = 2,
           labeller = labeller(Country = c("Algeria" = "Algeria (PC1 = 29.2%; PC2 = 19.3%)", 
                                           "Libya" = "Libya (PC1 = 36.5%; PC2 = 14.2%)", 
                                           "Türkiye" = "Türkiye (PC1 = 28.9%; PC2 = 16.5%)", 
                                           "United Kingdom" = "United Kingdom (PC1 = 44.7%; PC2 = 12.0%)")))

ggsave(filename = "outputs/Environmental_Overlaps/Main_Figs/country_facet.png",
       plot = country_facet,
       height = 180, width = 180, units = "mm",
       dpi = 600)

################################## ---- MCP Area ---- #################################### 
# Add an extra line of code to the calculate_mcp function, to turn the MCP to a spatial object
calculate_minimum_convex_polygons <- function(pca_data, area, pct) {
  subset(pca_data, area == Area) %>% 
    dplyr::select(PC1, PC2) %>% 
    SpatialPointsDataFrame(coords = ., 
                           data = .,
                           proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% 
    mcp(., percent = pct) %>% 
    st_as_sf() # add sf to the MCP for this analysis
  
}

# Helper function to safely get area or assign 0 if empty, i.e. without an error
safe_area <- function(x) {
  if (nrow(x) == 0) {
    return(0)
  } else {
    return(as.numeric(st_area(x)))
  }
}
background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[2]]$PCA_df, area = "Background", pct = 100)

# Apply the calculate_minimum_convex_polygons function to calculate the area of each country's 
# various minimum convex polygons
mcp_area <- list()

for (i in 1:nrow(countries)) {
  # Calculate MCPs for country's background environmental space
  background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "Background", pct = 100)
  
  # Calculate MCPs for PAs, KBAs and IPAs
  pa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "PA", pct = 99)
  kba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "KBA", pct = 99)
  ipa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, area = "IPA", pct = 99)
  
  # Calculate MCPs for pairwise overlaps of PAs, KBAs and IPAs
  pa_kba_mcp <- st_intersection(pa_mcp, kba_mcp)
  pa_ipa_mcp <- st_intersection(pa_mcp, ipa_mcp)
  kba_ipa_mcp <- st_intersection(kba_mcp, ipa_mcp)
  
  # Calculate MCPs for KBAs and IPAs outside of PA's environmental space
  kba_not_pa <- st_difference(kba_mcp, pa_mcp)
  ipa_not_pa <- st_difference(ipa_mcp, pa_mcp)
  
  # Calculate the MCP for PAs, KBAs and IPAs combined (i.e. all networks)
  all <- st_union(st_union(pa_mcp, kba_mcp),
                  ipa_mcp)
  
  # Calculate area of each MCP and then convert area to proportion of background environmental space
  dat <- data.frame(country = countries$IUCN_name[i],
                    designation = c("Background", "PCA", "KBA", "IPA", 
                                    "PCA_KBA", "PCA_IPA", "KBA_IPA",
                                    "KBA_not_PCA", "IPA_not_PCA",
                                    "All"),
                    area = rbind(safe_area(background_mcp), 
                                 safe_area(pa_mcp), 
                                 safe_area(kba_mcp), 
                                 safe_area(ipa_mcp),
                                 safe_area(pa_kba_mcp),
                                 safe_area(pa_ipa_mcp),
                                 safe_area(kba_ipa_mcp),
                                 safe_area(kba_not_pa),
                                 safe_area(ipa_not_pa),
                                 safe_area(all))) %>% 
    mutate(
      background_area = area[designation == "Background"],
      area_pct = (area / background_area) * 100
    ) %>%
    dplyr::select(-area, -background_area)
  
  mcp_area[[i]] <- dat
}

# # Plotting to visualise these calculations and overlaps
# # Calculate MCPs for country's background environmental space
# background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[4]]$PCA_df, area = "Background", pct = 100)
# 
# # Calculate MCPs for PAs, KBAs and IPAs
# pa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[4]]$PCA_df, area = "PA", pct = 99)
# kba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[4]]$PCA_df, area = "KBA", pct = 99)
# ipa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[4]]$PCA_df, area = "IPA", pct = 99)
# 
# # Calculate MCPs for pairwise overlaps of PAs, KBAs and IPAs
# pa_kba_mcp <- st_intersection(pa_mcp, kba_mcp)
# pa_ipa_mcp <- st_intersection(pa_mcp, ipa_mcp)
# kba_ipa_mcp <- st_intersection(kba_mcp, ipa_mcp)
# 
# # Calculate MCPs for KBAs and IPAs outside of PA's environmental space
# kba_not_pa <- st_difference(kba_mcp, pa_mcp)
# ipa_not_pa <- st_difference(ipa_mcp, pa_mcp)
# 
# # Calculate the MCP for PAs, KBAs and IPAs combined (i.e. all networks)
# all <- st_union(st_union(pa_mcp, kba_mcp),
#                 ipa_mcp)
# 
# plot(background_mcp %>% as_Spatial())
# plot(pa_mcp %>% as_Spatial(), add = T, border = "#0072B2")
# plot(kba_mcp %>% as_Spatial(), add = T, border = "#E69F00")
# plot(ipa_mcp %>% as_Spatial(), add = T, border = "#009E73")
# 
# plot(ipa_not_pa %>% as_Spatial(), add = T, col = "black")

# Combine into a dataframe and pivot wider to have a column for each area pct category
mcp_area_df <- do.call(rbind, mcp_area) %>% 
  pivot_wider(names_from = designation, values_from = area_pct) 

# Find mean across countries and add to dataframe
means <- c("Mean", as.numeric(colMeans(mcp_area_df[, -1])))
mcp_area_df <- rbind(mcp_area_df, means) %>% 
  dplyr::select(-Background)

# Make values numeric and round to 2 decimal places
mcp_area_df <- mcp_area_df %>% 
  mutate(across(2:10, ~ as.numeric(.))) %>% 
  mutate(across(2:10, ~ round(., 2))) 

# # Save this as a csv
# write.csv(mcp_area_df,
#           file = "data/output_data/env_space_mcp_area.csv",
#           row.names = F)

# Read in this environmental space MCP area data
mcp_area_df <- read.csv("data/output_data/env_space_mcp_area.csv") %>% 
  pivot_longer(PCA:All, names_to = "area_category", values_to = "pct") %>% 
  mutate(global = "Global")

global_boxplot <- ggplot(data = mcp_area_df,
       mapping = aes(x = pct,
                     y = factor(area_category,
                                level = c("IPA_not_PCA",
                                          "KBA_not_PCA",
                                          "KBA_IPA",
                                          "PCA_IPA",
                                          "PCA_KBA",
                                          "IPA",
                                          "KBA",
                                          "PCA",
                                          "All"),
                                labels = c("IPA (not PCA)",
                                           "KBA (not PCA)",
                                           "KBA & IPA",
                                           "PCA & IPA",
                                           "PCA & KBA",
                                           "IPA",
                                           "KBA",
                                           "PCA",
                                           "All"))))+
  geom_boxplot_pattern(mapping = aes(fill = area_category,
                                     pattern = area_category,
                                     pattern_fill = area_category),
                       outliers = FALSE, show.legend = FALSE, pattern_density = 0.5, color = "black")+
  geom_hline(yintercept = 5.5, linetype = "dashed", colour = "grey")+
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "grey")+
  geom_hline(yintercept = 8.5, linetype = "dashed", colour = "grey")+
  geom_jitter(height = 0.25, size = 0.7, show.legend = F)+
  scale_x_continuous(expand = c(0,0.5,0,5))+
  scale_fill_manual(values = c(
    "All" = "#666666",
    "PCA" = "#0072B2",
    "KBA" = "#E69F00",
    "IPA" = "#009E73",
    "PCA_KBA" = "#0072B2",
    "PCA_IPA" = "#0072B2",
    "KBA_IPA" = "#E69F00",
    "KBA_not_PCA" = "#E69F00",
    "IPA_not_PCA" = "#009E73"
  )) +
  scale_pattern_manual(values = c(
    "All" = "none",
    "PCA" = "none",
    "KBA" = "none",
    "IPA" = "none",
    "PCA_KBA" = "stripe",
    "PCA_IPA" = "stripe",
    "KBA_IPA" = "stripe",
    "KBA_not_PCA" = "none",
    "IPA_not_PCA" = "none"
  )) +
  scale_pattern_fill_manual(values = c(
    "All" = NA,
    "PCA" = NA,
    "KBA" = NA,
    "IPA" = NA,
    "PCA_KBA" = "#E69F00",
    "PCA_IPA" = "#009E73",
    "KBA_IPA" = "#009E73",
    "KBA_not_PCA" = NA,
    "IPA_not_PCA" = NA
  )) +
  labs(x = "Coverage of background environmental space (%)",
       y = NULL)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10))+
  facet_wrap(~global)

ggsave(filename = "outputs/Environmental_Overlaps/Main_Figs/global_boxplot.png",
       plot = global_boxplot,
       height = 140, width = 180, units = "mm",
       dpi = 600)             

global_boxplot_and_country_facet <- ggarrange(global_boxplot, country_facet, 
                                              ncol = 1,
                                              labels = c("A", "B"),
                                              hjust = -1,
                                              heights = c(1, 1.25))

ggsave(filename = "outputs/Environmental_Overlaps/Main_Figs/global_boxplot_and_country_facet.png",
       plot = global_boxplot_and_country_facet,
       width = 180, height = 240, units = "mm",
       dpi = 600)

################################ ---- PCA Loadings ---- ################################## 
# Create a vector of real and simple variable names
real_variable_names <- c("Aspect", "Clay content", "Nitrogen content", "Sand content",
                           "Slope", "Topographic position index", "Precipitation seasonality",
                           "Precipitation of coldest ¼", "Mean diurnal range",
                           "Isothermality", "Mean temperature of wettest ¼", "Mean temperature of driest ¼",
                           "Elevation")

simple_variable_names <- c("Aspect", "Clay", "Nitrogen", "Sand", "Slope", "TPI", "Bio15",
                           "Bio19", "Bio2", "Bio3", "Bio8", "Bio9", "Elevation")

loading_list <- list()
for (i in 1:nrow(countries)) {
  # Get variable loadings for PC1 and PC2
  loading <- pca_list[[i]]$PCA_model$rotation[,1:2] %>% 
    as.data.frame() %>% 
    mutate(Variable = rownames(.), # add variable column
           real_variable = real_variable_names,
           simple_variable = simple_variable_names, # add simple variable name column
           Country = countries$IUCN_name[i]) # add country column
  
  # Convert to long format
  loading_longer <- loading %>% 
    pivot_longer(PC1:PC2, names_to = "PC", values_to = "Loading") 
  
  loading_list[[i]] <- loading_longer
}

# Combine scree_list into one dataframe
loading_dat_a <- do.call(rbind, loading_list[1:15]) # data for first figure panel
loading_dat_b <- do.call(rbind, loading_list[16:length(loading_list)]) # data for second figure panel

# Create the first loading plot panel
loading_a <- ggplot(data = loading_dat_a,
       mapping = aes(x = Loading, y = simple_variable, fill = PC))+
  geom_col(position = "dodge", width = 0.6, color = "black", linewidth = 0.2)+
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             color = "darkgray")+
  labs(x =  "PCA loading value",
       y = NULL,
       fill = NULL)+
  scale_fill_manual(values = c("PC1" = "#ffdb00", "PC2" = "#400b0b"))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3)

# Create the second loading plot panel
loading_b <- ggplot(data = loading_dat_b,
       mapping = aes(x = Loading, y = simple_variable, fill = PC))+
  geom_col(position = "dodge", width = 0.6, color = "black", linewidth = 0.2)+
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             color = "darkgray")+
  labs(x =  "PCA loading value",
       y = NULL,
       fill = NULL)+
  scale_fill_manual(values = c("PC1" = "#ffdb00", "PC2" = "#400b0b"))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3)

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/pca_loadings_a.png",
       plot = loading_a,
       width = 180, height = 297, units = "mm",
       dpi = 600)

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/pca_loadings_b.png",
       plot = loading_b,
       width = 180, height = 240, units = "mm",
       dpi = 600)

################################ ---- Scree Plots ---- ################################### 
# Get eigenvalues and variance explained by each principal component in each country's PCA
scree_list <- list()
for (i in 1:nrow(countries)) {
  # Calculate eigenvalues and explained variance
  eigenvalues <- pca_list[[i]]$PCA_model$sdev^2
  variance_explained <- (eigenvalues / sum(eigenvalues)) * 100
  
  # Calculate cumulative variance explained
  cumulative_variance_explained <- cumsum(variance_explained)
  
  # Create a data frame with PC, eigenvalue, variance explained and cumulative variance explained
  dat <- data.frame(
    PC = factor(1:length(eigenvalues)),  # Principal component numbers
    Eigenvalue = eigenvalues,            # Eigenvalues
    Variance = variance_explained,       # Variance explained by each PC
    CumulativeVariance = cumulative_variance_explained, # Cumulative variance explained
    Country = countries$IUCN_name[i]     # Country name
  )
  
  # Store the data frame in the list
  scree_list[[i]] <- dat
}
scree_eig <- do.call(rbind, scree_list) %>% 
  filter(PC == 3)
# Combine scree_list into one dataframe
scree_dat_a <- do.call(rbind, scree_list[1:15]) # data for first figure panel
scree_dat_b <- do.call(rbind, scree_list[16:length(scree_list)]) # data for second figure panel

# Create the first scree plot panel
scree_a <- ggplot(scree_dat_a, aes(x = PC)) +
  geom_line(aes(y = Eigenvalue, group = 1, color = "Eigenvalue"), linewidth = 0.7) +
  geom_point(aes(y = Eigenvalue, color = "Eigenvalue")) +
  geom_line(aes(y = CumulativeVariance / 10, group = 1, color = "Cumulative Variance Explained"), 
            linewidth = 0.7) +  
  geom_hline(yintercept = 1, color = "black", linewidth = 0.5, linetype = "dashed")+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 1),
    sec.axis = sec_axis(~. * 10, name = "Cumulative variance explained (%)")  
  ) +
  scale_color_manual(values = c("Eigenvalue" = "#ffdb00", "Cumulative Variance Explained" = "#400b0b"),
                     labels = c("Eigenvalue" = "Eigenvalue", "Cumulative Variance Explained" = "Cumulative variance explained")) +
  labs(x = "Principal component",
       y = "Eigenvalue",
       color = NULL)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3)

# Create the second scree plot panel
scree_b <- ggplot(scree_dat_b, aes(x = PC)) +
  geom_line(aes(y = Eigenvalue, group = 1, color = "Eigenvalue"), linewidth = 0.7) +
  geom_point(aes(y = Eigenvalue, color = "Eigenvalue")) +
  geom_line(aes(y = CumulativeVariance / 10, group = 1, color = "Cumulative Variance Explained"), 
            linewidth = 0.7) +  
  geom_hline(yintercept = 1, color = "black", linewidth = 0.5, linetype = "dashed")+
  scale_y_continuous(
    labels = scales::number_format(accuracy = 1),
    sec.axis = sec_axis(~. * 10, name = "Cumulative variance explained (%)")  
  ) +
  scale_color_manual(values = c("Eigenvalue" = "#ffdb00", "Cumulative Variance Explained" = "#400b0b"),
                     labels = c("Eigenvalue" = "Eigenvalue", "Cumulative Variance Explained" = "Cumulative variance explained")) +
  labs(x = "Principal component",
       y = "Eigenvalue",
       color = NULL)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~Country, ncol = 3)

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/scree_a.png",
       plot = scree_a,
       width = 180, height = 297, units = "mm",
       dpi = 600)

ggsave(filename = "outputs/Environmental_Overlaps/Supplementary_Materials/scree_b.png",
       plot = scree_b,
       width = 180, height = 240, units = "mm",
       dpi = 600)

########################## ---- Spatial Autocorrelation ---- ############################# 
# 1) EACH COUNTRY
## Function to buffer points in XY space:
## Returns the original data table with buffered points removed.
## Runs numerous iterations, as the random point selection can result in more/fewer output points.
# foo - a data.frame to select from with columns x, y
# buffer - the minimum distance between output points
# reps - the number of repetitions for the points selection
buffer.f <- function(foo, buffer, reps){
  # Make list of suitable vectors
  suitable <- list()
  for(k in 1:reps){
    # Make the output vector
    outvec <- as.numeric(c())
    # Make the vector of dropped (buffered out) points
    dropvec <- c()
    for(i in 1:nrow(foo)){
      # Stop running when all points exhausted
      if(length(dropvec)<nrow(foo)){
        # Set the rows to sample from
        if(i>1){
          rowsleft <- (1:nrow(foo))[-c(dropvec)]
        } else {
          rowsleft <- 1:nrow(foo)
        }
        # Randomly select point
        outpoint <- as.numeric(sample(as.character(rowsleft),1))
        outvec[i] <- outpoint
        # Remove points within buffer
        outcoord <- foo[outpoint,c("x","y")]
        dropvec <- c(dropvec, which(sqrt((foo$x-outcoord$x)^2 + (foo$y-outcoord$y)^2)<buffer))
        # Remove unnecessary duplicates in the buffered points
        dropvec <- dropvec[!duplicated(dropvec)]
      } 
    } 
    # Populate the suitable points list
    suitable[[k]] <- outvec
  }
  # Go through the iterations and pick a list with the most data
  best <- unlist(suitable[which.max(lapply(suitable,length))])
  foo[best,]
}

thinned_env_extractions <- list()

for (i in 1:nrow(countries)) {
  # Select the country's environmental data and keep the x and y columns only
  env_dat <- standardised_extractions[[i]] %>% 
    dplyr::select(x, y)
  
  # Apply the spatial thinning function with a buffer of 10km and 5 repetitions
  out <- buffer.f(foo = env_dat, buffer = 10000, reps = 5)
  
  # Turn the output into a spatial object
  out <- st_as_sf(out, coords = c("x", "y"), crs = "+proj=moll")
  
  # Turn the environmental data into a spatial object and filter to only keep rows with
  # the retained geometries
  env_dat <- standardised_extractions[[i]] %>% 
    st_as_sf(., coords = c("x", "y"), crs = "+proj=moll") %>% 
    filter(geometry %in% out$geometry)
  
  # Find the points that intersect with PCAs, KBAs and IPAs
  pa_intersections <- st_intersects(env_dat, pa_list[[i]], sparse = FALSE)
  pa_intersections <- env_dat[apply(pa_intersections, 1, any), ] %>%
    mutate(Area = "PCA") # Add area name
  
  # Make sure at least one PCA has been included
  if (nrow(pa_intersections) == 0) {
    env_dat <- standardised_extractions[[i]] %>% 
      st_as_sf(., coords = c("x", "y"), crs = "+proj=moll")
    
    pa_intersections <- st_intersects(env_dat, pa_list[[i]], sparse = FALSE)
    pa_intersections <- env_dat[apply(pa_intersections, 1, any), ] %>%
      mutate(Area = "PCA") %>% # Add area name 
      sample_n(size = 1) # take one row 
  }
  
  kba_intersections <- st_intersects(env_dat, kba_list[[i]], sparse = FALSE)
  kba_intersections <- env_dat[apply(kba_intersections, 1, any), ] %>%
    mutate(Area = "KBA") # Add area name
  
  # Make sure at least one KBA has been included
  if (nrow(kba_intersections) == 0) {
    env_dat <- standardised_extractions[[i]] %>% 
      st_as_sf(., coords = c("x", "y"), crs = "+proj=moll")
    
    kba_intersections <- st_intersects(env_dat, kba_list[[i]], sparse = FALSE)
    kba_intersections <- env_dat[apply(kba_intersections, 1, any), ] %>%
      mutate(Area = "KBA") %>% # Add area name 
      sample_n(size = 1) # take one row 
  }
  
  ipa_intersections <- st_intersects(env_dat, ipa_list[[i]], sparse = FALSE)
  ipa_intersections <- env_dat[apply(ipa_intersections, 1, any), ] %>%
    mutate(Area = "IPA") # Add area name
  
  # Make sure at least one IPA has been included
  if (nrow(ipa_intersections) == 0) {
    env_dat <- standardised_extractions[[i]] %>% 
      st_as_sf(., coords = c("x", "y"), crs = "+proj=moll")
    
    ipa_intersections <- st_intersects(env_dat, ipa_list[[i]], sparse = FALSE)
    ipa_intersections <- env_dat[apply(ipa_intersections, 1, any), ] %>%
      mutate(Area = "IPA") %>% # Add area name 
      sample_n(size = 1) # take one row 
  }
  
  thinned_env_extractions[[i]] <- rbind(pa_intersections,
                                        kba_intersections,
                                        ipa_intersections)
}

# Save the thinned standardised environmental intersections in excel
sheet_names <- countries$IUCN_name
named_thinned_intersections <- setNames(thinned_env_extractions, sheet_names)
named_thinned_intersections <- lapply(named_thinned_intersections, function(df) {
  df %>%
    mutate(
      x = st_coordinates(.)[, 1], # Extract X coordinates
      y = st_coordinates(.)[, 2]  # Extract Y coordinates
    ) %>%
    st_drop_geometry() # Optionally drop the geometry column if it's no longer needed
})
    
write_xlsx(named_thinned_intersections, path = "data/output_data/thinned_env_extractions.xlsx",
           col_names = TRUE, format = TRUE)

# Read the thinned standardised environmental intersections back in
thinned_env_extractions <- lapply(excel_sheets("data/output_data/thinned_env_extractions.xlsx"), 
                                  function(sheet) {
                                    # Read the sheet
                                    df <- read_excel("data/output_data/thinned_env_extractions.xlsx", sheet = sheet)
                                    
                                    # Convert x and y to an sf object with POINT geometry
                                    st_as_sf(df, coords = c("x", "y"), crs = "+proj=moll") 
                                  })

# Create variograms to assess whether the spatial thinning has reduced spatial autocorrelation
variogram_list <- list()
for (i in 1:nrow(countries)) {
  intersect_results <- thinned_env_extractions[[i]] %>% 
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  
  coordinates(intersect_results) = ~x + y
  intersect_results$Area <- as.factor(intersect_results$Area)
  
  output <- variogram(Area ~ x+y, intersect_results) %>% 
    mutate(country = countries$IUCN_name[i])
  
  variogram_list[[i]] <- output
}

# Combine into one dataframe
variograms <- do.call(rbind, variogram_list)

# Plot the variograms
variograms <- ggplot(data = variograms, aes(x = dist, y = gamma))+
  geom_point(size = 1)+
  labs(x = "Distance (m)", y = "Semivariance")+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, by = 0.5))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 500000), breaks = seq(0, 500000, by = 250000))+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~country, ncol = 4)

# Save the variograms
ggsave("outputs/Environmental_Overlaps/Supplementary_Materials/Semi-Variograms.png",
       plot = variograms,
       height = 297, width = 180, units = "mm",
       dpi = 600)

################################# ---- PERMANOVA ---- #################################### 
adonis_list <- list()

for (i in 1:nrow(countries)) {
  # Dataframe with Area only
  area_grouping <- data.frame(Area = thinned_env_extractions[[i]]$Area) %>% 
    mutate(Area = as.factor(Area)) %>% 
    st_drop_geometry()
  
  # Remove Area column from environmental data df 
  env_dat <- thinned_env_extractions[[i]] %>% 
    dplyr::select(-c(Area, ID)) %>%
    st_drop_geometry()
  
  # Run PERMANOVA on the subset data
  adonis_out <- adonis2(
    env_dat ~ Area, 
    data = area_grouping, 
    method = "eu",
    permutations = 999
  )
  
  # Run pairwise adonis
  pairwise_out <- pairwise.adonis(
    x = env_dat, 
    factors = area_grouping$Area,
    sim.method = "euclidean", 
    p.adjust.m = "bonferroni",
    perm = 999
  )
  
  adonis_list[[i]] <- list(adonis = adonis_out, pairwise = pairwise_out)
}

# Function to extract the F and p values
extract_values <- function(x, country_name) {
  data.frame(
    Country = country_name,
    F_value = x$adonis$F[1],
    p_value = x$adonis$`Pr(>F)`[1],
    PCA_KBA_F = x$pairwise$F.Model[1],
    PCA_KBA_p = x$pairwise$p.adjusted[1],
    PCA_IPA_F = x$pairwise$F.Model[2],
    PCA_IPA_p = x$pairwise$p.adjusted[2],
    KBA_IPA_F = x$pairwise$F.Model[3],
    KBA_IPA_p = x$pairwise$p.adjusted[3]
  )
}

# Apply the extraction function to each element in the list
permanova_results <- mapply(extract_values, adonis_list, countries$IUCN_name, SIMPLIFY = FALSE)
permanova_results <- do.call(rbind, permanova_results)

write.csv(permanova_results, "data/output_data/permanova_output.csv", row.names = FALSE)

# 2) GLOBAL
all_env_intersections <- list()

for (i in 1:nrow(countries)) {
  env_dat <- st_as_sf(standardised_extractions[[i]], coords = c("x", "y"), crs = "+proj=moll")
  
  # Find the points that intersect with PCAs, KBAs and IPAs
  pa_intersections <- st_intersects(env_dat, pa_list[[i]], sparse = FALSE)
  pa_intersections <- env_dat[apply(pa_intersections, 1, any), ] %>%
    mutate(Area = "PCA") # Add area name
  
  kba_intersections <- st_intersects(env_dat, kba_list[[i]], sparse = FALSE)
  kba_intersections <- env_dat[apply(kba_intersections, 1, any), ] %>%
    mutate(Area = "KBA") # Add area name
  
  ipa_intersections <- st_intersects(env_dat, ipa_list[[i]], sparse = FALSE)
  ipa_intersections <- env_dat[apply(ipa_intersections, 1, any), ] %>%
    mutate(Area = "IPA") # Add area name
  
  all_env_intersections[[i]] <- rbind(pa_intersections,
                                      kba_intersections,
                                      ipa_intersections) %>% 
    mutate(country = countries$IUCN_name[i]) %>% 
    st_drop_geometry()
}

# Combine the standardised environmental intersections into one dataframe
all_env_intersections_df <- do.call(rbind, all_env_intersections)

# Take a sub-sample of the stratified sample, maintaining original data structure, i.e.
# with the same proportion of points between countries and the same proportion of points
# per area (PCA:KBA:IPA) within each country
total_samples <- 10000 # define sub-sample size

sampled_env_intersections <- all_env_intersections_df %>% 
  group_by(country, Area) %>% 
  mutate(group_size = n()) %>% 
  sample_frac(total_samples / nrow(.), weight = group_size) %>% 
  ungroup()

# Conduct one global PERMANOVA 
env_dat <- sampled_env_intersections[, 1:13] %>% as.matrix() # select the 13 environmental variables and convert to matrix

area_grouping <- data.frame(Area = sampled_env_intersections$Area, # select the area and country variables and convert them to factors
                            country = sampled_env_intersections$country) %>% 
  mutate(Area = as.factor(Area),
         country = as.factor(country))

# Run the PERMANOVA with permutations free across the dataset, i.e. environmental variation across area without controlling for country-level differences
global_permanova_out <- adonis2(
  env_dat ~ Area, 
  data = area_grouping, 
  method = "eu",
  permutations = 999
)

global_pairwise_permanova_out <- pairwise.adonis(
  x = env_dat, factors = area_grouping$Area,
  sim.method = "euclidean", 
  p.adjust.m = "bonferroni",
  perm = 999
)

# Run the PERMANOVA with permutations restricted within each country, i.e. environmental variation across area after controlling for country-level differences 
global_permanova_out_constrained <- adonis2(
  env_dat ~ Area, 
  data = area_grouping, 
  method = "eu", 
  permutations = how(blocks = area_grouping$country)  # Constrain permutations - usually I would use the strata = argument, but it doesn't work/change anything here (I think due to countries with too few samples), so this is an explicit alternative
)

# NB: pairwise.adonis does not support constrained permutations

# Read in each country's PERMANOVA outputs
national_permanova_out <- read.csv("data/output_data/permanova_output.csv")

# Add the global PERMANOVA outputs to the national ones and overwrite the file
all_permanova_out <- rbind(
  national_permanova_out,
  c("Global", global_permanova_out$F[1], global_permanova_out$`Pr(>F)`[1], 
    global_pairwise_permanova_out$F.Model[3], global_pairwise_permanova_out$p.adjusted[3],
    global_pairwise_permanova_out$F.Model[2], global_pairwise_permanova_out$p.adjusted[2],
    global_pairwise_permanova_out$F.Model[1], global_pairwise_permanova_out$p.adjusted[1]),
  c("Global_constrained", global_permanova_out_constrained$F[1], global_permanova_out_constrained$`Pr(>F)`[1],
    rep(NA, 6))
)

write.csv(all_permanova_out, "data/output_data/permanova_output.csv", row.names = FALSE)

############################# ---- Environmental Differentiation ---- ################################ 
# Read in the raw environmental extractions for each country
env_dat <- "Downloads/wetransfer_env_extractions-xlsx_2025-03-04_1653/env_extractions.xlsx"
env_dat <- lapply(excel_sheets(env_dat), function(sheet){
  read_excel(env_dat, sheet = sheet)
})

# Read in the standardised environmental intersection data
intersections <- "Downloads/standardised_env_intersections.xlsx"
intersections <- lapply(excel_sheets(intersections), function(sheet){
  read_excel(intersections, sheet = sheet)
})

# We want to use the raw environmental data but add a column showing whether each point relates to 
# PCAs, KBAs or IPAs (or the background).
for (i in 1:length(env_dat)) {
  env_dat[[i]] <- left_join(intersections[[i]][, 14:15], env_dat[[i]]) %>% 
    filter(Area != "IBA") %>% # we don't need the IBA data here
    dplyr::select(-ID, -x, -y) # only retain the Area column and the environmental variables
}

countries <- c(unique(country_results$country), "Virgin Islands, British")

# Assign simple names to each variable
for (i in 1:length(env_dat)) {
  env_dat[[i]] <- env_dat[[i]] %>% 
    rename(Area = Area, Aspect = aspect, Clay = `clay_5-15cm`, Nitrogen = `nitrogen_5-15cm`, Sand = `sand_5-15cm`,
           Slope = slope, TPI = TPI, Bio15 = wc2.1_30s_bio_15, Bio19 = wc2.1_30s_bio_19, Bio2 = wc2.1_30s_bio_2,
           Bio3 = wc2.1_30s_bio_3, Bio8 = wc2.1_30s_bio_8, Bio9 = wc2.1_30s_bio_9, Elevation = wc2.1_30s_elev) %>%
    mutate(Area = if_else(Area == "PA", "PCA", Area)) %>% # change PA to PCA
    mutate(Country = countries[i])
}

# Density plots
density_plot <- do.call(rbind, env_dat) %>% 
  dplyr::select(-Country) %>% 
  mutate(Aspect = cos(Aspect * pi / 180)) %>% 
  rename("cosine(Aspect)" = Aspect) %>% 
  mutate(Area = factor(Area, levels = c("Background", "PCA", "KBA", "IPA"))) %>%  
  pivot_longer(cols = -Area, names_to = "variable", values_to = "value") %>%
  arrange(variable, Area) %>% 
  mutate(value = round(value, 4)) %>% 
  ggplot(aes(x = value, fill = Area, color = Area)) +
  geom_density(alpha = 0.05, linewidth = 0.7) +
  scale_fill_manual(values = c("Background" = "grey",
                               "PCA" = "#0072B2",
                               "KBA" = "#E69F00",
                               "IPA" = "#009E73"))+
  scale_color_manual(values = c("Background" = "grey",
                                "PCA" = "#0072B2",
                                "KBA" = "#E69F00",
                                "IPA" = "#009E73"))+
  facet_wrap(~ variable, scales = "free", ncol = 3)+
  labs(x = "Value",
       y = "Density",
       fill = NULL,
       color = NULL) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank())

density_plot <- density_plot+
  facetted_pos_scales(
    x = list(
      variable == "TPI" ~ scale_x_continuous(limits = c(-5, 5)),
      variable == "Bio19" ~ scale_x_continuous(limits = c(0, 100)),
      variable == "Slope" ~ scale_x_continuous(limits = c(0, 2.5)),
      variable == "Nitrogen" ~ scale_x_continuous(limits = c(0, 5))
    )
  )

plot(density_plot)
ggsave("OneDrive - University of Cambridge/External_Projects/Kew/EnvDifferentiation_DensityPlot.png",
       density_plot, 
       width = 180, height = 260, units = "mm")

env_geo <- readRDS("Downloads/geo_vs_env_space_slopes.rds")

fig <- ggarrange(env_geo, density_plot,
                 ncol = 1,
                 labels = c("A", "B"),
                 hjust = -1,
                 heights = c(1, 1.75))

ggsave(plot = fig, "OneDrive - University of Cambridge/External_Projects/Kew/slopes_and_densityplot.png",
       width = 180, height = 280, units = "mm", dpi = 600, device = "png")

######################### ---- National Environmental Differentiation ---- #########################
# env_dat[[2]] %>% 
#   pivot_longer(cols = -Area, names_to = "variable", values_to = "value") %>%
#   arrange(variable, Area) %>% 
#   mutate(value = round(value, 4)) %>% 
#   ggplot(aes(x = value, y = variable, fill = Area)) +
#   geom_boxplot(outliers = FALSE)+
#   scale_fill_manual(values = c("Background" = "grey",
#                                "PCA" = "#0072B2",
#                                "KBA" = "#E69F00",
#                                "IPA" = "#009E73"))+
#   labs(x = "Value",
#        y = NULL,
#        fill = NULL) +
#   theme_bw()+
#   theme(legend.position = "bottom",
#         panel.background = element_blank(),
#         panel.grid = element_blank())
# 
# env_dat[[2]] %>%
#   mutate(across(where(is.numeric), ~ decostand(.x, method = "standardize"))) %>% 
#   pivot_longer(cols = -Area, names_to = "variable", values_to = "value") %>%
#   arrange(variable, Area) %>% 
#   mutate(value = round(value, 4)) %>% 
#   ggplot(aes(x = value, y = variable, fill = Area)) +
#   geom_boxplot(outliers = FALSE)+
#   scale_fill_manual(values = c("Background" = "grey",
#                                "PCA" = "#0072B2",
#                                "KBA" = "#E69F00",
#                                "IPA" = "#009E73"))+
#   labs(x = "Value",
#        y = NULL,
#        fill = NULL) +
#   theme_bw()+
#   theme(legend.position = "bottom",
#         panel.background = element_blank(),
#         panel.grid = element_blank())

rbind_env_dat <- do.call(rbind, env_dat) %>%
  filter(Area != "Background") %>% 
  mutate(across(where(is.numeric), ~ decostand(.x, method = "standardize"))) %>% 
  pivot_longer(cols = -c(Area, Country), names_to = "variable", values_to = "value") %>% 
  arrange(variable, Area) %>% 
  mutate(value = round(value, 4)) %>% 
  group_by(Area, variable, Country) %>% 
  summarise(mean = mean(value),
            se = sd(value)/sqrt(n()))

x <- do.call(rbind, env_dat) %>%
  filter(Area != "Background") %>% 
  mutate(across(where(is.numeric), ~ decostand(.x, method = "standardize")))

sd(x$Bio8)

country_names <- unique(rbind_env_dat$Country)

p1 <- subset(rbind_env_dat, Country %in% country_names[1:9]) %>% 
  ggplot(aes(x = mean, y = variable, fill = Area)) +
  geom_col(position = "dodge", width = 0.7)+
  geom_errorbar(aes(xmin = mean-se, xmax = mean+se),
                position = position_dodge(width = 0.9), width = 0.3)+
  scale_fill_manual(values = c("Background" = "grey",
                               "PCA" = "#0072B2",
                               "KBA" = "#E69F00",
                               "IPA" = "#009E73"))+
  labs(x = "Mean Standardised Value",
       y = NULL,
       fill = NULL) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~Country, ncol = 3)

ggsave("OneDrive - University of Cambridge/External_Projects/Kew/EnvDifferentiation_Facet1.png", p1, 
       width = 180, height = 260, units = "mm")

p2 <- subset(rbind_env_dat, Country %in% country_names[10:18]) %>% 
  ggplot(aes(x = mean, y = variable, fill = Area)) +
  geom_col(position = "dodge", width = 0.7)+
  geom_errorbar(aes(xmin = mean-se, xmax = mean+se),
                position = position_dodge(width = 0.9), width = 0.3)+
  scale_fill_manual(values = c("Background" = "grey",
                               "PCA" = "#0072B2",
                               "KBA" = "#E69F00",
                               "IPA" = "#009E73"))+
  labs(x = "Mean Standardised Value",
       y = NULL,
       fill = NULL) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~Country, ncol = 3)

ggsave("OneDrive - University of Cambridge/External_Projects/Kew/EnvDifferentiation_Facet2.png", p2, 
       width = 180, height = 260, units = "mm")

p3 <- subset(rbind_env_dat, Country %in% country_names[19:length(country_names)]) %>% 
  ggplot(aes(x = mean, y = variable, fill = Area)) +
  geom_col(position = "dodge", width = 0.7)+
  geom_errorbar(aes(xmin = mean-se, xmax = mean+se),
                position = position_dodge(width = 0.9), width = 0.3)+
  scale_fill_manual(values = c("Background" = "grey",
                               "PCA" = "#0072B2",
                               "KBA" = "#E69F00",
                               "IPA" = "#009E73"))+
  labs(x = "Mean Standardised Value",
       y = NULL,
       fill = NULL) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(~Country, ncol = 3)

ggsave("OneDrive - University of Cambridge/External_Projects/Kew/EnvDifferentiation_Facet3.png", p3, 
       width = 180, height = 260, units = "mm")

##################################### ---- Stats and Tests ---- ####################################
# Combine the data from all the countries into one global dataframe
env_dat_all <- do.call(rbind, env_dat) %>% 
  filter(Area != "Background")

total_samples <- 10000

sampled_env_dat_all <- env_dat_all %>% 
  group_by(Country, Area) %>% 
  mutate(group_size = n()) %>% 
  sample_frac(total_samples/nrow(.), weight = group_size) %>% 
  ungroup()

table(sampled_env_dat_all$Country, sampled_env_dat_all$Area)

vars <- colnames(sampled_env_dat_all)[-c(1,15,16)]

res <- lapply(vars, function(x) {
  out <- kruskal.test(as.formula(paste(x, "~ Area")), data = sampled_env_dat_all)
  out <- data.frame(variable = x, stat = out$statistic, p = out$p.value)
})
res <- do.call(rbind, res)
rownames(res) <- NULL

pair <- lapply(vars, function(x) {
  rstatix::dunn_test(formula = as.formula(paste(x, "~ Area")), 
                     data = sampled_env_dat_all, 
                     p.adjust.method = "bonferroni",
                     detailed = FALSE) %>% 
    mutate(group = paste0(group1, "/", group2)) %>% 
    dplyr::select(variable = .y., group, statistic, p.adj) %>% 
    pivot_wider(names_from = "group", values_from = c("statistic", "p.adj"))
})
pair <- do.call(rbind, pair)


?DunnTest
diff <- left_join(res, pair, by = "variable")

write.csv(diff, file = "OneDrive - University of Cambridge/External_Projects/Kew/EnvVarDifferentiation.csv",
          row.names = FALSE)

# Summary data
summary <- env_dat_all %>%
  group_by(Area) %>%
  summarise(across(
    where(is.numeric), 
    list(
      median = median
    ), 
    na.rm = TRUE
  ))

summary <- t(summary)
colnames(summary) <- summary[1, ]
summary <- summary[-1, ] %>% 
  as.data.frame() 

summary$IPA <- as.numeric(summary$IPA)
summary$KBA <- as.numeric(summary$KBA)
summary$PCA <- as.numeric(summary$PCA)

summary <- summary %>% 
  mutate(across(where(is.numeric), ~ round(., 2))) %>% 
  dplyr::select(PCA, KBA, IPA)

write.csv(summary, "OneDrive - University of Cambridge/External_Projects/Kew/env_summary_dat.csv",
          row.names = FALSE)


summary_mean <- env_dat_all %>%
  group_by(Area) %>%
  summarise(across(
    where(is.numeric), 
    list(
      mean = mean
    ), 
    na.rm = TRUE
  ))

summary_mean <- t(summary_mean)
colnames(summary_mean) <- summary_mean[1, ]
summary_mean <- summary_mean[-1, ] %>% 
  as.data.frame() 

