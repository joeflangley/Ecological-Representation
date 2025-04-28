# Load libraries
library(sf)
library(tidyverse)
library(terra)
library(raster)
library(readxl)
library(writexl)
library(lubridate)
library(vegan)
library(adehabitatHR)
library(units)

############################# ---- Data Preparation ---- ################################# 
# Read in the country names vector and then the country boundaries
countries <- st_read("data/raw_data/country_boundaries/countries.shp") %>% 
  filter(NAME_ENGL != "China" & NAME_ENGL != "India")

# Read in the prepared PCA, KBA and IPA data, add a year to the missing KBA rows and
# remove India and China's networks for this analysis
pa_data <- st_read("data/output_data/prepared_areas/wdpa_wdoecm_prepared.shp") %>% 
  filter(NAME_ENGL != "China" & NAME_ENGL != "India")

ipa_data <- st_read("data/output_data/prepared_areas/important_plant_areas_prepared.shp") %>% 
  filter(NAME_EN != "China" & NAME_EN != "India") %>% 
  drop_na(Dt_f_p_) %>% # drop one row that says it's in Nepal but it's an India PA
  dplyr::select(NAME_EN, IUCN_nm, Dt_f_p_) %>% 
  rename(NAME_ENGL = NAME_EN, IUCN_name = IUCN_nm, STATUS_YR = Dt_f_p_)

kba_data <- st_read("data/output_data/prepared_areas/key_biodiversity_areas_prepared.shp") %>% 
  filter(NAME_ENGL != "China" & NAME_ENGL != "India") %>% 
  mutate(STATUS_YR = year(AddedDate)) %>% 
  group_by(NAME_ENGL) %>%
  mutate(STATUS_YR = ifelse(STATUS_YR == -1, sample(STATUS_YR[STATUS_YR != -1], size = 1), STATUS_YR)) %>% # randomly assign year from a KBA with the same NAME_ENGL if it doesn't have a year
  ungroup() %>% 
  dplyr::select(NAME_ENGL, IUCN_name, STATUS_YR)

any(is.na(pa_data$STATUS_YR))
any(is.na(ipa_data$STATUS_YR))
any(is.na(kba_data$STATUS_YR))

# Break the PCA, KBA, IBA and IPA data down into lists, where each element is a country
pa_list <- split(pa_data, pa_data$IUCN_name)
kba_list <- split(kba_data, kba_data$IUCN_name)
ipa_list <- split(ipa_data, ipa_data$IUCN_name)

# Read the standardised environmental extractions back in
standardised_extractions <- lapply(excel_sheets("data/output_data/standardised_env_extractions.xlsx"), 
                                   function(sheet) {
                                     read_excel("data/output_data/standardised_env_extractions.xlsx", sheet = sheet)
                                   })

# # Remove China and India from the standardised extractions list
# standardised_extractions <- standardised_extractions[-c(6, 12)]

# Function to perform intersection check for all points with each polygon dataset
check_intersections <- function(points_sf, pa_polygons, kba_polygons, ipa_polygons) {
  # Find intersections using st_intersects (returns logical matrices)
  pa_intersections <- st_join(points_sf, pa_polygons, join = st_intersects)
  kba_intersections <- st_join(points_sf, kba_polygons, join = st_intersects)
  ipa_intersections <- st_join(points_sf, ipa_polygons, join = st_intersects)
  
  # Convert logical matrices to dfs and add area names
  pa_intersect_points <- pa_intersections %>% 
    filter(!is.na(STATUS_YR)) %>% 
    group_by(ID) %>% 
    slice_min(STATUS_YR) %>%         
    ungroup() %>% 
    dplyr::select(1:15, STATUS_YR) %>% 
    distinct() %>% 
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "PCA") # Add area name
  
  kba_intersect_points <- kba_intersections %>% 
    filter(!is.na(STATUS_YR)) %>% 
    group_by(ID) %>% 
    slice_min(STATUS_YR) %>%         
    ungroup() %>% 
    dplyr::select(1:15, STATUS_YR) %>% 
    distinct() %>% 
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "KBA") # Add area name
  
  ipa_intersect_points <- ipa_intersections %>% 
    filter(!is.na(STATUS_YR)) %>% 
    group_by(ID) %>% 
    slice_min(STATUS_YR) %>%         
    ungroup() %>% 
    dplyr::select(1:15, STATUS_YR) %>% 
    distinct() %>% 
    st_drop_geometry() %>%   # Drop geometry
    mutate(Area = "IPA") # Add area name
  
  # Combine all intersecting points into a single dataframe
  combined_df <- rbind(pa_intersect_points, kba_intersect_points, ipa_intersect_points)
  
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
  ipa_dat <- ipa_list[[i]]
  
  # Covert the environmental extractions df to an sf object
  points_sf <- st_as_sf(points_data, coords = c("x", "y"), crs = "+proj=moll") 
  
  # Use the function above to find sampled points inside PAs, KBAs, IBAs and IPAs
  intersect_results <- check_intersections(points_sf, pa_dat, kba_dat, ipa_dat)
  
  # Prepare the background environmental extractions df to match the format of intersect_results
  background_df <- points_data %>% 
    mutate(STATUS_YR = "Background",
           Area = "Background") %>%
    dplyr::select(-x, -y)
  
  # Combine the background_df and intersect_results
  intersect_results <- rbind(background_df, intersect_results)
  
  env_intersections[[i]] <- intersect_results
}

# Save the standardised environmental intersections in excel
sheet_names <- countries$IUCN_name
named_env_intersections <- setNames(env_intersections, sheet_names)

write_xlsx(named_env_intersections, path = "data/output_data/standardised_env_intersections_withDates.xlsx",
           col_names = TRUE, format = TRUE)

# Read the standardised environmental intersections back in
env_intersections <- lapply(excel_sheets("data/output_data/standardised_env_intersections_withDates.xlsx"), 
                            function(sheet) {
                              x <- read_excel("data/output_data/standardised_env_intersections_withDates.xlsx", sheet = sheet)
                            })

######################## ---- Principal Component Analysis ---- ########################## 
time_periods <- c(seq(1960, 2020, by = 10), 2024)

# Conduct principal component analysis (PCA) on each country's environmental data
pca_list <- list()
for (i in 1:nrow(countries)) {
  # Conduct PCA on the background points
  background_pca <- subset(env_intersections[[i]], Area == "Background")[,-c(14:16)] %>% 
    prcomp()
  
  # Get summary and model data from the background_pca
  pca_model <- background_pca
  pca_summary <- summary(background_pca)
  
  # Filter the standardised data to exclude the background points
  others_data <- subset(env_intersections[[i]], Area != "Background")[,-c(14:16)]
  
  # Conduct PCA on the PCA, KBA, IBA and IPA points
  others_pca <- data.frame(predict(background_pca, newdata = others_data))
  
  # Add the Area and STATUS_YR columns to the others_pca
  others_pca <- cbind(
    STATUS_YR = subset(env_intersections[[i]], Area != "Background")[,15],
    Area = subset(env_intersections[[i]], Area != "Background")[,16],
    others_pca
  )
  
  pre1960 <- others_pca %>% filter(STATUS_YR <= 1960) %>% 
    mutate(STATUS_YR = 1960)
  pre1970 <- others_pca %>% filter(STATUS_YR <= 1970) %>% 
    mutate(STATUS_YR = 1970)
  pre1980 <- others_pca %>% filter(STATUS_YR <= 1980) %>% 
    mutate(STATUS_YR = 1980)
  pre1990 <- others_pca %>% filter(STATUS_YR <= 1990) %>% 
    mutate(STATUS_YR = 1990)
  pre2000 <- others_pca %>% filter(STATUS_YR <= 2000) %>% 
    mutate(STATUS_YR = 2000)
  pre2010 <- others_pca %>% filter(STATUS_YR <= 2010) %>% 
    mutate(STATUS_YR = 2010)
  pre2020 <- others_pca %>% filter(STATUS_YR <= 2020) %>% 
    mutate(STATUS_YR = 2020)
  pre2024 <- others_pca %>% filter(STATUS_YR <= 2024) %>% 
    mutate(STATUS_YR = 2024)
  
  # Add the Area and STATUS_YR columns to the background_pca
  background_pca <- cbind(
    STATUS_YR = subset(env_intersections[[i]], Area == "Background")[,15],
    Area = subset(env_intersections[[i]], Area == "Background")[,16],
    data.frame(background_pca$x)
  )
  
  # Combine the background and others PCA dataframes into one
  dat <- rbind(background_pca, pre1960, pre1970, pre1980, pre1990, pre2000, pre2010, pre2020, pre2024)
  
  pca_list[[i]] <- list(PCA_df = dat, PCA_model = pca_model, PCA_summary = pca_summary)
}

########################### ---- Minimum Convex Polygons ---- ############################ 
time_periods <- c(seq(1960, 2020, by = 10), 2024)

# Modified function to filter based on time period (STATUS_YR)
calculate_minimum_convex_polygons <- function(pca_data, area, pct, time_period) {
  # Filter the data for the current time period and area
  pca_data_filtered <- subset(pca_data, STATUS_YR == time_period & Area == area)
  
  # Check if we have enough data
  if (nrow(pca_data_filtered) >= 5) {
    # Create the MCP for the filtered data
    mcp_result <- pca_data_filtered %>% 
      dplyr::select(PC1, PC2) %>% 
      SpatialPointsDataFrame(coords = ., 
                             data = ., 
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% 
      mcp(., percent = pct)
    
    return(mcp_result)
  } else {
    return(NULL)  # Return NULL if not enough data for MCP
  }
}

# List to store results for each time period
mcp_list <- list()

# Loop over each country
for (i in 1:nrow(countries)) {
  country_mcp_list <- list()
  
  # Calculate the background MCP for the specified time period (2024)
  background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                      area = "Background", pct = 100, 
                                                      time_period = "Background")
  
  if (!is.null(background_mcp)) {
    background_data <- as.data.frame(background_mcp@polygons[[1]]@Polygons[[1]]@coords) %>%
      mutate(Area = "Background", STATUS_YR = "Background") %>%   # Label it as "Background"
      relocate(STATUS_YR, Area)
    country_mcp_list[["Background"]] <- background_data
  }
  
  # Loop over each time period
  for (time_period in time_periods) {
    # Create a list to hold MCP results for the given time period
    mcp_data_list <- list()
    
    # PCA MCP
    pca_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "PCA", pct = 99, 
                                                 time_period = time_period)
    if (!is.null(pca_mcp)) {
      mcp_data_list <- append(mcp_data_list, list(as.data.frame(pca_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                                                    mutate(Area = "PCA", STATUS_YR = time_period) %>% 
                                                    relocate(STATUS_YR, Area)))
    }
    
    # KBA MCP
    kba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "KBA", pct = 99, 
                                                 time_period = time_period)
    if (!is.null(kba_mcp)) {
      mcp_data_list <- append(mcp_data_list, list(as.data.frame(kba_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                                                    mutate(Area = "KBA", STATUS_YR = time_period) %>% 
                                                    relocate(STATUS_YR, Area)))
    }
    
    # IPA MCP
    ipa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "IPA", pct = 99, 
                                                 time_period = time_period)
    if (!is.null(ipa_mcp)) {
      mcp_data_list <- append(mcp_data_list, list(as.data.frame(ipa_mcp@polygons[[1]]@Polygons[[1]]@coords) %>% 
                                                    mutate(Area = "IPA", STATUS_YR = time_period) %>% 
                                                    relocate(STATUS_YR, Area)))
    }
    
    # Combine MCP results for the current time period
    if (length(mcp_data_list) > 0) {
      dat <- bind_rows(mcp_data_list)
      country_mcp_list[[as.character(time_period)]] <- dat
    }
  }
  
  # Store the MCPs for this country
  mcp_list[[i]] <- country_mcp_list
}

for (i in 1:length(mcp_list)) {
  mcp_list[[i]] <- do.call(rbind, mcp_list[[i]])
  rownames(mcp_list[[i]]) <- NULL
}

################################## ---- PCA Plots ---- ################################### 
# Combine the principal component analysis and CHULL data into one dataframe per country
ggplot_friendly_list <- list()
for (i in 1:nrow(countries)) {
  # Merge the PCA and CHULL data
  dat <- rbind(pca_list[[i]]$PCA_df[,1:4] %>% # retain only the Area, Year, PC1 and PC2 columns
                 mutate(Type = "point"), # assign Type value of "point" to denote PCA points
               mcp_list[[i]] %>% 
                 filter(Area %in% c("Background", "PCA")) %>% 
                 mutate(Type = "polygon") # assign "polygon" to denote MCP
  )
  
  # Add a country column
  dat <- dat %>% 
    mutate(Country = unique(countries$IUCN_name)[i])
  
  ggplot_friendly_list[[i]] <- dat
}

names(ggplot_friendly_list) <- countries$IUCN_name

countries_facet_a <- do.call(rbind, ggplot_friendly_list[1:15])
p1 <- ggplot() +
  geom_point(data = subset(countries_facet_a, Type == "point"),
             mapping = aes(x = PC1, y = PC2, 
                           color = as.factor(reorder(STATUS_YR, STATUS_YR == "Background", decreasing = TRUE))),
             alpha = 0.1, size = 0.8, shape = 1)+
  geom_polygon(data = subset(countries_facet_a, Type == "polygon"),
               mapping = aes(x = PC1, y = PC2, 
                             color = as.factor(reorder(STATUS_YR, STATUS_YR == "Background", decreasing = TRUE))),
               alpha = 0, linewidth = 0.7) +
  scale_color_manual(values = c("Background" = "grey80",
                                "2024" = "#440154",
                                "2020" = "#46327e",
                                "2010" = "#365c8d",
                                "2000" = "#277f8e",
                                "1990" = "#1fa187",
                                "1980" = "#4ac16d",
                                "1970" = "#a0da39",
                                "1960" = "#fde725"))+
  labs(x = "PC1", y = "PC2", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~ Country, scales = "free", ncol = 3)

ggsave(filename = "outputs/Environmental_Expansion/Supplementary_Materials/PCA_Plots/supp_pca_a.png",
       plot = p1,
       width = 180, height = 297, units = "mm",
       dpi = 600)

countries_facet_b <- do.call(rbind, ggplot_friendly_list[16:26])
p2 <- ggplot() +
  geom_point(data = subset(countries_facet_b, Type == "point"),
             mapping = aes(x = PC1, y = PC2, 
                           color = as.factor(reorder(STATUS_YR, STATUS_YR == "Background", decreasing = TRUE))),
             alpha = 0.1, size = 0.8, shape = 1)+
  geom_polygon(data = subset(countries_facet_b, Type == "polygon"),
               mapping = aes(x = PC1, y = PC2, 
                             color = as.factor(reorder(STATUS_YR, STATUS_YR == "Background", decreasing = TRUE))),
               alpha = 0, linewidth = 0.7) +
  scale_color_manual(values = c("Background" = "grey80",
                                "2024" = "#440154",
                                "2020" = "#46327e",
                                "2010" = "#365c8d",
                                "2000" = "#277f8e",
                                "1990" = "#1fa187",
                                "1980" = "#4ac16d",
                                "1970" = "#a0da39",
                                "1960" = "#fde725"))+
  labs(x = "PC1", y = "PC2", color = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~ Country, scales = "free", ncol = 3)

ggsave(filename = "outputs/Environmental_Expansion/Supplementary_Materials/PCA_Plots/supp_pca_b.png",
       plot = p2,
       width = 180, height = 240, units = "mm",
       dpi = 600)

################################## ---- MCP Area ---- #################################### 
# Modified function to filter based on time period (STATUS_YR)
calculate_minimum_convex_polygons <- function(pca_data, area, pct, time_period) {
  # Filter the data for the current time period and area
  pca_data_filtered <- subset(pca_data, STATUS_YR == time_period & Area == area)
  
  # Check if we have enough data
  if (nrow(pca_data_filtered) >= 5) {
    # Create the MCP for the filtered data
    mcp_result <- pca_data_filtered %>% 
      dplyr::select(PC1, PC2) %>% 
      SpatialPointsDataFrame(coords = ., 
                             data = ., 
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) %>% 
      mcp(., percent = pct) %>% 
      st_as_sf()
    
    return(mcp_result)
  } else {
    return(NULL)  # Return NULL if not enough data for MCP
  }
}

# Helper function to safely get area or assign 0 if empty, i.e. prevent error if area = 0
safe_area <- function(x) {
  if (length(x) == 0) {
    return(0)
  } else {
    return(as.numeric(st_area(x)))
  }
}
  
# List to store results for each time period
mcp_area <- list()

# Loop over each country
for (i in 1:nrow(countries)) {
  country_mcp_list <- list()
  
  # Calculate the background MCP for the specified time period
  background_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                      area = "Background", pct = 100, 
                                                      time_period = "Background")

  # Loop over each time period
  for (time_period in time_periods) {
    # Create a list to hold MCP results for the current time period
    mcp_data_list <- list()
    
    # PCA MCP
    pca_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "PCA", pct = 99, 
                                                 time_period = time_period)
    
    # KBA MCP
    kba_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "KBA", pct = 99, 
                                                 time_period = time_period)

    # IPA MCP
    ipa_mcp <- calculate_minimum_convex_polygons(pca_data = pca_list[[i]]$PCA_df, 
                                                 area = "IPA", pct = 99, 
                                                 time_period = time_period)

    if (length(kba_mcp) != 0 & length(ipa_mcp) != 0){
      kba_ipa_joint <- st_union(kba_mcp, ipa_mcp)
    } else if (length(kba_mcp) != 0) {
      kba_ipa_joint <- kba_mcp
    } else if (length(ipa_mcp) != 0) {
      kba_ipa_joint <- ipa_mcp
    } else {
      kba_ipa_joint <- NULL
    }
    
    if (length(kba_ipa_joint) != 0 & length(pca_mcp) != 0){
      kba_ipa_additional <- st_difference(kba_ipa_joint, pca_mcp)
    } else if (length(kba_ipa_joint) != 0) {
      kba_ipa_additional <- kba_ipa_joint
    } else {
      kba_ipa_additional <- NULL
    }
      
    if (length(kba_ipa_joint) != 0 & length(pca_mcp) != 0){
      kba_ipa_pca <- st_intersection(kba_ipa_joint, pca_mcp)
    } else {
      kba_ipa_pca <- NULL
    }

    if (length(kba_ipa_joint) != 0 & length(pca_mcp) != 0){
      pca_alone <- st_difference(pca_mcp, kba_ipa_joint)
    } else if (length(pca_mcp) != 0) {
      pca_alone <- pca_mcp
    } else {
      pca_alone <- NULL
    }
    
    if (NROW(pca_alone) == 0){
      pca_alone <- NULL
    }
  
    dat <- data.frame(country = countries$IUCN_name[i],
                      designation = c("Background", "PCA", "KBA", "IPA", "KBA_IPA_joint",
                                      "KBA_IPA_additional","KBA_IPA_PCA", "PCA_alone"),
                      time_period = time_period,
                      area = rbind(safe_area(background_mcp), 
                                   safe_area(pca_mcp), 
                                   safe_area(kba_mcp), 
                                   safe_area(ipa_mcp),
                                   safe_area(kba_ipa_joint),
                                   safe_area(kba_ipa_additional),
                                   safe_area(kba_ipa_pca),
                                   safe_area(pca_alone)
                                   )
    ) %>% 
      mutate(
        background_area = area[designation == "Background"],
        area_pct = (area / background_area) * 100
      ) %>%
      dplyr::select(-area, -background_area)
    
    

    country_mcp_list[[as.character(time_period)]] <- dat

  
  # Store the MCPs for this country
  mcp_area[[i]] <- country_mcp_list
  }
}

mcp_area_df <- list()
for (i in 1:length(mcp_area)) {
  dat <- do.call(rbind, mcp_area[[i]])
  rownames(dat) <- NULL
  mcp_area_df[[i]] <- dat
} 

mcp_area_df <- do.call(rbind, mcp_area_df)

pca_expansion_df <- mcp_area_df %>% 
  filter(designation == "PCA") %>% 
  dplyr::select(-designation) %>% 
  pivot_wider(names_from = time_period, values_from = area_pct)
write.csv(pca_expansion_df, file = "data/output_data/env_space_through_time.csv", row.names = FALSE)

##################################### ---- Area Plots ---- ###############################
area_plot_df <- mcp_area_df %>% 
  filter(designation %in% c("KBA_IPA_additional", "KBA_IPA_PCA", "PCA_alone"))

country_names <- unique(countries$IUCN_name)

p1 <- area_plot_df %>% 
  filter(country %in% country_names[1:15]) %>% 
  ggplot()+
  geom_area(aes(x = time_period, y = area_pct, group = designation, color = designation, 
                fill = designation),
            position = "stack", alpha = 0.5, linewidth = 0.7)+
  labs(x = "Year", y = "Coverage of environmental space (%)",
       fill = NULL, 
       color = NULL)+
  scale_color_manual(values = c(PCA_alone = "#0072B2",
                                KBA_IPA_PCA = "#440154",
                                KBA_IPA_additional = "#739F3A"),
                     breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                     labels = c("PCAs", "Important areas", "Both"))+
  scale_fill_manual(values = c(PCA_alone = "#0072B2",
                               KBA_IPA_PCA = "#440154",
                               KBA_IPA_additional = "#739F3A"),
                    breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                    labels = c("PCAs", "Important areas", "Both"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0, 0, 0, 5))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 0.25, vjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        panel.spacing = unit(0.5, "lines"))+
  facet_wrap(~ country, ncol = 3)

ggsave(filename = "outputs/Environmental_Expansion/Supplementary_Materials/supp_area_plot_a.png",
       plot = p1,
       width = 180, height = 297, units = "mm",
       dpi = 600)

p2 <- area_plot_df %>% 
  filter(country %in% country_names[16:26]) %>% 
  ggplot()+
  geom_area(aes(x = time_period, y = area_pct, group = designation, color = designation, 
                fill = designation),
            position = "stack", alpha = 0.5, linewidth = 0.7)+
  labs(x = "Year", y = "Coverage of environmental space (%)",
       fill = NULL, 
       color = NULL)+
  scale_color_manual(values = c(PCA_alone = "#0072B2",
                                KBA_IPA_PCA = "#440154",
                                KBA_IPA_additional = "#739F3A"),
                     breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                     labels = c("PCAs", "Important areas", "Both"))+
  scale_fill_manual(values = c(PCA_alone = "#0072B2",
                               KBA_IPA_PCA = "#440154",
                               KBA_IPA_additional = "#739F3A"),
                    breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                    labels = c("PCAs", "Important areas", "Both"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0, 0, 0, 5))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 0.25, vjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 10),
        panel.spacing = unit(0.5, "lines"))+
  facet_wrap(~ country, ncol = 3)

ggsave(filename = "outputs/Environmental_Expansion/Supplementary_Materials/supp_area_plot_b.png",
       plot = p2,
       width = 180, height = 240, units = "mm",
       dpi = 600)

mean_dat <- area_plot_df %>% 
  group_by(designation, time_period) %>% 
  summarise(mean = mean(area_pct)) %>% 
  mutate(country = "Mean")

mean_plot <- mean_dat %>% 
  ggplot()+
  geom_area(aes(x = time_period, y = mean, group = designation, color = designation, 
                fill = designation),
            position = "stack", alpha = 0.5, linewidth = 0.7)+
  labs(x = NULL, y = "Coverage of environmental space (%)",
       fill = NULL, 
       color = NULL)+
  scale_color_manual(values = c(PCA_alone = "#0072B2",
                                KBA_IPA_PCA = "#440154",
                                KBA_IPA_additional = "#739F3A"),
                     breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                     labels = c("PCAs", "Important areas", "Both"))+
  scale_fill_manual(values = c(PCA_alone = "#0072B2",
                               KBA_IPA_PCA = "#440154",
                               KBA_IPA_additional = "#739F3A"),
                    breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                    labels = c("PCAs", "Important areas", "Both"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0, 0, 0, 5))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~ country, ncol = 1)

rep_plot <- area_plot_df %>% 
  filter(country %in% c("Algeria", "Libya", "Türkiye", "United Kingdom")) %>% 
  ggplot()+
  geom_area(aes(x = time_period, y = area_pct, group = designation, color = designation, 
                fill = designation),
            position = "stack", alpha = 0.5, linewidth = 0.7)+
  labs(x = "Year", y = "Coverage of environmental space (%)",
       fill = NULL, 
       color = NULL)+
  scale_color_manual(values = c(PCA_alone = "#0072B2",
                                KBA_IPA_PCA = "#440154",
                                KBA_IPA_additional = "#739F3A"),
                     breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                     labels = c("PCAs", "Important areas", "Both"))+
  scale_fill_manual(values = c(PCA_alone = "#0072B2",
                               KBA_IPA_PCA = "#440154",
                               KBA_IPA_additional = "#739F3A"),
                    breaks = c("PCA_alone", "KBA_IPA_additional", "KBA_IPA_PCA"),
                    labels = c("PCAs", "Important areas", "Both"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0, 0, 0, 5))+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.text = element_text(size = 10))+
  facet_wrap(~ country, ncol = 2)

library(ggpubr)
fig <- ggarrange(mean_plot, rep_plot, 
          ncol = 1,
          labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom",
          hjust = -1,
          heights = c(1, 1.25))


ggsave(plot = fig, "outputs/Environmental_Expansion/fig.png", 
       width = 180, height = 240, units = "mm", dpi = 300,
       device = "png")
