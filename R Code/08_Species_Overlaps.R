# Load Libraries
library(sf)
library(tidyverse)
library(raster)
library(terra)
library(data.table)
library(pryr)
library(rstatix)

# Read in boundaries of each of the studied countries
countries <- st_read("data/raw_data/country_boundaries/countries.shp")

# Generate a 5000m Molleweide spatial raster template
r <- raster()
x <- projectRaster(r, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
r2 <- x
res(r2) <- 5000
raster_template <- resample(x, r2)
spat_rast_template <- rast(raster_template) # turn from raster to spatial raster

# Read in and rasterize the PCA, KBA and IPA polygons
pa_data <- st_read("data/output_data/prepared_areas/wdpa_wdoecm_prepared.shp")
osm_pa_data <- st_read("data/output_data/prepared_areas/IND_CHN_osm_prepared.shp") %>% 
  rename("NAME" = "name")
pa_data <- bind_rows(pa_data, osm_pa_data) # combine the PCA data into one dataframe
ipa_data <- st_read("data/output_data/prepared_areas/important_plant_areas_prepared.shp") %>% 
  rename(IUCN_name = IUCN_nm)
kba_data <- st_read("data/output_data/prepared_areas/key_biodiversity_areas_prepared.shp")

# Function to rasterize polygons without missing small polygons
rasterize_process <- function(input_data, template_data) {
  vv <- vect(input_data)
  cnt <- centroids(vv)
  x <- terra::rasterize(vv, template_data, touches = T)
  y <- terra::rasterize(cnt, template_data, touches = T)
  cover(x, y)
}

# Apply for each country boundary
country_rasters <- lapply(1:nrow(countries), function(i){
  rasterize_process(countries[i, ], spat_rast_template)
})

# Apply for the PCA, KBA & IPA polygons
pca_rasters <- rasterize_process(pa_data, spat_rast_template)
kba_rasters <- rasterize_process(kba_data, spat_rast_template)
ipa_rasters <- rasterize_process(ipa_data, spat_rast_template)

rm(pa_data, kba_data, osm_pa_data, ipa_data, r, x, r2)

# Read in the species spatial range data (very large so possibly better to run in chunks)
species_files <- list.files("data/output_data/prepared_species_ranges/",
                            pattern = ".gpkg",
                            full.names = TRUE)
species_sf <- lapply(species_files, function(x) {
  st_read(x)
})

# Combine the species data into one dataframe
species_sf <- do.call(rbind, species_sf)

# Set the output directory
output_dir <- "data/output_data/prepared_species_ranges/country_species_results"

# Set the cell area - this doesn't matter as long as it stays consistent
cell_area_km2 <- 1

# For each country & species, find the species' overlap area with each country and its PCAs/KBAs/IPAs
for (i in 1:nrow(countries)) {
  country <- countries[i, ]
  
  # Crop the country/PCA/KBA/IPA rasters by the country's boundary
  country_raster <- crop(country_rasters[[i]], country)
  pca_raster_national <- crop(pca_rasters, country)
  kba_raster_national <- crop(kba_rasters, country)
  ipa_raster_national <- crop(ipa_rasters, country)
  
  # Compute country-specific masks for PCAs, KBAs, and IPAs
  country_pca_mask <- country_raster * pca_raster_national
  country_kba_mask <- country_raster * kba_raster_national
  country_ipa_mask <- country_raster * ipa_raster_national
  
  # Loop over each species
  for (j in 1:nrow(species_sf)) {
    
    # Rasterize each species using the function that does not miss small polygons
    species_raster <- rasterize_process(species_sf[j, ], spat_rast_template)
    
    # Crop the species raster by the country
    species_raster <- crop(species_raster, country)
    species_raster_national <- species_raster * country_raster
    
    # Calculate the national area of the species' range
    national_area <- sum(!is.na(values(species_raster_national))) * cell_area_km2
    
    # Compute the species overlap area with PCAs/KBAs/IPAs, only if it has a national area
    if (national_area > 0) {
      pca_area <- sum(!is.na(values(species_raster_national * country_pca_mask))) * cell_area_km2
      kba_area <- sum(!is.na(values(species_raster_national * country_kba_mask))) * cell_area_km2
      ipa_area <- sum(!is.na(values(species_raster_national * country_ipa_mask))) * cell_area_km2
      
      # Prepare results for the current species
      results <- data.frame(
        country = country$NAME_ENGL,
        species = species_sf[j,]$SCI_NAME,
        national_area = national_area,
        pca_area = pca_area,
        kba_area = kba_area,
        ipa_area = ipa_area,
        stringsAsFactors = FALSE
      )
      
      # File path for the current country's CSV         
      country_csv <- file.path(
        output_dir,
        paste0(country$NAME_ENGL, "_species_results.csv") 
      )
      
      # Write header row if the CSV does not exist
      if (!file.exists(country_csv)) {
        write.csv(
          data.frame(
            country = character(),
            species = character(),
            national_area = numeric(),
            pca_area = numeric(),
            kba_area = numeric(),
            ipa_area = numeric(),
            stringsAsFactors = FALSE
          ),
          file = country_csv,
          row.names = FALSE
        )
      }
      
      # Append results to the CSV
      fwrite(results, file = country_csv, append = TRUE)
      
      # Use garbage collection to reduce memory usage every 200 species
      if (j %% 200 == 0) { 
        gc()
      }
    }
  }
  
  # Use garbage collection after each country iteration
  gc()
}

######################## ---- Calculate Global Species Range ---- ########################
# Vector of file names for each prepared species range geopackage
species_files <- list.files("data/output_data/prepared_species_ranges/",
                            pattern = ".gpkg",
                            full.names = TRUE)

# Set the output directory
output_file <- "data/output_data/prepared_species_ranges/country_species_results/species_global_area.csv"

# Write header row if the CSV does not exist
if (!file.exists(output_file)) {
  write.csv(
    data.frame(
      species = character(),
      global_area = numeric(),
      stringsAsFactors = FALSE
    ),
    file = output_file,
    row.names = FALSE
  )
}

for (i in 1:length(species_files)){
  species_sf <- st_read(species_files[i])
  
  for (j in 1:nrow(species_sf)) {
    # Rasterize each species using the function that does not miss small polygons
    species_raster <- rasterize_process(species_sf[j, ], spat_rast_template)
    
    # Calculate the species' global area
    global_species_area <- sum(!is.na(values(species_raster))) * cell_area_km2
    
    # Prepare results for the current species
    results <- data.frame(
      species = species_sf$SCI_NAME[j],
      global_area = as.numeric(global_species_area),
      stringsAsFactors = FALSE
    )
    
    # Append results to the CSV
    fwrite(results, file = output_file, append = TRUE)
    
    # Use garbage collection to reduce memory usage every 200 species
    if (j %% 200 == 0) { 
      gc()
    }
  }
  
  # Use garbage collection after each prepared species range file iteration
  gc()
}

########################## ---- Calculate IPA/KBA/PCA Area ---- ##########################
network_area <- list()
for (i in 1:nrow(countries)){
  country <- countries[i, ]
  
  # Crop the country/PCA/KBA/IPA rasters by the country's boundary
  country_raster <- crop(country_rasters[[i]], country)
  pca_raster_national <- crop(pca_rasters, country)
  kba_raster_national <- crop(kba_rasters, country)
  ipa_raster_national <- crop(ipa_rasters, country)
  
  # Compute country-specific masks for PCAs, KBAs, and IPAs
  country_pca_mask <- country_raster * pca_raster_national
  country_kba_mask <- country_raster * kba_raster_national
  country_ipa_mask <- country_raster * ipa_raster_national
  
  # Calculate the raster area values
  pca_area <- sum(!is.na(values(country_pca_mask))) * cell_area_km2
  kba_area <- sum(!is.na(values(country_kba_mask))) * cell_area_km2
  ipa_area <- sum(!is.na(values(country_ipa_mask))) * cell_area_km2
  
  # Calculate the actual area values in km2 for reference
  country_pca_mask <- st_intersection(pa_data, country)
  country_kba_mask <- st_intersection(kba_data, country)
  country_ipa_mask <- st_intersection(ipa_data, country)
  
  pca_area_sf <- st_area(st_union(country_pca_mask)) %>% units::set_units(km2)
  kba_area_sf <- st_area(st_union(country_kba_mask)) %>% units::set_units(km2)
  ipa_area_sf <- st_area(st_union(country_ipa_mask)) %>% units::set_units(km2)
  
  # Calculate the country area in km2
  country_area_sf <- st_area(st_union(country)) %>% units::set_units(km2)
  
  # Add all these area values to a dataframe
  dat <- data.frame(country = country$NAME_ENGL, pca_network_area = pca_area,
                    kba_network_area = kba_area, ipa_network_area = ipa_area,
                    pca_network_area_sf = pca_area_sf, kba_network_area_sf = kba_area_sf, 
                    ipa_network_area_sf = ipa_area_sf, country_area_sf = country_area_sf)
  
  network_area[[i]] <- dat
}

network_area <- do.call(rbind, network_area)

write.csv(network_area, file = "data/output_data/prepared_areas/area_by_country.csv",
          row.names = FALSE)

############################### ---- Format Overlap Results ---- ################################
# Read in the data for each country on species ranges in PCAs, KBAs and IPAs
country_results <- list.files("data/output_data/prepared_species_ranges/country_species_results/",
                              pattern = "species_results",
                              full.names = TRUE)
country_results <- lapply(country_results, read.csv)
country_results <- do.call(rbind, country_results)

# Read in summary data for each species
species_files <- list.files("data/output_data/prepared_species_ranges/",
                              pattern = ".gpkg",
                              full.names = TRUE)
species_summary <- lapply(species_files, function(x) {
  st_read(x) %>% 
    st_drop_geometry()
})
species_summary <- do.call(rbind, species_summary)

# Combine the species areas and species summaries
country_results <- left_join(country_results, species_summary, by = c("species" = "SCI_NAME"))

# Add in the species global range area data
global_ranges <- list.files("data/output_data/prepared_species_ranges/country_species_results/",
                            pattern = "species_global_area", full.names = TRUE)
global_ranges <- lapply(global_ranges, read.csv)
global_ranges <- do.call(rbind, global_ranges)
global_ranges <- unique(global_ranges)

country_results <- left_join(country_results, global_ranges, by = "species")

# Add in each country's PCA, KBA and IPA area values 
network_area <- read.csv("data/output_data/prepared_areas/area_by_country.csv")
country_results <- left_join(country_results, network_area[, 1:4], by = "country")

# Save the country results 
saveRDS(country_results, file = "data/output_data/prepared_species_ranges/country_species_results/results.RDS")

############################### ---- Overlap Plots and Stats ---- ################################
# Read in the country results and clean the dataframe
country_results <- readRDS("~/OneDrive - University of Cambridge/External_Projects/Kew/species_results.RDS") %>% 
  dplyr::select(ID_NO, species, taxon, status, redlistCategory, country, species, global_area, national_area,
                pca_area, kba_area, ipa_area, pca_network_area, kba_network_area, ipa_network_area) %>% 
  mutate(
    status = case_when(
      redlistCategory %in% c("Least Concern", "Near Threatened") ~ "Not Threatened",
      redlistCategory %in% c("Critically Endangered", "Endangered", "Vulnerable") ~ "Threatened",
      TRUE ~ status
    ),
    taxon = case_when(
      taxon == "Herps" ~ "Herptiles",
      TRUE ~ taxon
    )
  )

# Find the total number of species overlapping at least one country
country_results %>% 
  dplyr::select(species, status) %>% 
  distinct() %>% 
  group_by(status) %>% 
  summarise(n = n())

# Find the number of species, by status and taxon, in each country and across all countries
species_by_country <- country_results %>% 
  group_by(country, status, taxon) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = c("status", "taxon"), values_from = "n")

species_globally <- country_results %>% 
  select(species, status, taxon) %>% 
  distinct() %>% 
  group_by(status, taxon) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pivot_wider(names_from = c("status", "taxon"), values_from = "n")

all_sp <- rbind(species_by_country, species_globally)
all_sp[29,1] <- "Global"

write.csv(all_sp, file = "OneDrive - University of Cambridge/External_Projects/Kew/species_n_by_country.csv",
          row.names = FALSE)

# Find the number of species in each combination of PCA/KBA/IPA
countries <- unique(country_results$country)
species_overlaps_by_country <- list()
for (i in seq_along(countries)) {
  dat <- country_results %>% 
    filter(country == countries[i])
  
  out <- c(country = countries[i],
           "PCA" = nrow(subset(dat, pca_area != 0)),
           "PCA/KBA" = nrow(subset(dat, pca_area != 0 & kba_area != 0)),
           "PCA/IPA" = nrow(subset(dat, pca_area != 0 & ipa_area != 0)),
           "KBA/IPA" = nrow(subset(dat, kba_area != 0 & ipa_area != 0)),
           "KBA" = nrow(subset(dat, pca_area == 0 & kba_area != 0)),
           "IPA" = nrow(subset(dat, pca_area == 0 & ipa_area != 0)),
           "Total" = nrow(subset(dat, pca_area != 0 | kba_area != 0 | ipa_area != 0))
  )
  species_overlaps_by_country[[i]] <- out
}

global_overlaps <- c(country = "Summed Overall",
                     "PCA" = nrow(subset(country_results, pca_area != 0)),
                     "PCA/KBA" = nrow(subset(country_results, pca_area != 0 & kba_area != 0)),
                     "PCA/IPA" = nrow(subset(country_results, pca_area != 0 & ipa_area != 0)),
                     "KBA/IPA" = nrow(subset(country_results, kba_area != 0 & ipa_area != 0)),
                     "KBA" = nrow(subset(country_results, pca_area == 0 & kba_area != 0)),
                     "IPA" = nrow(subset(country_results, pca_area == 0 & ipa_area != 0)),
                     "Total" = nrow(subset(country_results, pca_area != 0 | kba_area != 0 | ipa_area != 0))
)

species_overlaps_by_country <- do.call(rbind, species_overlaps_by_country) %>% 
  rbind(., global_overlaps)

rownames(species_overlaps_by_country) <- NULL
write.csv(species_overlaps_by_country, "OneDrive - University of Cambridge/External_Projects/Kew/species_overlaps_by_country.csv",
          row.names = FALSE)

# Calculate overlap percentages
overlap_results <- country_results %>% 
  # Calculate overlap percentage for each species in each country's conservation areas
  mutate(
    pca_overlap = (pca_area / national_area) * 100,
    kba_overlap = (kba_area / national_area) * 100,
    ipa_overlap = (ipa_area / national_area) * 100,
  ) %>% 
  
  # Find overlap percentage per unit of conservation area, i.e. weighted by area the country's PCA/KBA/IPA network
  mutate(
    pca_overlap_per_unit = pca_overlap / pca_network_area,
    kba_overlap_per_unit = kba_overlap / kba_network_area,
    ipa_overlap_per_unit = ipa_overlap / ipa_network_area
  ) %>% 
  
  # Find overlap percentage weighted by the proportion of the species' global area in that country
  mutate(
    pca_overlap_weighted = pca_overlap * (national_area/global_area),
    kba_overlap_weighted = kba_overlap * (national_area/global_area),
    ipa_overlap_weighted = ipa_overlap * (national_area/global_area),
  ) %>% 
  
  # Find proportion of species' global range in given country
  mutate(
    global_national_pct = (national_area/global_area)*100,
    endemic = case_when(
      global_area != national_area ~ FALSE,
      global_area == national_area ~ TRUE
    )
  )

####################################### ---- Bar Plots ---- ########################################
bar_dat <- rbind(
  overlap_results %>% 
    pivot_longer(pca_overlap:ipa_overlap, names_to = "Network", values_to = "OverlapPct") %>% 
    mutate(
      Network = case_when(
        Network == "ipa_overlap" ~ "Important Plant Areas",
        Network == "kba_overlap" ~ "Key Biodiversity Areas",
        Network == "pca_overlap" ~ "Protected & Conserved Areas",
      )) %>% 
    group_by(taxon, status, Network) %>% 
    summarise(mean = mean(OverlapPct),
              se = sd(OverlapPct) / sqrt(n())) %>% 
    mutate(OverlapVar = "Raw") %>% 
    rename(Group = status) %>% 
    dplyr::select(taxon, Group, OverlapVar, Network, mean, se),
  
  overlap_results %>% 
    filter(endemic == TRUE) %>% 
    pivot_longer(pca_overlap:ipa_overlap, names_to = "Network", values_to = "OverlapPct") %>% 
    mutate(
      Network = case_when(
        Network == "ipa_overlap" ~ "Important Plant Areas",
        Network == "kba_overlap" ~ "Key Biodiversity Areas",
        Network == "pca_overlap" ~ "Protected & Conserved Areas",
      )) %>% 
    group_by(taxon, Network) %>% 
    summarise(mean = mean(OverlapPct),
              se = sd(OverlapPct) / sqrt(n())) %>% 
    mutate(Group = "Endemic", 
           OverlapVar = "Raw") %>% 
    dplyr::select(taxon, Group, OverlapVar, Network, mean, se),
  
  overlap_results %>% 
    pivot_longer(pca_overlap_per_unit:ipa_overlap_per_unit, names_to = "Network", values_to = "OverlapPct") %>% 
    mutate(
      Network = case_when(
        Network == "ipa_overlap_per_unit" ~ "Important Plant Areas",
        Network == "kba_overlap_per_unit" ~ "Key Biodiversity Areas",
        Network == "pca_overlap_per_unit" ~ "Protected & Conserved Areas"
      )) %>% 
    group_by(taxon, status, Network) %>% 
    summarise(mean = mean(OverlapPct),
              se = sd(OverlapPct) / sqrt(n())) %>% 
    mutate(OverlapVar = "Per Unit Area") %>% 
    rename(Group = status) %>% 
    dplyr::select(taxon, Group, OverlapVar, Network, mean, se),
  
  overlap_results %>% 
    filter(endemic == TRUE) %>% 
    pivot_longer(pca_overlap_per_unit:ipa_overlap_per_unit, names_to = "Network", values_to = "OverlapPct") %>% 
    mutate(
      Network = case_when(
        Network == "ipa_overlap_per_unit" ~ "Important Plant Areas",
        Network == "kba_overlap_per_unit" ~ "Key Biodiversity Areas",
        Network == "pca_overlap_per_unit" ~ "Protected & Conserved Areas"
      )
    ) %>% 
    group_by(taxon, Network) %>% 
    summarise(mean = mean(OverlapPct),
              se = sd(OverlapPct) / sqrt(n())) %>% 
    mutate(Group = "Endemic", 
           OverlapVar = "Per Unit Area") %>% 
    dplyr::select(taxon, Group, OverlapVar, Network, mean, se)
)

# Manipulate the chunk of code below to get pairwise differences in overlap percentages
dat_subset <- bar_dat %>% 
  filter(OverlapVar == "Per Unit Area" & Group == "Endemic") %>% 
  arrange(taxon, Group)
# & Network == "Protected & Conserved Areas"

dat_subset %>% 
  mutate(diff = lag(mean, default = NA) - mean) %>% 
  filter(row_number() %% 2 == 0)

# Get the mean overlap percentages by area
mean(overlap_results$pca_overlap)
mean(overlap_results$kba_overlap)
mean(overlap_results$ipa_overlap)

# Create two plots: one with the raw overlap %s and one with the %s per unit area
p1 <- subset(bar_dat, OverlapVar == "Raw") %>% 
  mutate(Group = factor(Group, levels = c("Endemic", "Threatened", "Not Threatened"))) %>%
  ggplot(mapping = aes(x = mean, y = taxon, fill = Group))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(xmin = mean-se, xmax = mean+se),
                position = position_dodge(0.9), width = 0.3)+
  labs(x = "Proportion of species' national ranges inside PCAs, KBAs and IPAs (%)",
       y = NULL,
       fill = NULL)+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.75, "lines"))+
  guides(fill = guide_legend(reverse=T))+
  facet_wrap(~factor(Network, c("Protected & Conserved Areas", 
                                "Key Biodiversity Areas", 
                                "Important Plant Areas")))

p2 <- subset(bar_dat, OverlapVar == "Per Unit Area") %>% 
  mutate(Group = factor(Group, levels = c("Endemic", "Threatened", "Not Threatened"))) %>%
  ggplot(mapping = aes(x = mean, y = taxon, fill = Group))+
  geom_col(position = "dodge")+
  geom_errorbar(aes(xmin = mean-se, xmax = mean+se),
                position = position_dodge(0.9), width = 0.3)+
  labs(x = expression(paste("Proportion of species' national ranges inside PCAs, KBAs and IPAs, per ", 25^2, " km cell")),
       y = NULL,
       fill = NULL)+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0.75, "lines"))+
  guides(fill = guide_legend(reverse=T))+
  facet_wrap(~factor(Network, c("Protected & Conserved Areas", 
                                "Key Biodiversity Areas", 
                                "Important Plant Areas")))

# Arrange the plots into one figure
library(ggpubr)
species_plot <- ggarrange(p1, p2, 
                          nrow = 2, ncol = 1,
                          labels = c("A", "B"),
                          common.legend = TRUE, legend = "bottom")

ggsave(filename = "OneDrive - University of Cambridge/External_Projects/Kew/species_plot_correct_new.png",
       plot = species_plot,
       height = 160, width = 180, units = "mm",
       dpi = 600)

################################## ---- Additional Species ---- ####################################
# NATIONAL: species not in PCAs but in KBAs or IPAs
additional_species <- overlap_results %>% 
  filter(pca_area == 0) %>% 
  filter(kba_area != 0 | ipa_area != 0)

ipa_only <- subset(additional_species, kba_area == 0) # number of species only in IPAs (not KBAs or PCAs)
kba_only <- subset(additional_species, ipa_area == 0) # number of species only in KBAs (not IPAs or PCAs)
kba_ipa_only <- subset(additional_species, ipa_area != 0 | kba_area != 0) # number of species only in KBAs and IPAs (not PCAs)

table(ipa_only$taxon, ipa_only$redlistCategory)
table(kba_only$taxon, kba_only$redlistCategory)
table(ipa_only$taxon, ipa_only$endemic)
table(kba_only$taxon, kba_only$endemic)

table(kba_ipa_only$taxon, kba_ipa_only$redlistCategory) %>% 
  as.data.frame() %>% 
  pivot_wider(values_from = "Freq", names_from = "Var2")
table(kba_ipa_only$taxon, kba_ipa_only$endemic)

# GLOBAL: species not in PCAs but in KBAs or IPAs
additional_species_global <- overlap_results %>%
  group_by(species) %>% 
  filter(all(pca_area == 0)) %>%  # Ensure species has no PCA overlap in any country
  ungroup() %>%
  filter(kba_area != 0 | ipa_area != 0)  # Ensure species has overlap in KBAs/IPAs

additional_species_global <- additional_species_global %>% 
  distinct(species, .keep_all = TRUE)

ipa_only <- subset(additional_species_global, kba_area == 0) # number of species only in IPAs (not KBAs or PCAs)
nrow(subset(additional_species_global, ipa_area != 0))
table(ipa_only$status)
table(ipa_only$endemic)

kba_only <- subset(additional_species_global, ipa_area == 0) # number of species only in KBAs (not IPAs or PCAs)
nrow(subset(additional_species_global, kba_area != 0))
table(kba_only$status)
table(kba_only$endemic)      

kba_ipa_only <- subset(additional_species_global, ipa_area != 0 | kba_area != 0) # number of species only in KBAs and IPAs (not PCAs)

table(kba_ipa_only$status)
table(kba_ipa_only$endemic)

table(kba_ipa_only$taxon, kba_ipa_only$status)
table(kba_ipa_only$taxon, kba_ipa_only$endemic)

table(kba_ipa_only$taxon, kba_ipa_only$redlistCategory) %>% 
  as.data.frame() %>% 
  pivot_wider(values_from = "Freq", names_from = "Var2")
table(kba_ipa_only$taxon, kba_ipa_only$endemic, kba_ipa_only$status)

################################# ---- Global Species Values ---- ##################################
n_sp_pca <- overlap_results %>%
  filter(pca_area != 0) %>%
  group_by(status) %>% 
  summarise(n_species = n_distinct(species)) 
n_sp_pca
sum(n_sp_pca$n_species)

n_sp_kba <- overlap_results %>%
  filter(kba_area != 0) %>%
  group_by(status) %>% 
  summarise(n_species = n_distinct(species)) 
n_sp_kba
sum(n_sp_kba$n_species)

n_sp_ipa <- overlap_results %>%
  filter(ipa_area != 0) %>%
  group_by(status) %>% 
  summarise(n_species = n_distinct(species)) 
n_sp_ipa
sum(n_sp_ipa$n_species)

overlap_results %>%
  filter(pca_area != 0 & kba_area != 0) %>%
  summarise(n_species = n_distinct(species)) %>% 
  pull(n_species)

overlap_results %>%
  filter(pca_area != 0 & ipa_area != 0) %>%
  summarise(n_species = n_distinct(species)) %>% 
  pull(n_species)

overlap_results %>%
  filter(kba_area != 0 & ipa_area != 0) %>%
  summarise(n_species = n_distinct(species)) %>% 
  pull(n_species)

overlap_results %>%
  filter(kba_area != 0 | ipa_area != 0 | pca_area != 0) %>%
  summarise(n_species = n_distinct(species)) %>% 
  pull(n_species)

((25493-24574)/24574)*100 # % additional species in all three than PCAs alone

################################### ---- Chi Square Tests ---- #####################################
# Find total number of distinct species in each
n_sp_pca <- overlap_results %>% # In PCAs
  filter(pca_area != 0 & status == "Threatened") %>%
  group_by(taxon) %>% 
  summarise(n_sp = n_distinct(species)) %>% 
  pull(n_sp)

n_sp_kba <- overlap_results %>% # In KBAs
  filter(kba_area != 0 & status == "Threatened") %>%
  group_by(taxon) %>% 
  summarise(n_sp = n_distinct(species)) %>% 
  pull(n_sp)

n_sp_ipa <- overlap_results %>% # In IPAs
  filter(ipa_area != 0 & status == "Threatened") %>%
  group_by(taxon) %>% 
  summarise(n_sp = n_distinct(species)) %>% 
  pull(n_sp)

# Define and prepare contingency table
cont_table <- matrix(c(n_sp_pca, n_sp_kba, n_sp_ipa),
                     nrow = 4)

rownames(cont_table) <- c("Birds", "Herptiles", "Mammals", "Plants")
colnames(cont_table) <- c("PCA", "KBA", "IPA")

cont_table

# Run the chi square test
chisq.test(cont_table)

# Get the cramer's v value
cramerV(cont_table, conf = 0.95)

# Run pairwise tests
pairwise.prop.test(cont_table[, c("PCA", "KBA")], p.adjust.method = "bonferroni")  # PCA vs KBA
pairwise.prop.test(cont_table[, c("PCA", "IPA")], p.adjust.method = "bonferroni")  # PCA vs IPA
pairwise.prop.test(cont_table[, c("KBA", "IPA")], p.adjust.method = "bonferroni")  # KBA vs IPA
