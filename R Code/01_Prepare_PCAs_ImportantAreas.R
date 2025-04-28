# Load libraries
library(sf)
library(tidyverse)
library(osmdata)
library(giscoR)
library(stringdist)
library(terra)
library(units)
sf_use_s2(FALSE)

############################ ---- Important Plant Areas ---- ############################# 
# Read in the IPA data
ipa_files <- list.files("data/raw_data/IPAs/", pattern = ".shp$", full.names = TRUE)
ipa_shapefiles <- lapply(ipa_files, function(x) {
  dat <- st_read(x)
  dat <- st_transform(dat, crs = "EPSG:4326")
  return(dat)
})

# Make the column names and order the same for each file
ipa_shapefiles[[1]] <- ipa_shapefiles[[1]] %>% 
  rename(Name = IPAName) %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[2]] <- ipa_shapefiles[[2]] %>% 
  rename(Name = IPA_NAME, Country = COUNTRYNAM) %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[3]] <- ipa_shapefiles[[3]] %>% 
  rename(Name = IPA_name) %>% 
  mutate(Country = "United Kingdom") %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[4]] <- ipa_shapefiles[[4]] %>% 
  rename(Name = IPAName) %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[5]] <- ipa_shapefiles[[5]] %>% 
  rename(Name = IPA_NAME, Country = COUNTRYNAM) %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[6]] <- ipa_shapefiles[[6]] %>% 
  rename(Name = IPAName) %>% 
  dplyr::select(Country, Name)
ipa_shapefiles[[7]] <- ipa_shapefiles[[7]] %>% 
  rename(Name = IPA_NAME, Country = COUNTRYNAM) %>% 
  dplyr::select(Country, Name)

ipa_shapefiles <- do.call(rbind, ipa_shapefiles)

# Split into point and polygon geometries - we will only use the polygons
ipa_point <- ipa_shapefiles[st_geometry_type(ipa_shapefiles) == "POINT", ] %>% 
  st_zm(drop = TRUE) %>% # drop z geometries
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

ipa_poly <- ipa_shapefiles[st_geometry_type(ipa_shapefiles) != "POINT", ]

######################## ---- Tropical Important Plant Areas ---- ######################## 
# Read in the TIPAs data from the RBG Kew TIPAs Explorer
tipa_point <- st_read("data/raw_data/TIPAs/TIPA_Composite_POINT.shp") %>% 
  filter(Country != "New Guinea") %>% 
  st_zm(drop = TRUE) %>% # drop z geometries
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

tipa_poly <- st_read("data/raw_data/TIPAs/TIPA_Composite_POLYGON.shp") %>% 
  dplyr::select(Country, Name) %>% 
  filter(Country != "New Guinea")

# Read in the Colombia TIPAs data from Kor & Diazgranados (2023)
col_tipas <- list.files("data/raw_data/TIPAs/COL_top_priority/", pattern = ".shp$", 
                        full.names = TRUE)
col_tipas <- lapply(col_tipas, st_read) 
col_tipas <- do.call(rbind, col_tipas) %>% 
  mutate(Country = "Colombia", Name = NA) %>% 
  dplyr::select(Country, Name)

# Combine all TIPA and IPA polygons into one file
ipa_data <- rbind(ipa_poly, tipa_poly, col_tipas) %>% 
  mutate(Country = ifelse(Country == "BVI", "British Virgin Islands", Country)) %>% 
  st_zm(drop = TRUE) %>% # drop z geometries
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

# Check for and fix invalid geometries
invalid_geom <- !st_is_valid(ipa_data) # identify invalid geometries 
table(invalid_geom) 
ipa_data$geometry[invalid_geom] <- st_make_valid(ipa_data$geometry[invalid_geom]) # fix the invalid geometries

# Check geometry types - sometimes fixing geometries can throw different geometry types in, and we only want polygons
geom_types <- st_geometry_type(ipa_data)
unique(geom_types) # there seems to be a GEOMETRYCOLLECTION record

problematic_rows <- which(geom_types != "POLYGON" & geom_types != "MULTIPOLYGON") # identify problematic rows with the undesired geometries
ipa_data$geometry[problematic_rows] <- st_collection_extract(ipa_data$geometry[problematic_rows], "POLYGON") # extract the polygons from these problematic rows

# Dissolve geometries of sites with the same name and from the same country
ipa_data_no_na <- ipa_data[!is.na(ipa_data$Name), ] # step 1: filter out rows with NA in the Name column

ipa_data_dissolved <- ipa_data_no_na %>% # step 2: dissolve geometries by country and site names
  group_by(Country, Name) %>%
  summarize(geometry = st_union(geometry), .groups = 'drop')

ipa_data_with_na <- ipa_data[is.na(ipa_data$Name), ] # step 3: find rows with NA in the Name column

ipa_data_dissolved <- rbind(ipa_data_dissolved, ipa_data_with_na) # step 4: combine into one dataframe

# Use st_intersection to allocate each row a country name
countries <- sort(unique(ipa_data_dissolved$Country))# vector of country names
write.csv(countries, "data/output_data/country_names.csv", row.names = FALSE)

countries <- read.csv("data/output_data/country_names.csv")[,1]
countries <- gisco_get_countries(country = countries, year = "2020", resolution = "01") %>% # country boundaries
  dplyr::select(ISO3_CODE, NAME_ENGL) %>% 
  st_transform(., crs = "+proj=moll") %>%  # reproject to Mollweide equal-area projection
  mutate(IUCN_name = NAME_ENGL) %>% 
  mutate(IUCN_name = case_when(
    IUCN_name == "Bolivia" ~ "Bolivia, Plurinational State of",
    IUCN_name == "Palestine" ~ "Palestine, State of",
    IUCN_name == "Syria" ~ "Syrian Arab Republic",
    IUCN_name == "British Virgin Islands" ~ "Virgin Islands, British",
    TRUE ~ IUCN_name
  ))

countries <- countries %>%
  arrange(IUCN_name)

write_sf(countries, "data/raw_data/country_boundaries/countries.shp")

# Read in the countries data
countries <- st_read("data/raw_data/country_boundaries/countries.shp")

# Crop the data to only include geometries inside the countries - this neatly separates records 
# if they are in multiple countries
ipa_data_in_countries <- st_intersection(ipa_data_dissolved, countries)
ipa_points_in_countries <- st_intersection(ipa_point, countries) # n = 2
tipa_points_in_countries <- st_intersection(tipa_point, countries) # n = 179

# 179 point IPAs, 1178 polygon IPA... 13.2% of IPAs are points in studied countries

# Save the prepared IPA data
write_sf(ipa_data_in_countries, "data/output_data/prepared_areas/important_plant_areas_prepared.shp")

# Read in the IPA dates of identification csv
ipa_dates <- read.csv("data/raw_data/IPAs/IPA_country_programmes_year.csv") %>% 
  dplyr::select(Country, Date.of.publication.in)

# Combine the IPA data with the years
ipa_data_in_countries <- st_read("data/output_data/prepared_areas/important_plant_areas_prepared.shp")
ipa_data_in_countries <- left_join(ipa_data_in_countries, ipa_dates, by = "Country")

# Overwrite the prepared IPA data
write_sf(ipa_data_in_countries, "data/output_data/prepared_areas/important_plant_areas_prepared.shp")

########################### ---- Key Biodiversity Areas ---- ############################# 
kba_polygons <- st_read("data/raw_data/KBAsGlobal_2024_June_02_project/KBAsGlobal_2024_June_02_POL.shp")
kba_points <- st_read("data/raw_data/KBAsGlobal_2024_June_02_project/KBAsGlobal_2024_June_02_PNT.shp")

kba_polygons <- kba_polygons %>% 
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection
kba_points <- kba_points %>% 
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

# Crop the data to only include geometries inside the countries - this neatly separates records 
# if they are in multiple countries
kba_polygons_in_countries <- st_intersection(kba_polygons, countries) # n = 3167
kba_points_in_countries <- st_intersection(kba_points, countries) # n = 111, so proportion that are points = 3.38%

write_sf(kba_polygons_in_countries, "data/output_data/prepared_areas/key_biodiversity_areas_prepared.shp")

############################ ---- Important Bird Areas ---- ############################## 
iba_polygons <- st_read("data/raw_data/IBAsGlobal_2024_March_01/IBAsGlobal_2024_March_POL_01.shp")
iba_points <- st_read("data/raw_data/IBAsGlobal_2024_March_01/IBAsGlobal_2024_March_PNT_01.shp")

iba_polygons <- iba_polygons %>% 
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection
iba_points <- iba_points %>% 
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

invalid_geom <- !st_is_valid(iba_polygons) # identify invalid geometries
table(invalid_geom) # n = 9
iba_polygons$geometry[invalid_geom] <- st_make_valid(iba_polygons$geometry[invalid_geom]) # fix only the invalid geometries

iba_polygons_in_countries <- st_intersection(iba_polygons, countries) # n = 2526
iba_points_in_countries <- st_intersection(iba_points, countries) # n = 62, so proportion that are points = 2.39%

write_sf(iba_polygons_in_countries, "data/output_data/prepared_areas/important_bird_areas_prepared.shp")

############################## ---- KBA/IBA Overlaps ---- ################################ 
iba_area <- st_area(st_union(iba_polygons_in_countries))
kba_area <- st_area(st_union(kba_polygons_in_countries))
ibas_in_kbas_area <- st_area(st_union(st_intersection(iba_polygons_in_countries, 
                                                      kba_polygons_in_countries)))

ibas_in_kbas_area/iba_area*100 # 92.7% of the area identified as IBAs is also within KBAs

######################### ---- Protected & Conserved Areas ---- ########################## 
# Read in the WDPA/WDOECM data
wdpa_polygons <- rbind(st_read("data/raw_data/WDPA_WDOECM/Part_1/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-polygons.shp"),
                       st_read("data/raw_data/WDPA_WDOECM/Part_2/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-polygons.shp"),
                       st_read("data/raw_data/WDPA_WDOECM/Part_3/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-polygons.shp"))

wdpa_points <- rbind(st_read("data/raw_data/WDPA_WDOECM/Part_1/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-points.shp"),
                     st_read("data/raw_data/WDPA_WDOECM/Part_2/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-points.shp"),
                     st_read("data/raw_data/WDPA_WDOECM/Part_3/WDPA_WDOECM_Aug2024_Public_e3a33d4e19793936b59ba7b3c658ef1e9f2e950a625b6b79fd1957eb12891b00_shp-points.shp"))

# Filter and prepare the WDPA/WDOECM data
wdpa_polygons <- wdpa_polygons %>% 
  filter(MARINE != 2) %>% # remove entirely marine PAs (n = 635)
  filter(STATUS != "Proposed") %>% # remove Proposed PAs (n = 212), which are PAs in the process of being designated
  mutate(STATUS_YR = if_else(STATUS_YR == 0, NA_integer_, STATUS_YR)) %>% # assign NA to records without a STATUS_YR
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

wdpa_points <- wdpa_points %>% 
  filter(MARINE != 2) %>% # remove entirely marine PAs (n = 57)
  filter(STATUS != "Proposed") %>% # remove Proposed PAs (n = 17), which are PAs in the process of being designated
  mutate(STATUS_YR = if_else(STATUS_YR == 0, NA_integer_, STATUS_YR)) %>% # assign NA to records without a STATUS_YR
  st_transform(., crs = "+proj=moll") # reproject to Mollweide equal-area projection

# 7.9% of all WDPA/WDOECM geometries were points 

points_without_REP_AREA <- wdpa_points %>% 
  filter(REP_AREA == 0) # 3.5% of WDPA/WDOECM geometries were points that were removed, the rest were points we buffered

# Turn the points into polygons using a buffer equal to the reported area (REP_AREA) of the PA 
wdpa_points <- wdpa_points %>% 
  filter(REP_AREA != 0) # remove points without a reported area

wdpa_points <- wdpa_points %>%
  mutate(radius = sqrt(REP_AREA / pi) * 1000) %>% # convert area from km to meters, so the radius is in meters
  mutate(geometry = st_buffer(geometry, dist = radius)) %>% # buffer each point by the calculated radius (in meters)
  dplyr::select(-radius)

# Visconti et al (2013) suggest 1) clipping buffered centroids known to be only marine or 
# terrestrial to their respective realm (we do this when we crop to country boundaries); 
# 2) clip buffered centroids known from only one country to its border...
unique(wdpa_points$ISO3) # no points have multiple ISO3 codes, so we can clip them by their sole country boundary

parent_iso3_codes <- unique(wdpa_points$PARENT_ISO) # get the parent ISO3 codes of the countries with point PAs

wdpa_points_intersect <- data.frame()
for (i in seq_along(parent_iso3_codes)) {
  dat <- wdpa_points %>% 
    filter(PARENT_ISO == parent_iso3_codes[i]) %>% # filter by ISO3 code
    st_intersection(., subset(countries, ISO3_CODE == parent_iso3_codes[i])) # find the intersection between those filtered points and the boundary of the corresponding country
  
  wdpa_points_intersect <- rbind(wdpa_points_intersect, dat) # there are fewer rows in this output compared to wdpa_points as some points are outside of the country's borders
}

# Combine the polygons and buffered, intersected points into one data frame
wdpa_polygons <- wdpa_polygons %>% 
  dplyr::select(-GIS_M_AREA, -GIS_AREA) # remove these columns as they are not relevant and are not in the points data frame, ibhibiting rbind

wdpa_points_intersect <- wdpa_points_intersect %>% 
  dplyr::select(-ISO3_CODE, -NAME_ENGL, -IUCN_name) # remove these columns as they were unneccesarily imported during the st_intersection

wdpa_polygons <- rbind(wdpa_polygons, wdpa_points_intersect)

# If STATUS_YR contains missing data, randomly assign a year from another PA in the same country
missing_years <- subset(wdpa_polygons, is.na(STATUS_YR)) # records with missing years for reference

wdpa_polygons <- wdpa_polygons %>% 
  group_by(PARENT_ISO) %>%
  mutate(STATUS_YR = ifelse(is.na(STATUS_YR), sample(STATUS_YR[!is.na(STATUS_YR)], size = 1), STATUS_YR)) %>% # randomly assign year from a PA with the same PARENT_ISO
  ungroup() %>% 
  st_zm(drop = TRUE) # drop z geometries

# Check number of OECMs vs PAs
table(wdpa_polygons$PA_DEF) # 409 OECMs; 16617 PAs

# Check for and fix invalid geometries
invalid_geom <- !st_is_valid(wdpa_polygons) # identify invalid geometries
table(invalid_geom) # n = 4
wdpa_polygons$geometry[invalid_geom] <- st_make_valid(wdpa_polygons$geometry[invalid_geom]) # fix only the invalid geometries

# Check geometry types - sometimes fixing geometries can throw different geometry types in, and we only want polygons
geom_types <- st_geometry_type(wdpa_polygons)
unique(geom_types) # good news! there's only POLYGON and MULTIPOLYGON geometries

# Crop the data to only include geometries inside the countries - this neatly separates records 
# if they are in multiple countries
wdpa_data_in_countries <- st_intersection(wdpa_polygons, countries)

geom_types <- st_geometry_type(wdpa_data_in_countries)
unique(geom_types)
wdpa_data_in_countries <- wdpa_data_in_countries[st_geometry_type(wdpa_data_in_countries) %in% c("POLYGON", "MULTIPOLYGON"), ] # some geometries change following intersection

# Save the prepared WDPA/WDOECM data
write_sf(wdpa_data_in_countries, "data/output_data/prepared_areas/wdpa_wdoecm_prepared.shp")

# India and China have incomplete datasets on the WDPA, but we can obtain polygons from
# OpenStreetMap using the osmdata package
# 1) India
india_bb <- getbb("India") # get India's bounding box coordinates

india_pa <- india_bb %>%
  opq() %>% # initialise overpass query
  add_osm_features(features = c( # select features
    "boundary" = "protected_area",
    "leisure" = "nature_reserve",
    "boundary" = "national_park"
  )) %>%
  osmdata_sf() # return as sf object

# Prepare India's protected area polygons
india_poly <- india_pa$osm_polygons %>% 
  st_transform(., crs = st_crs(countries)) # reproject polygons to Mollweide projection

invalid_geom <- !st_is_valid(india_poly) # identify invalid geometries
table(invalid_geom) # n = 0, no need to fix any

india_poly <- india_poly %>% 
  st_intersection(., subset(countries, NAME_ENGL == "India")) # find polygons within India's borders

india_poly <- india_poly %>% 
  dplyr::select(name, name.en, start_date, ref.WDPA, wikidata) # retain relevant columns only

# Prepare India's protected area multipolygons
india_multipoly <- india_pa$osm_multipolygons %>% 
  st_transform(., crs = st_crs(countries)) # reproject polygons to Mollweide projection

invalid_geom <- !st_is_valid(india_multipoly) # identify invalid geometries
table(invalid_geom) # n = 75
india_multipoly$geometry[invalid_geom] <- st_make_valid(india_multipoly$geometry[invalid_geom]) # fix only the invalid geometries

india_multipoly <- india_multipoly %>% 
  st_intersection(., subset(countries, NAME_ENGL == "India")) # find multipolygons within India's borders

india_multipoly <- india_multipoly %>% 
  dplyr::select(name, name.en, start_date, ref.WDPA, wikidata) # retain relevant columns only

# Combine into one dataframe
india_poly <- rbind(india_poly, india_multipoly) %>% 
  mutate_at(vars(-geometry), ~ na_if(., "")) %>% # replace "" with NA in each column, excluding the geometry column
  mutate(name = ifelse(!is.na(name.en), name.en, name)) %>% # replace name with name.en if name.en is not missing
  dplyr::select(name) %>%  # no longer need name.en
  st_intersection(., subset(countries, NAME_ENGL == "India")) # find polygons within India's borders

# 2) China
china_bb <- getbb("China") # get China's bounding box coordinates

china_pa <- china_bb %>%
  opq() %>% # initialise overpass query
  add_osm_features(features = c( # select features
    "boundary" = "protected_area",
    "leisure" = "nature_reserve",
    "boundary" = "national_park"
  )) %>%
  osmdata_sf() # return as sf object

# Prepare China's protected area polygons
china_poly <- china_pa$osm_polygons %>% 
  st_transform(., crs = st_crs(countries)) # reproject polygons to Mollweide projection

invalid_geom <- !st_is_valid(china_poly) # identify invalid geometries
table(invalid_geom) # n = 0, no need to fix any

china_poly <- china_poly %>% 
  dplyr::select(name) # retain relevant columns only

# Prepare China's protected area multipolygons
china_multipoly <- china_pa$osm_multipolygons %>% 
  st_transform(., crs = st_crs(countries)) # reproject polygons to Mollweide projection

invalid_geom <- !st_is_valid(china_multipoly) # identify invalid geometries
table(invalid_geom) # n = 312
china_multipoly$geometry[invalid_geom] <- st_make_valid(china_multipoly$geometry[invalid_geom]) # fix only the invalid geometries

china_multipoly <- china_multipoly %>% 
  dplyr::select(name) # retain relevant columns only

# Combine into one dataframe
china_poly <- rbind(china_poly, china_multipoly) %>% 
  mutate_at(vars(-geometry), ~ na_if(., "")) %>% # replace "" with NA in each column, excluding the geometry column
  dplyr::select(name) %>% # no longer need name.en
  st_intersection(., subset(countries, NAME_ENGL == "China")) # find polygons within China's borders

india_china_poly <- rbind(india_poly, china_poly)

# write_sf(india_china_poly, "data/output_data/prepared_areas/IND_CHN_osm_prepared.shp")

# # WE EXCLUDE INDIA AND CHINA FROM THE THROUGH TIME PRINCIPAL COMPONENT ANALYSIS DUE TO THE
# # UNCERTAINTY ASSOCIATED OF THEIR YEAR OF IDENTIFICATION, BUT BELOW IS THE LONG PROCESS WE
# # UNDERTOOK TO GLEAN YEARS OF IDENTIFICATION
# india_poly <- st_read("data/raw_data/OSM/india_protected_areas.shp")
# 
# india_poly$start_date <- str_remove(india_poly$start_date, "-07-07")
# india_poly$start_date <- str_remove(india_poly$start_date, "-02-09")
# india_poly$name <- str_replace_all(india_poly$name, "Aravalli Biodiversity Park", "Aravali Biodiversity Park")
# 
# india_poly <- india_poly %>% 
#   mutate(start_date = as.numeric(start_date), # make date numeric, so we can retain it using min/max in the summarise below
#          ref_WDPA = as.numeric(ref_WDPA)) # make ref.WDPA numeric, so we can retain it using min/max in the summarise below
# 
# # Dissolve geometries of sites with the same name
# india_poly_no_na <- india_poly[!is.na(india_poly$name), ] # step 1: filter out rows with NA in the Name column - we won't use sites without a name
# 
# india_poly_dissolved <- india_poly_no_na %>% # step 2: dissolve geometries by name
#   group_by(name) %>%
#   summarise(
#     geometry = st_union(geometry),  # union the geometries
#     start_date = min(start_date), # earliest start_year
#     ref_WDPA = min(ref_WDPA), # WDPA code
#     wikidata_first = first(wikidata),
#     wikidata_last = last(wikidata)) %>% # retain as much wikidata as possible in case the dissolved geometries have different wikidata
#   relocate(geometry, .after = wikidata_last)
# 
# # Add year to india_poly based on corresponding name in india_pa_dates
# india_pa_dates <- read.csv("data/raw_data/india_pa_dates.csv") # a dataset from the National Data Repository, Directorate General of Hydrocarbons, Government of India
# 
# common_rows <- india_poly_dissolved %>%
#   filter(name %in% india_pa_dates$Name) # find common rows with identical names in both india_poly and india_pa_dates
# 
# non_common_names <- setdiff(india_poly_dissolved$name, common_rows$name) # find names of rows that are not common
# non_common_rows <- india_poly_dissolved %>% # find non-common rows
#   filter(name %in% non_common_names)
# 
# find_similar_names <- function(df1, df2, threshold) { # function to find similar names between dataframes
#   results <- data.frame()
#   
#   for (name1 in df1$name) { # loop through each name
#     distances <- stringdist(name1, df2$Name, method = "lv") # compute the string distances between the current name and all names in the second data frame
#     min_distance <- min(distances) # find the minimum distance from the computed distances
#     
#     if (min_distance < threshold) { # check if the minimum distance is below the threshold
#       similar_names <- df2$Name[distances == min_distance] # identify names in df2 that have the minimum distance
#       results <- rbind(results, data.frame( # combine into a dataframe
#         osm_name = name1,
#         similar_name = similar_names,
#         distance = min_distance
#       ))
#     }
#   }
#   results
# }
# 
# # Apply the function with a threshold
# similar_names_results <- find_similar_names(non_common_rows, india_pa_dates, threshold = 3.5)
# 
# similar_names_results <- similar_names_results %>% # manually filter out names that are clearly wrong
#   filter(!osm_name %in% c("Kailam WLS", "Land", "land", "Nangal WLS", "Trees", "Nandhaur WLS",
#                           "Nellai WLS", "Rampara WLS", "SITE", "Tamhini WLS", "Udanti WLS"))
# 
# similar_rows <- non_common_rows %>%
#   filter(name %in% similar_names_results$osm_name) %>% 
#   left_join(similar_names_results, by = c("name" = "osm_name")) %>%
#   mutate(name = ifelse(is.na(similar_name), name, similar_name)) %>%
#   select(-similar_name, -distance)
# 
# common_rows <- rbind(common_rows, similar_rows) # add the rows identified as being similar in name
# non_common_rows <- non_common_rows %>% 
#   filter(!(name %in% similar_names_results$osm_name))
#   
# # There are still 364 rows we haven't matched to india_pa_dates!
# india_pa_dates_not_matched <- india_pa_dates %>% 
#   filter(!(Name %in% common_rows$name))
# 
# wikidata_rows <- non_common_rows %>% 
#   filter(!is.na(wikidata_first)) %>% # find the rows with wikidata - we will manually match dates
#   dplyr::select(name, wikidata_first, wikidata_last) %>% 
#   st_drop_geometry()
# 
# write.csv(wikidata_rows, "data/output_data/wikidata_india_rows.csv", row.names = FALSE)
# 
# # Let's assign dates to the rows we've matched. 
# # 1) Common Rows
# common_rows <- common_rows %>% # match these with dates from india_pa_dates 
#   dplyr::select(name, geometry) %>% 
#   group_by(name) %>% # combine the geometries of rows with the same name
#   summarise(geometry = st_union(geometry)) %>% 
#   left_join(., india_pa_dates, by = c("name" = "Name")) %>% 
#   relocate(geometry, .after = "Year")
# 
# # Sort the data frame by 'name' and 'year' (ascending order)
# common_rows <- common_rows[order(common_rows$name, common_rows$Year), ]
# 
# # Remove duplicates, keeping only the first occurrence of each 'name'
# common_rows <- common_rows[!duplicated(common_rows$name), ]
# 
# # 2) Wikidata Rows
# wikidata_rows <- read.csv("data/output_data/wikidata_india_rows_with_dates.csv") %>% 
#   left_join(., non_common_rows, by = "name") %>% # join to add geometry back in
#   dplyr::select(name, Year, geometry) %>% 
#   st_as_sf(., crs = st_crs(common_rows))
# 
# # 3) Non Common Rows
# non_common_rows <- non_common_rows %>% 
#   filter(!(name %in% wikidata_rows$name)) %>% 
#   rename(Year = start_date) %>%
#   mutate(Year = ifelse(wikidata_last == "Q235878", 1971, Year)) %>% # manually assign correct year to any records with wikidata_last values
#   mutate(Year = ifelse(wikidata_last == "Q7764533", 1987, Year))
# 
# non_common_rows_manual <- non_common_rows %>% 
#   st_drop_geometry() %>% 
#   dplyr::select(name, Year)
# 
# write.csv(non_common_rows_manual, "data/output_data/non_common_india_rows.csv", row.names = FALSE)
# 
# # Combine the dated rows together
# india_poly <- rbind(common_rows, wikidata_rows)
# 
# write_sf(india_poly, "data/output_data/prepared_areas/osm_india.shp")

########################### ---- Area in Geographic Space ---- ########################### 
# Read in the country names vector and then the country boundaries
countries <- st_read("data/raw_data/country_boundaries/countries.shp")

# Read in the prepared PCA, KBA and IPA data
pa_data <- st_read("data/output_data/prepared_areas/wdpa_wdoecm_prepared.shp")
osm_pa_data <- st_read("data/output_data/prepared_areas/IND_CHN_osm_prepared.shp") %>% 
  rename("NAME" = "name")
pa_data <- bind_rows(pa_data, osm_pa_data) # combine the PCA data into one dataframe

ipa_data <- st_read("data/output_data/prepared_areas/important_plant_areas_prepared.shp") %>% 
  rename(IUCN_name = IUCN_nm)
kba_data <- st_read("data/output_data/prepared_areas/key_biodiversity_areas_prepared.shp")

# Break the PCA, KBA and IPA data down into lists, where each element is a country
pa_list <- split(pa_data, pa_data$IUCN_name)
kba_list <- split(kba_data, kba_data$IUCN_name)
ipa_list <- split(ipa_data, ipa_data$IUCN_name)

# Calculate area of each country's PCAs, KBAs and IPAs in geographic space
area_list <- list()

for (i in 1:nrow(countries)) {
  country_area <- st_area(st_union(countries[i, ])) %>% 
    set_units(km2) %>% 
    drop_units()
  pa_area <- st_area(st_union(pa_list[[i]])) %>% 
    set_units(km2) %>% 
    drop_units()
  kba_area <- st_area(st_union(kba_list[[i]])) %>% 
    set_units(km2) %>% 
    drop_units()
  ipa_area <- st_area(st_union(ipa_list[[i]])) %>% 
    set_units(km2) %>% 
    drop_units()
  
  area_list[[i]] <- data.frame(Country = country_area, 
                               PCA = pa_area, 
                               KBA = kba_area, 
                               IPA = ipa_area) %>% 
    mutate(PCA_pct = PCA/Country*100,
           KBA_pct = KBA/Country*100,
           IPA_pct = IPA/Country*100,
           Country_Name = countries$IUCN_name[i]) %>% 
    dplyr::select(-Country) %>% 
    relocate(Country_Name)
}

geo_cov <- do.call(rbind, area_list)

# Find mean across countries and add to dataframe
means <- c("Mean", as.numeric(colMeans(geo_cov[, -1])))
geo_cov <- rbind(geo_cov, means)

# Find the global values (different from the mean)
country_area <- st_area(st_union(countries)) %>% 
  set_units(km2) %>% 
  drop_units()
pa_area <- st_area(st_union(pa_data)) %>% 
  set_units(km2) %>% 
  drop_units()
kba_area <- st_area(st_union(kba_data)) %>% 
  set_units(km2) %>% 
  drop_units()
ipa_area <- st_area(st_union(ipa_data)) %>% 
  set_units(km2) %>% 
  drop_units()

global_area <- data.frame(Country = country_area,
                          PCA = pa_area, 
                          KBA = kba_area, 
                          IPA = ipa_area) %>% 
  mutate(PCA_pct = PCA/Country*100,
         KBA_pct = KBA/Country*100,
         IPA_pct = IPA/Country*100,
         Country_Name = "Overall") %>% 
  dplyr::select(-Country) %>% 
  relocate(Country_Name)

geo_cov <- rbind(geo_cov, global_area)

# Make values numeric and round to 2 decimal places
geo_cov <- geo_cov %>% 
  mutate(across(2:7, ~ as.numeric(.))) %>% 
  mutate(across(2:7, ~ round(., 2))) 

write.csv(geo_cov, "data/output_data/geo_space_area.csv", row.names = FALSE)

# Also calculate global PCA area in 1960
a <- pa_data %>% 
  filter(STATUS_YR <= 1960)
b <- pa_data %>% 
  filter(is.na(STATUS_YR))
