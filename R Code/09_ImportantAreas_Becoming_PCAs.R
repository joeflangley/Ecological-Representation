# Load libraries
library(tidyverse)
library(sf)
library(raster)
library(terra)
library(fasterize)
library(epitools)
library(rcompanion)

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
any(is.na(ipa_data$Dt_f_p_))
any(is.na(kba_data$year))

# Combine the KBA and IPA data into one dataframe
important_areas <- rbind(kba_data, ipa_data) 

######################### ---- Rasterize PCAs, KBAs & IPAs ---- ##########################
# Generate a 5000m Molleweide spatial raster template
r <- raster()
x <- projectRaster(r, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
r2 <- x
res(r2) <- 5000
raster_template <- resample(x, r2)
spat_rast_template <- rast(raster_template) # turn from raster to spatial raster

# Function to rasterize polygons without missing small polygons
rasterize_process <- function(input_data, template_data) {
  vv <- vect(input_data)
  cnt <- centroids(vv)
  x <- terra::rasterize(vv, template_data, touches = T, field = STATUS_YR, background = 0)
  y <- terra::rasterize(cnt, template_data, touches = T, field = STATUS_YR, background = 0)
  cover(x, y)
}

# Apply for each country boundary
country_rasters <- lapply(1:nrow(countries), function(i){
  rasterize_process(countries[i, ], spat_rast_template)
})

# Apply for the PCA & Important Area polygons
pca_rasters <- rasterize_process(pa_data, spat_rast_template)
important_area_rasters <- rasterize_process(important_areas, spat_rast_template)

all_rasters <- list()
for (i in 1:nrow(countries)) {
  # Get the country's spatial polygon
  country <- countries[i, ]

  # Crop the country/PCA/KBA/IPA rasters by the country's boundary
  country_raster <- crop(country_rasters[[i]], country)
  pca_raster_national <- crop(pca_rasters, country)
  important_area_raster_national <- crop(important_area_rasters, country)

  # Compute country-specific masks for PCAs, KBAs, and IPAs
  country_pca_mask <- country_raster * pca_raster_national
  country_important_area_mask <- country_raster * important_area_raster_national
  
  # Combine raster values from both raster outputs into a dataframe
  raster_df <- data.frame(PCA = values(pca_raster),
                          KBA_IPA = values(important_areas_raster)) %>% 
    drop_na() %>% # remove cells outside the country's boundary
    filter(PCA == 0 | PCA > 2010) # exclude PCAs that were designated in 2010 or before, while keeping cells that aren't PCAs (PCA = 0)
  
  # Condition 1: Outside KBA/IPA in 2010 and outside PCA by 2024
  condition_1 <- (raster_df$KBA_IPA == 0) & (raster_df$PCA == 0)
  
  # Condition 2: Inside KBA/IPA in 2010 but outside PCA by 2024
  condition_2 <- (raster_df$KBA_IPA != 0) & (raster_df$KBA_IPA <= 2010) & (raster_df$PCA == 0)
  
  # Condition 3: Outside KBA/IPA in 2010 but inside PCA by 2024
  condition_3 <- (raster_df$KBA_IPA == 0) & (raster_df$PCA != 0)
  
  # Condition 4: Inside KBA/IPA in 2010 but inside PCA by 2024
  condition_4 <- (raster_df$KBA_IPA != 0) & (raster_df$KBA_IPA <= 2010) & (raster_df$PCA != 0)
  
  # Count the number of cells for each condition
  count_condition_1 <- sum(condition_1)
  count_condition_2 <- sum(condition_2)
  count_condition_3 <- sum(condition_3)
  count_condition_4 <- sum(condition_4)
  
  # Create the contingency table
  contingency_table <- matrix(c(count_condition_1, count_condition_2,
                                count_condition_3, count_condition_4),
                              nrow = 2, byrow = FALSE)
  
  rownames(contingency_table) <- c("Outside Important Area in 2010", "Inside Important Area in 2010")
  colnames(contingency_table) <- c("Outside PCA in 2024", "Inside PCA in 2024")
  
  all_rasters[[i]] <- list(df = raster_df, contingency_table = contingency_table)
}

# 1) By Country ----
chisq_outputs <- list()

for (i in 1:nrow(countries)) {
  
  if (any(all_rasters[[i]]$contingency_table == 0)) {
    res <- chisq.test(all_rasters[[i]]$contingency_table)
    cramer <- cramerV(all_rasters[[i]]$contingency_table, conf = 0.95)
    critical_value <- qchisq(p = res$p.value, df = res$parameter, lower.tail = FALSE)

    dat <- cbind.data.frame(
      country = countries$IUCN_name[i],
      `X^2` = res$statistic, # chi squared statistic
      p = res$p.value, # chi squared p-value
      std_res = res$stdres[2,2], # standardized residual value for "Inside Important Area in 2010, Inside PCA in 2024"
      critical_val = critical_value,
      cramers_v = cramer, # cramer's v
      odds_ratio = NA,
      odds_lower = NA,
      odds_upper = NA,
      odds_p = NA
    )
    
  } else {
    res <- chisq.test(all_rasters[[i]]$contingency_table)
    cramer <- cramerV(all_rasters[[i]]$contingency_table, conf = 0.95)
    critical_value <- qchisq(p = res$p.value, df = res$parameter, lower.tail = FALSE)
    odds <- oddsratio(all_rasters[[i]]$contingency_table, correction = TRUE)
    
    dat <- cbind.data.frame(
      country = countries$IUCN_name[i],
      `X^2` = res$statistic, # chi squared statistic
      p = res$p.value, # chi squared p-value
      std_res = res$stdres[2,2], # standardized residual value for "Inside Important Area in 2010, Inside PCA in 2024"
      critical_val = critical_value,
      cramers_v = cramer, # cramer's v
      odds_ratio = odds$measure[2,1], # odds ratio
      odds_lower = odds$measure[2,2], # odds ratio lower bound
      odds_upper = odds$measure[2,3], # odds ratio upper bound
      odds_p = odds$p.value[2,3] # odds ratio chi-squared p-value
    )
  }

  chisq_outputs[[i]] <- dat
}

chisq_outputs_df <- do.call(rbind, chisq_outputs)
rownames(chisq_outputs_df) <- NULL

# Replace NaN with NA across all columns of the dataframe
chisq_outputs_df <- chisq_outputs_df %>%
  mutate(across(everything(), ~ replace(., is.nan(.), NA)))

# 2) Global ----

# Sum the contingency tables in the all_rasters list to create a global contingency table
global_cont_table <- Reduce("+", lapply(all_rasters, function(x) x$contingency_table))

# Save the contingency table
write.csv(global_cont_table, "data/output_data/chi_squared_cont_table.csv")

# Run the global chi-squared analyses
res <- chisq.test(global_cont_table)
cramer <- cramerV(global_cont_table, conf = 0.95)
critical_value <- qchisq(p = res$p.value, df = res$parameter, lower.tail = FALSE)
odds <- oddsratio(global_cont_table, correction = TRUE)
global_chisq_output_df <- cbind.data.frame(
  country = "Global",
  `X^2` = res$statistic, # chi squared statistic
  p = res$p.value, # chi squared p-value
  std_res = res$stdres[2,2], # standardized residual value for "Inside Important Area in 2010, Inside PCA in 2024"
  critical_val = critical_value,
  cramers_v = cramer, # cramer's v
  odds_ratio = odds$measure[2,1], # odds ratio
  odds_lower = odds$measure[2,2], # odds ratio lower bound
  odds_upper = odds$measure[2,3], # odds ratio upper bound
  odds_p = odds$p.value[2,3] # odds ratio chi-squared p-value
)

rownames(global_chisq_output_df) <- NULL

# Add the global output values to the nationmal output values and save the data
chisq_outputs_df <- rbind(chisq_outputs_df, global_chisq_output_df)
write.csv(chisq_outputs_df, "data/output_data/chi_squared_outut.csv", row.names = FALSE)

