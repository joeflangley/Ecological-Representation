# Load libraries
# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(lsmeans)
library(ggh4x)

# Read in the geographical and environmental coverage datasets
geo_cov <- read.csv("data_old/output_data/geographical_coverage_per_country.csv") 
env_cov <- read.csv("data_old/output_data/pa_bpa/pct_overlaps_with_background.csv")

geo_cov <- read.csv("data/output_data/geo_space_area.csv") 
env_cov <- read.csv("data/output_data/env_space_mcp_area.csv")

# Pivot both of them into longer formats
geo_cov <- geo_cov %>% 
  rename(country = Country_Name) %>% # give country column the same name as env_cov
  dplyr::select(country, PCA_pct, KBA_pct, IPA_pct) %>% # retain the relevant columns
  rename("PCA" = PCA_pct, "KBA" = KBA_pct, "IPA" = IPA_pct) %>% # give the percentage columns the same names as env_cov
  pivot_longer(PCA:IPA, values_to = "geo_pct", names_to = "category") # make longer

env_cov <- env_cov %>% 
  dplyr::select(country, PCA, KBA, IPA) %>% # retain the relevant columns
  pivot_longer(PCA:IPA, values_to = "env_pct", names_to = "category") # make longer

# Combine into one dataframe
df <- left_join(env_cov, geo_cov, by = c("country", "category")) %>% 
  filter(country != "Mean") # remove mean values

n_distinct(df$country) # Check we have all studied countries

################################ ---- Linear Models ---- ################################# 
# Subset the data by category and log transform the geo and env percentages
pca_df <- subset(df, category == "PCA") %>% 
  mutate(log_env_pct = log(env_pct),
         log_geo_pct = log(geo_pct)) # to remove 0s

kba_df <- subset(df, category == "KBA") %>% 
  mutate(log_env_pct = log(env_pct),
         log_geo_pct = log(geo_pct)) # to be consistent with pa_df

ipa_df <- subset(df, category == "IPA") %>% 
  mutate(log_env_pct = log(env_pct),
         log_geo_pct = log(geo_pct)) # to be consistent with pa_df

# Fit linear models 
pca_lm <- lm(log_env_pct ~ log_geo_pct, data = pca_df)
kba_lm <- lm(log_env_pct ~ log_geo_pct, data = kba_df)
ipa_lm <- lm(log_env_pct ~ log_geo_pct, data = ipa_df)

# Add predictions from respective models to each dataframe
pca_df$pred <- predict(pca_lm)
kba_df$pred <- predict(kba_lm)
ipa_df$pred <- predict(ipa_lm)

# Combine into one dataframe
all_df <- rbind(pca_df, kba_df, ipa_df)

# Plot the linear models on one figure
env_geo_slopes <- ggplot(data = all_df,
                         mapping = aes(x = log_geo_pct, 
                                       y = log_env_pct,
                                       color = category))+
  geom_point(shape = 1)+
  geom_line(mapping = aes(y = pred), linewidth = 1)+
  scale_color_manual(values = c("PCA" = "#0072B2",
                                "KBA" = "#E69F00",
                                "IPA" = "#009E73")) +
  labs(color = NULL,
       x = "Log geographical space (%)",
       y = "Log environmental space (%)")+
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 9.5),
        legend.text = element_text(size = 10),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.35))

# Save the plot
ggsave(plot = env_geo_slopes, "outputs/Environmental_Expansion/geo_vs_env_space_slopes.png",
       device = "png",
       width = 180,
       height = 140,
       dpi = 600,
       units = "mm")

# Test for differences between slopes of PA, KBA and IPA linear models
lm_all <- lm(log_env_pct ~ log_geo_pct*category, data = all_df) # linear model with interaction of geo_pct by conservation network
anova(lm_all) # interaction term does not have a significant effect
lm_lst <- lstrends(lm_all, "category", var="log_geo_pct") # compute slopes of geo_pct for each network
pairs(lm_lst) # compare slopes - no significant differences

summary(pca_lm)
summary(kba_lm)
summary(ipa_lm)

saveRDS(env_geo_slopes, "outputs/geo_vs_env_space_slopes.rds")
density_plot <- readRDS("data/output_data/dens_plot.rds")
plot(density_plot)

library(ggpubr)
fig <- ggarrange(env_geo_slopes, density_plot, 
                 ncol = 1,
                 labels = c("A", "B"),
                 hjust = -1,
                 heights = c(1, 1.75))

ggsave(plot = fig, "outputs/Environmental_Expansion/slopes_and_densityplot.png", 
       width = 180, height = 280, units = "mm", dpi = 600,
       device = "png")

