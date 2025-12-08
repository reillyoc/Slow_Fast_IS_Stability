#Calculating a Mean Body Size from All Databases

#Author(s): Reilly O'Connor
#Version: 2023-11-21

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)

#load all data
df_metab <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/slow_fast_metabolism resolved taxonomy.csv", header = T)
df_growth <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/slow_fast_growth resolved taxonomy.csv", header = T)
df_body_sizes <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/Amniote database resolved taxonomy.csv", header = T)

##### Code #####
#reduce dataframes to just body mass and species... 

df_metab_bs <- df_metab %>% filter(Major_taxa == "Mammal") %>%
  dplyr::select(Species, GBIF_ID, Mass_g)
df_growth_bs <- df_growth %>% filter(Major_taxa == "Mammal") %>%
  dplyr::select(Species, GBIF_ID, Mass_g)
df_body_sizes <- df_body_sizes %>%
  dplyr::select(Species, GBIF_ID, Mass_g)

df_bs <- rbind(df_metab_bs, df_growth_bs, df_body_sizes)
unique(df_bs$GBIF_ID)

df_bs_sp <- df_bs %>% group_by(Species, GBIF_ID) %>%
  summarize(mean_body_mass = mean(Mass_g))

unique(df_bs_sp$Species)
unique(df_bs_sp$GBIF_ID)

options(scipen = 999)

write.csv(df_bs_sp, "../Slow_Fast_IS_Stability/R/Outputs/mean species body size.csv")


