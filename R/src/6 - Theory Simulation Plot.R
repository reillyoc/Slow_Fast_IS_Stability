# Description

# Author(s): 
# Version: YYYY-MM-DD

# Load Pkgs
library(tidyverse)
library(cowplot)

source("../Slow_Fast_IS_Stability/R/src/0 - Functions.R")

# load data
df_cr_is_flux <- read.csv("../Slow_Fast_IS_Stability/Julia/Outputs/CR_IS_Flux_rmax.csv")
cr_dir <- "/Users/reillyoconnor/Desktop/R Projects/Slow_Fast_IS_Stability/Julia/Outputs/CR"
cr_csv_files <- list.files(path = cr_dir, pattern = "*.csv", full.names = TRUE)

df_cr <- cr_csv_files %>%
  lapply(read.csv) %>%
  bind_rows() %>%
  mutate(rmax_increase = as.character(rmax_increase)) %>%
  group_by(rmax_increase) %>%
  summarise(se = standard_error(cv),
            mean = mean(mean, na.rm = T),
            sd = mean(std, na.rm = T),
            sd_cv = sd(cv, na.rm = T),
            cv = mean(cv, na.rm = T)) %>%
  mutate(rmax_increase = as.numeric(rmax_increase),
         k_cv = sqrt(cv^2/(1+cv^2)))

##### Consumer-Resource Interaction Strength, Growth Potential, CV2 Plots #####
gg_cr_is_flux <- ggplot(df_cr_is_flux, aes(x = rmax, y = interaction_strength_flux)) +
  #geom_smooth(se = T, aes(color = NA), color = "black", linetype = "solid", linewidth = 1.5) +
  geom_point(alpha = 0.5,  fill = "black", stroke = 0.3, size = 3, shape = 21) +
  #scale_fill_gradient(low = "#3CA373", high = "#EB8F00") +
  ylab("Interaction Strength (Flux)") +
  xlab("Growth Potential (rmax)") +
  theme_classic(base_size = 16) +
  theme(
        legend.position = "none",
        text = element_text(family = "Arial")) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) 
  

gg_cr_is_flux

ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 2 - Consumer IS-Flux-Rmax Theory.jpeg", plot = gg_cr_is_flux, width = 8, height = 6, dpi = 300)

gg_cr_cv_rmax <- ggplot(df_cr, aes(x = rmax_increase, y = k_cv)) +
  #geom_smooth(se = T, aes(color = NA), color = "black", linetype = "solid", linewidth = 1.5) +
  geom_point(fill = "black", stroke = 0.3, size = 3, shape = 21) +
  #scale_fill_gradient(low = "#3CA373", high = "#EB8F00") +
  ylab("Temporal Variabiltiy (CV)") +
  xlab("Growth Potential (rmax)") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    text = element_text(family = "Arial")) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) 

gg_cr_cv_rmax

ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 2 - Consumer CV-Rmax Theory.jpeg", plot = gg_cr_cv_rmax, width = 8, height = 6, dpi = 300)

gg_is_cv_grid <- plot_grid(gg_cr_is_flux, gg_cr_cv_rmax, nrow = 1, align = "hv")
gg_is_cv_grid

ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 1 - Consumer IS, CV vs Rmax theory.jpeg", plot = gg_is_cv_grid, width = 8, height = 4, dpi = 300)


