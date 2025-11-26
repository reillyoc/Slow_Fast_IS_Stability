# Bayesian Mediation Analysis

# Author(s): Reilly O'Connor
# Version: 2025-07-23

# Load Pkgs
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(beepr)
library(easystats)
library(ggpubr)
library(extrafont)
library(visreg)
library(MuMIn)
library(tidyverse)
library(ape)
library(caper)
library(phylolm)
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)

#Set BRMS cmdstanr
options(brms.backend = "cmdstanr")

#load data
#Mammal Phylogeny from Beccari et al 2024
df_tree <- readRDS("../R/Data/PhilogenyReady.rds")
# df_lpi <- read.csv("../R/Data/lpi mammal population cv detrended resolved taxonomy lifespan windows.csv", header = T)
df_lpi <- read.csv("../R/Data/lpi mammal population cv detrended resolved taxonomy.csv", header = T)

unique(df_lpi$ID)
unique(df_lpi$Binomial)

df_mammal_groups <- read.csv("../R/Data/Suggested Mammal Groups.csv", header = T)

unique(df_mammal_groups$Order)
unique(df_mammal_groups$Family)

#Hatton Data
df_metab <- read.csv("../R/Data/slow_fast_metabolism resolved taxonomy.csv", header = T)

df_growth <- read.csv("../R/Data/slow_fast_growth resolved taxonomy.csv", header = T)

df_growth_pc <- read.csv("../R/Data/slow_fast_pc1_rmax_tax_resolved.csv", header = T)

#Beccari et al 2024
df_mammal_fs <- read.csv("../R/Data/Beccari Fast-Slow resolved taxonomy.csv", header = T)

#Body Size Data for CV - See Mean Body Size Calculation.R Script
df_body_mass <- read.csv("../R/Data/mean species body size.csv")

##### Code #####
#Subset final species list
df_lpi_species_list <- df_lpi %>% 
  dplyr::select(Binomial, Order, Family, Genus, Species)

df_lpi_species_list <- unique(df_lpi_species_list)

#write.csv(df_lpi_species_list, "../R/Data/LPI CV Final Species List.csv")

df_lpi$Latitude <- as.numeric(df_lpi$Latitude)
df_lpi$Longitude <- as.numeric(df_lpi$Longitude)

df_lpi_groups <- merge(df_lpi, df_mammal_groups, by = "Binomial")

df_lpi_mammal_bs <- merge(df_lpi_groups, df_body_mass, by = "GBIF_ID") 

unique(df_lpi_mammal_bs$ID)

df_lpi_mbs <- df_lpi_mammal_bs %>%
  dplyr::select(Group, GBIF_ID, species_names, Species.y, Binomial, Common_name, cv_window, cv_value, k_cv_value, mean_body_mass, Latitude, Longitude) %>%
  rename(Species = Species.y)


fit_log_smr <- lm(log10(Mass_Spec_Meta_Watt) ~ log10(Mass_g),
                  data = df_metab)

df_metab$MISMR <- residuals(fit_log_smr)


##### CV ~ Body Size + PC1 #####
unique(df_lpi_mbs$GBIF_ID)

df_mammal_cv_bs <- df_lpi_mbs %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value)) %>%
  group_by(Group, Binomial, GBIF_ID, cv_window) %>%
  summarise(mean_cv = mean(cv_value, na.rm = T),
            mean_log_cv = mean(log_cv, na.rm = T),
            # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
            mean_logit_k_cv = mean(logit_k_cv, na.rm = T),
            mean_k_cv = mean(k_cv_value, na.rm = T),
            mean_body_mass = mean(mean_body_mass, na.rm = T),
            mean_latitude = mean(Latitude, na.rm = T)) %>%
  mutate(log_mean_body_mass = log10(mean_body_mass),
         abs_latitude = abs(mean_latitude))

unique(df_mammal_cv_bs$GBIF_ID)

df_lpi_fs <- merge(df_lpi_mbs, df_mammal_fs, by = "GBIF_ID")

df_mammal_cv_bs_fs <- df_lpi_fs %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value)) %>%
  # filter(cv_window == 5) %>%
  group_by(Group, GBIF_ID, cv_window) %>%
  summarise(mean_cv = mean(cv_value, na.rm = T),
            mean_log_cv = mean(log_cv, na.rm = T),
            # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
            mean_logit_k_cv = mean(logit_k_cv, na.rm = T),
            PC1 = mean(PC1, na.rm = T),
            PC2 = mean(PC2, na.rm = T),
            mean_k_cv = mean(k_cv_value, na.rm = T),
            mean_body_mass = mean(mean_body_mass, na.rm = T),
            mean_latitude = mean(Latitude, na.rm = T)) %>%
  mutate(log_mean_body_mass = log10(mean_body_mass),
         abs_latitude = abs(mean_latitude))

unique(df_mammal_cv_bs_fs$GBIF_ID)

df_mammal_cv_bs_fs_w <- df_mammal_cv_bs_fs %>%
  filter(cv_window == 5)

# Model 2: k_cv as a function of rmax and body size
bf_kcv_bs <- bf(mean_k_cv ~ log_mean_body_mass + PC1)

priors_kcv_bs <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("exponential(1)", class = "phi")
  # set_prior("cauchy(0, 1)", class = "sigma")
)

brms_k_cv_bs <- brm(
  formula = bf_kcv_bs,
  data = df_mammal_cv_bs_fs_w,
  prior = priors_kcv_bs,
  family = Beta(link = "logit"),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95)
)

loo_k_cv_bs <- loo(brms_k_cv_bs)
loo_k_cv_bs

pp_check(brms_k_cv_bs, ndraws = 100)
summary(brms_k_cv_bs)
plot(brms_k_cv_bs)

# Create new data for prediction
df_pred_bs <- tibble(
  log_mean_body_mass = seq(min(df_mammal_cv_bs_fs_w$log_mean_body_mass),
                           max(df_mammal_cv_bs_fs_w$log_mean_body_mass), length.out = 100),
  PC1 = mean(df_mammal_cv_bs_fs_w$PC1),
  PC2 = mean(df_mammal_cv_bs_fs_w$PC2)
)

# Get expected posterior predictions
epred_draws_bs <- add_epred_draws(brms_k_cv_bs, newdata = df_pred_bs)

gg_epred_bs <- ggplot(epred_draws_bs, aes(x = log_mean_body_mass, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.5, linewidth = 1.5, color = "black", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_epred_bs

df_pred_pc1 <- tibble(
  log_mean_body_mass = mean(df_mammal_cv_bs_fs$log_mean_body_mass),
  PC1 = seq(min(df_mammal_cv_bs_fs$PC1),
            max(df_mammal_cv_bs_fs$PC1), length.out = 100),
  PC2 = mean(df_mammal_cv_bs_fs$PC2)
)

# Get expected posterior predictions
epred_draws_pc1 <- add_epred_draws(brms_k_cv_bs, newdata = df_pred_pc1)

gg_epred_pc1 <- ggplot(epred_draws_pc1, aes(x = PC1, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.5, linewidth = 1.5, color = "black", fill = "grey") +
  labs(x = "PC1 (Slow-Fast Traits)", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_epred_pc1

gg_epred_bs_pc1 <- plot_grid(gg_epred_bs, gg_epred_pc1, nrow = 1, align = 'hv')
gg_epred_bs_pc1

ggsave("../R/Figures/Figure 2 - Bayesian CV - Body Size.jpeg", plot = gg_epred_bs_pc1, width = 8, height = 4)


###### CV, rmax, Body Size, PC1 ######
df_growth_pc_cv <- merge(df_lpi_mbs, df_growth_pc, by = "GBIF_ID")
unique(df_growth_pc_cv$Binomial)

df_mammal_cv_rmax_pc <- df_growth_pc_cv %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value),
         se_logit_k_cv = 0.10 / (1 - k_cv_value)) %>%
  group_by(Group, Binomial, GBIF_ID, cv_window) %>%
  summarise(mean_cv = mean(cv_value, na.rm = T),
          mean_log_cv = mean(log_cv, na.rm = T),
          # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
          k_cv = mean(k_cv_value, na.rm = T),
          mean_logit_k_cv = mean(logit_k_cv, na.rm = T),
          se_logit_k_cv = mean(se_logit_k_cv, na.rm = TRUE),
          mean_body_mass = mean(mean_body_mass, na.rm = T),
          PC1 = mean(PC1, na.rm = T),
          PC2 = mean(PC2, na.rm = T),
          rmax = mean(rmax, na.rm = T),
          residual_rmax = mean(residual_rmax, na.rm = T)) %>%
  ungroup() %>%
  mutate(log_mean_body_mass = log10(mean_body_mass),
         log_rmax = log10(rmax),
         log_rmax_se = standard_error(log_rmax)) #%>%
#filter(! (Group == "Bats" | Group == "Marine Mammals"))

df_mammal_cv_fs <- df_lpi_fs %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value)) %>%
  # filter(cv_window == 5) %>%
  group_by(Group, Binomial.x, GBIF_ID, cv_window) %>%
  reframe(mean_cv = mean(cv_value, na.rm = T),
          mean_k_cv = mean(k_cv_value, na.rm = T),
          mean_log_cv = mean(log_cv, na.rm = T),
          # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
          mean_logit_k_cv = mean(logit_k_cv, na.rm = T),
          mean_body_mass = mean(mean_body_mass, na.rm = T),
          mean_body_mass_fs = mean(bm, na.rm = T),
          PC1 = mean(PC1, na.rm = T),
          PC2 = mean(PC2, na.rm = T)) %>%
  rename(Binomial = Binomial.x) # %>%
  # filter(! (Group == "Bats" | Group == "Marine Mammals"))

##### 5 Year Window - Body Mass vs CV #####
df_cv_5 <- df_mammal_cv_rmax_pc  %>% filter(cv_window == 5) #%>%
#filter(! (Group == "Marine Mammals" | Group == "Bats"))

unique(df_cv_5$Binomial)

# df_cv_5_fs <- df_mammal_cv_fs %>% filter(cv_window == 5) #%>% 
# #filter(! (Group == "Marine Mammals" | Group == "Bats"))

# Model 1: rmax as a function of body size and PC1
bs_rmax_1 <- bf(log_rmax | se(log_rmax_se, sigma = TRUE) ~ log_mean_body_mass + PC1)

# Model 2: k_cv as a function of rmax and body size
bf_kcv_rbs_pc <- bf(k_cv ~ log_rmax + log_mean_body_mass + PC1)

bf_kcv_rbs <- bf(k_cv ~ log_rmax + log_mean_body_mass)

bf_kcv_rmax <- bf(k_cv ~ log_rmax)

priors_updated <- c(
  # Priors for rmax model
  set_prior("normal(-0.24, 0.02)", class = "b", coef = "log_mean_body_mass", resp = "logrmax"),
  set_prior("normal(0.15, 0.02)", class = "b", coef = "PC1", resp = "logrmax"),
  set_prior("normal(0.60, 0.05)", class = "Intercept", resp = "logrmax"),
  set_prior("cauchy(0, 1)", class = "sigma", resp = "logrmax"),
  
  # Priors for k_cv model
  set_prior("normal(0, 1)", class = "b", resp = "kcv"),
  set_prior("exponential(1)", class = "phi", resp = "kcv")
  # set_prior("cauchy(0, 1)", class = "sigma", resp = "kcv")
)


fit_joint_rmax_rbs_pc <- brm(
  bs_rmax_1 + bf_kcv_rbs_pc + set_rescor(FALSE),
  data = df_cv_5,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit_joint_rmax_rbs <- brm(
  bs_rmax_1 + bf_kcv_rbs + set_rescor(FALSE),
  data = df_cv_5,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit_joint_rmax <- brm(
  bs_rmax_1 + bf_kcv_rmax + set_rescor(FALSE),
  data = df_cv_5,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
                ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

loo_bf_kcv_rmax_rbs_pc <- loo(fit_joint_rmax_rbs_pc)
loo_bf_kcv_rmax_rbs_pc

loo_bf_kcv_rmax_rbs <- loo(fit_joint_rmax_rbs)
loo_bf_kcv_rmax_rbs

loo_bf_kcv_rmax <- loo(fit_joint_rmax)
loo_bf_kcv_rmax

loo_compare(loo_bf_kcv_rmax, loo_bf_kcv_rmax_rbs, loo_bf_kcv_rmax_rbs_pc)

pp_check(fit_joint_rmax, resp = "logrmax")
pp_check(fit_joint_rmax, resp = "kcv")

summary(fit_joint_rmax)

# For overall model (R-squared for each response variable)
bayes_R2(fit_joint_rmax)

# To get R-squared for a specific response variable (e.g., 'logrmax' or 'kcv')
bayes_R2(fit_joint_rmax, resp = "logrmax")
bayes_R2(fit_joint_rmax, resp = "kcv")

options(scipen = 999)

posterior_summary(fit_joint_rmax, probs = c(0.055, 0.945))

plot(fit_joint_rmax)

mcmc_plot(fit_joint_rmax, 
          type = "intervals", 
          prob = 0.5, 
          prob_outer = .89,
          point_est = "median")

conditional_effects(fit_joint_rmax_rbs, "log_mean_body_mass", resp = "logrmax")
conditional_effects(fit_joint_rmax_rbs, "PC1", resp = "logrmax")
conditional_effects(fit_joint_rmax_rbs, "log_rmax", resp = "kcv")
conditional_effects(fit_joint_rmax_rbs, "log_mean_body_mass", resp = "kcv")

posterior <- as_draws_df(fit_joint_rmax_rbs)

posterior_ind <- posterior %>%
  mutate(
    indirect_bs = b_logrmax_log_mean_body_mass * b_kcv_log_rmax + b_kcv_log_mean_body_mass,
    indirect_pc1 = b_logrmax_PC1 * b_kcv_log_rmax)

rmax_median_ci <- posterior_ind %>%
  summarise(
    lower = quantile(b_kcv_log_rmax, 0.025),
    median = median(b_kcv_log_rmax),
    upper = quantile(b_kcv_log_rmax, 0.975),
    pred = "rmax")

bs_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_bs, 0.025),
    median = median(indirect_bs),
    upper = quantile(indirect_bs, 0.975),
    pred = "body size")

pc1_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_pc1, 0.025),
    median = median(indirect_pc1),
    upper = quantile(indirect_pc1, 0.975),
    pred = "fast-slow continuum")

# pc2_median_ci<- posterior_ind %>%
#   summarise(
#     lower = quantile(indirect_pc2, 0.025),
#     median = median(indirect_pc2),
#     upper = quantile(indirect_pc2, 0.975),
#     pred = "lifetime reproductive effort")

df_median_cis <- rbind(rmax_median_ci, bs_median_ci, pc1_median_ci)

ggplot(df_median_cis, aes(x = median, y = pred)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_point(shape = 21, stroke = 0.5, fill = "aliceblue", color = "black", size = 3) +
  # xlim(-0.3, 1.0) +
  labs(x = "Direct or Indirect Effect on CV", 
       y = "Covariate") +
  theme_classic(base_size = 14)


df_pred_pc1 <- tibble(
  log_mean_body_mass = mean(df_cv_5$log_mean_body_mass, na.rm = TRUE),
  PC1 = mean(df_cv_5$PC1, na.rm = TRUE),
  PC2 = mean(df_cv_5$PC2, na.rm = TRUE),
  log_rmax = seq(min(df_cv_5$log_rmax, na.rm = TRUE),
                 max(df_cv_5$log_rmax, na.rm = TRUE), length.out = 100),
  log_rmax_se = mean(df_cv_5$log_rmax_se, na.rm = TRUE)  # if needed
)

epred_pc1_cv <- add_epred_draws(fit_joint_rmax_rbs,
                                newdata = df_pred_pc1,
                                resp = "kcv")

gg_rmax_cv <- ggplot(epred_pc1_cv, aes(x = log_rmax, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.5, linewidth = 1.5,
                  color = "black", fill = "grey") +
  labs(x = "log rmax", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_rmax_cv


ggsave("../R/Figures/Figure 3 - Bayesian CV - Rmax.jpeg", plot = gg_rmax_cv, width = 5, height = 5)


###### Metabolism - CV #####
df_metcv <- merge(df_cv_5, df_metab, by = c("GBIF_ID")) 

# Model 1: rmax as a function of body size and PC1
bs_rmax_metab <- bf(log_rmax | se(log_rmax_se, sigma = TRUE) ~ log_mean_body_mass + PC1 + MISMR)

# Model 2: k_cv as a function of rmax and body size
bf_kcv_rbs_pc_met <- bf(k_cv ~ log_rmax + log_mean_body_mass + PC1 + MISMR)

bf_kcv_rbs_pc <- bf(k_cv ~ log_rmax + log_mean_body_mass + PC1)

bf_kcv_rbs <- bf(k_cv ~ log_rmax + log_mean_body_mass)

bf_kcv_rmax <- bf(k_cv ~ log_rmax)

priors_updated <- c(
  # Priors for rmax model
  set_prior("normal(-0.25, 0.1)", class = "b", coef = "log_mean_body_mass", resp = "logrmax"),
  set_prior("normal(0.15, 0.1)", class = "b", coef = "PC1", resp = "logrmax"),
  set_prior("normal(0.15, 0.15)", class = "b", coef = "MISMR", resp = "logrmax"),
  set_prior("normal(0.67, 0.2)", class = "Intercept", resp = "logrmax"),
  set_prior("cauchy(0, 1)", class = "sigma", resp = "logrmax"),
  
  # Priors for k_cv model
  set_prior("normal(0, 1)", class = "b", resp = "kcv"),
  set_prior("exponential(1)", class = "phi", resp = "kcv")
  # set_prior("cauchy(0, 1)", class = "sigma", resp = "kcv")
)


fit_joint_rmax_met_rbs_pc_met <- brm(
  bs_rmax_metab + bf_kcv_rbs_pc_met + set_rescor(FALSE),
  data = df_metcv,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit_joint_rmax_met_rbs_pc <- brm(
  bs_rmax_metab + bf_kcv_rbs_pc + set_rescor(FALSE),
  data = df_metcv,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit_joint_rmax_met_rbs <- brm(
  bs_rmax_metab + bf_kcv_rbs + set_rescor(FALSE),
  data = df_metcv,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

fit_joint_rmax_met <- brm(
  bs_rmax_metab + bf_kcv_rmax + set_rescor(FALSE),
  data = df_metcv,
  prior = priors_updated,
  family = list(gaussian(), Beta(link = "logit")
  ),  
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

loo_bf_kcv_rmax_met_rbs_pc_met <- loo(fit_joint_rmax_met_rbs_pc_met)
loo_bf_kcv_rmax_met_rbs_pc_met

loo_bf_kcv_rmax_met_rbs_pc <- loo(fit_joint_rmax_met_rbs_pc)
loo_bf_kcv_rmax_met_rbs_pc

loo_bf_kcv_rmax_met_rbs <- loo(fit_joint_rmax_met_rbs)
loo_bf_kcv_rmax_met_rbs

loo_bf_kcv_rmax_met <- loo(fit_joint_rmax_met)
loo_bf_kcv_rmax_met

loo_compare(loo_bf_kcv_rmax_met, loo_bf_kcv_rmax_met_rbs, loo_bf_kcv_rmax_met_rbs_pc, loo_bf_kcv_rmax_met_rbs_pc_met)

pp_check(fit_joint_rmax_met, resp = "logrmax")
pp_check(fit_joint_rmax_met, resp = "kcv")

summary(fit_joint_rmax_met)

# For overall model (R-squared for each response variable)
bayes_R2(fit_joint_rmax_met)

# To get R-squared for a specific response variable (e.g., 'logrmax' or 'kcv')
bayes_R2(fit_joint_rmax_met, resp = "logrmax")
bayes_R2(fit_joint_rmax_met, resp = "kcv")

options(scipen = 999)

posterior_summary(fit_joint_rmax_met, probs = c(0.055, 0.945))

plot(fit_joint_rmax_met)

mcmc_plot(fit_joint_rmax_met, 
          type = "intervals", 
          prob = 0.5, 
          prob_outer = .89,
          point_est = "median")

conditional_effects(fit_joint_rmax_met, "log_mean_body_mass", resp = "logrmax")
conditional_effects(fit_joint_rmax_met, "PC1", resp = "logrmax")
conditional_effects(fit_joint_rmax_met, "MISMR", resp = "logrmax")
conditional_effects(fit_joint_rmax_met, "log_rmax", resp = "kcv")
# conditional_effects(fit_joint_rmax_met, "log_mean_body_mass", resp = "kcv")


posterior <- as_draws_df(fit_joint_rmax_met)

posterior_ind <- posterior %>%
  mutate(
    indirect_bs = b_logrmax_log_mean_body_mass * b_kcv_log_rmax, #+ b_kcv_log_mean_body_mass,
    indirect_pc1 = b_logrmax_PC1 * b_kcv_log_rmax,
    indirect_mismr = b_logrmax_MISMR * b_kcv_log_rmax)

rmax_median_ci <- posterior_ind %>%
  summarise(
    lower = quantile(b_kcv_log_rmax, 0.025),
    median = median(b_kcv_log_rmax),
    upper = quantile(b_kcv_log_rmax, 0.975),
    pred = "rmax")

bs_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_bs, 0.025),
    median = median(indirect_bs),
    upper = quantile(indirect_bs, 0.975),
    pred = "body size")

pc1_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_pc1, 0.025),
    median = median(indirect_pc1),
    upper = quantile(indirect_pc1, 0.975),
    pred = "fast-slow continuum")

met_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_mismr, 0.025),
    median = median(indirect_mismr),
    upper = quantile(indirect_mismr, 0.975),
    pred = "Mass-Independent Metabolic Rate")

df_median_cis <- rbind(rmax_median_ci, bs_median_ci, pc1_median_ci, met_median_ci)

ggplot(df_median_cis, aes(x = median, y = pred)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_point(shape = 21, stroke = 0.5, fill = "aliceblue", color = "black", size = 3) +
  # xlim(-0.3, 1.0) +
  labs(x = "Direct or Indirect Effect on CV", 
       y = "Covariate") +
  theme_classic(base_size = 14)


df_pred_pc1 <- tibble(
  log_mean_body_mass = mean(df_metcv$log_mean_body_mass, na.rm = TRUE),
  PC1 = mean(df_metcv$PC1, na.rm = TRUE),
  PC2 = mean(df_metcv$PC2, na.rm = TRUE),
  log_rmax = seq(min(df_metcv$log_rmax, na.rm = TRUE),
                 max(df_metcv$log_rmax, na.rm = TRUE), length.out = 100),
  log_rmax_se = mean(df_metcv$log_rmax_se, na.rm = TRUE)  # if needed
)

epred_pc1_cv <- add_epred_draws(fit_joint_rmax_rbs,
                                newdata = df_pred_pc1,
                                resp = "kcv")

gg_rmax_cv <- ggplot(epred_pc1_cv, aes(x = log_rmax, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.5, linewidth = 1.5,
                  color = "black", fill = "grey") +
  labs(x = "log rmax", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_rmax_cv







