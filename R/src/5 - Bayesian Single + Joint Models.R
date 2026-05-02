# Exploring Relationships between rmax, body size and fast-slow life history strategies and metabolism in mammals

#Author(s): Reilly O'Connor
#Version: 2024-06-11

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(Hmisc)
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
library(emmeans)

source("../Slow_Fast_IS_Stability/R/src/0 - Functions.R")

#Set BRMS cmdstanr
options(brms.backend = "cmdstanr")

#load all data
#Mammal Phylogeny from Beccari et al 2024
df_tree <- readRDS("../Slow_Fast_IS_Stability/R/Data/PhilogenyReady.rds")

#Living Planet Index CV Calculations
df_lpi <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/lpi mammal population cv detrended resolved taxonomy.csv", header = T)
df_lpi <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/lpi mammal population cv detrended resolved taxonomy lifespan windows.csv", header = T)

unique(df_lpi$ID)
unique(df_lpi$Binomial)

#CV species to subset
df_lpi_sp <- df_lpi %>%
  distinct(Binomial, GBIF_ID)

df_mammal_groups <- read.csv("../Slow_Fast_IS_Stability/R/Data/Suggested Mammal Groups.csv", header = T)

unique(df_mammal_groups$Order)
unique(df_mammal_groups$Family)

#Body Size Data for CV - See Mean Body Size Calculation.R Script
df_body_mass <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/mean species body size.csv")

#Hatton Data
df_metab <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/slow_fast_metabolism resolved taxonomy.csv", header = T)
df_growth <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/slow_fast_growth resolved taxonomy.csv", header = T)
df_growth_count <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/slow_fast_growth resolved taxonomy - count method.csv", header = T)

#Beccari et al 2024
df_mammal_fs <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/Beccari Fast-Slow resolved taxonomy.csv", header = T)

##### Code #####
##### Relationship between rmax ~ Body Size + PC1 #####

#Line up Species names with Phylogeny for rmax, PC1 data frame
df_rmax_fs <- merge(df_mammal_fs, df_growth, by = c("GBIF_ID")) %>%
  dplyr::select(- X.x, -X.1, -X.y, - Species.y, - Binomial, Binomial = Species.x)

df_rmax_fs_sub <- df_rmax_fs %>%
  anti_join(df_lpi_sp, by = "GBIF_ID") %>%
  group_by(Binomial, GBIF_ID) %>%
  summarise(rmax = mean(rmax), 
         Mass_g = mean(Mass_g),
         PC1 = mean(PC1),
         PC2 = mean(PC2)) %>%
  ungroup() %>%
  mutate(log_rmax = log10(rmax),
         log_mean_body_mass = log10(Mass_g))

global_mean_log_body_mass <- mean(df_rmax_fs_sub$log_mean_body_mass)
global_sd_log_body_mass <- sd(df_rmax_fs_sub$log_mean_body_mass)

df_rmax_fs_sub <- df_rmax_fs_sub %>%
  mutate(sc_log_mean_body_mass = (log_mean_body_mass - global_mean_log_body_mass)/global_sd_log_body_mass)

unique(df_rmax_fs_sub$Binomial)

#Match Phylogeny up
sp  <- unique(df_rmax_fs_sub$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_sub) <- df_rmax_fs_sub$Binomial

#Phylogeny correlation matrix
rmax_fs_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_sub$phylo <- factor(df_rmax_fs_sub$Binomial, levels = rownames(rmax_fs_corrma))

#Model Priors
priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("student_t(3, 0, 10)", class = "sigma")
)

priors_null <- c(
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("student_t(3, 0, 10)", class = "sigma")
)

brms_rmax_null <- brm(
  formula = log_rmax ~ 1 + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_sub,
  data2 = list(A = rmax_fs_corrma),
  prior = priors_null,
  family = gaussian(),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_rmax_null <- brms::loo(brms_rmax_null)
loo_rmax_null

brms_rmax_bs <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_sub,
  data2 = list(A = rmax_fs_corrma),
  prior = priors,
  family = gaussian(),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_rmax_bs <- brms::loo(brms_rmax_bs)
loo_rmax_bs

brms_rmax_bs_pc1 <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_sub,
  data2 = list(A = rmax_fs_corrma),
  prior = priors,
  family = gaussian(),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_rmax_bs_pc1 <- brms::loo(brms_rmax_bs_pc1)
loo_rmax_bs_pc1

r2_bayes(brms_rmax_null)
r2_bayes(brms_rmax_bs)
r2_bayes(brms_rmax_bs_pc1)

brms_rmax_bs_pc1
pp_check(brms_rmax_bs_pc1, ndraws = 100)
bayes_R2(brms_rmax_bs_pc1, probs = c(0.055, 0.945))
summary(brms_rmax_bs_pc1, prob = 0.89)
plot(brms_rmax_bs_pc1)

#Estimates of body size unscaled
draws <- as_draws_df(brms_rmax_bs_pc1) %>%
  mutate(
    beta_x  = b_sc_log_mean_body_mass / sd(df_rmax_fs_sub$log_mean_body_mass),
    alpha_x = b_Intercept - b_sc_log_mean_body_mass * mean(df_rmax_fs_sub$log_mean_body_mass) / sd(df_rmax_fs_sub$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .055), u95 = quantile(beta_x, .945))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .055), u95 = quantile(alpha_x, .945))


#Create new data for prediction
df_pred_bs <- tibble(
  sc_log_mean_body_mass = seq(min(df_rmax_fs_sub$sc_log_mean_body_mass),
                           max(df_rmax_fs_sub$sc_log_mean_body_mass), length.out = 100),
  PC1 = mean(df_rmax_fs_sub$PC1),
  log_mean_body_mass = seq(min(df_rmax_fs_sub$log_mean_body_mass),
                           max(df_rmax_fs_sub$log_mean_body_mass), length.out = 100))


#Get expected posterior predictions
epred_draws_bs <- add_epred_draws(brms_rmax_bs_pc1, newdata = df_pred_bs, re_formula = NA)

gg_epred_bs <- ggplot(epred_draws_bs, aes(x = log_mean_body_mass, y = .epred)) +
  geom_point(data = df_rmax_fs_sub, aes(x = log_mean_body_mass, y = log_rmax), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted log(rmax)") +
  # scale_fill_gradient(low = "#3CA373", high = "#EB8F00") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2.2, 1.4), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_bs

#Create new data for prediction
df_pred_pc1 <- tibble(
  sc_log_mean_body_mass = mean(df_rmax_fs_sub$sc_log_mean_body_mass),
  PC1 = seq(min(df_rmax_fs_sub$PC1),
            max(df_rmax_fs_sub$PC1), length.out = 100)
)

#Get expected posterior predictions
epred_draws_pc1 <- add_epred_draws(brms_rmax_bs_pc1, newdata = df_pred_pc1, re_formula = NA)

gg_epred_pc1 <- ggplot(epred_draws_pc1, aes(x = PC1, y = .epred)) +
  geom_point(data = df_rmax_fs_sub, aes(x = PC1, y = log_rmax), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "PC1", y = "Predicted log(rmax)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2.25, 1.4), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_pc1

gg_epred_bs_pc1 <- plot_grid(gg_epred_bs, gg_epred_pc1, nrow = 1, align = 'hv', labels = c("a", "b"))
gg_epred_bs_pc1

# ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 1 - Bayesian Rmax - Body Size.jpeg", plot = gg_epred_bs_pc1, width = 8, height = 4)

###### Extract Posteriors for use as Priors in CV ~ rmax Models later ######
post_fixef_bspc <- data.frame(fixef(brms_rmax_bs_pc1))
post_sigma_bspc <- data.frame(posterior_summary(brms_rmax_bs_pc1, variable = "sigma"))
post_sd_phy_bspc <- data.frame(posterior_summary(brms_rmax_bs_pc1, variable = "sd_phylo__Intercept"))

##### Relationship between CV ~ Body Size + PC1 #####
###### Format and Combine CV, BS, FS Data #####
df_lpi$Latitude <- as.numeric(df_lpi$Latitude)
df_lpi$Longitude <- as.numeric(df_lpi$Longitude)

df_lpi_groups <- merge(df_lpi, df_mammal_groups, by = "Binomial")

df_lpi_mammal_bs <- merge(df_lpi_groups, df_body_mass, by = "GBIF_ID") 

df_lpi_groups_tbl_ave <- df_lpi_mammal_bs %>%
  filter(cv_window == 5) %>%
  group_by(Group, Binomial) %>%
  summarise(mean_body_mass = mean(mean_body_mass),
            k_cv_value = mean(k_cv_value)) %>%
  group_by(Group) %>%
  summarise(mean_body_size = mean(mean_body_mass),
            sd_body_size = sd(mean_body_mass),
            se_body_size = standard_error(mean_body_mass),
            mean_kcv = mean(k_cv_value),
            sd_kcv = sd(k_cv_value),
            se_kcv = standard_error(k_cv_value),
            count = n())

unique(df_lpi_mammal_bs$ID)

df_lpi_mbs <- df_lpi_mammal_bs %>%
  dplyr::select(Group, GBIF_ID, ID, species_names, Species.y, Binomial, Common_name, cv_window, cv_value, k_cv_value, mean_body_mass, Latitude, Longitude) %>%
  rename(Species = Species.y)

unique(df_lpi_mbs$GBIF_ID)

# df_mammal_cv_bs <- df_lpi_mbs %>%
#   mutate(log_cv = log10(cv_value),
#          # logit_p_cv = qlogis(p_cv_value),
#          logit_k_cv = qlogis(k_cv_value)) %>%
#   filter(cv_window == 10) %>%
#   group_by(Group, Binomial, GBIF_ID, cv_window) %>%
#   summarise(mean_cv = mean(cv_value, na.rm = T),
#             mean_log_cv = mean(log_cv, na.rm = T),
#             # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
#             mean_logit_k_cv = mean(logit_k_cv, na.rm = T),
#             mean_k_cv = mean(k_cv_value, na.rm = T),
#             mean_body_mass = mean(mean_body_mass, na.rm = T),
#             mean_latitude = mean(Latitude, na.rm = T)) %>%
#   mutate(log_mean_body_mass = log10(mean_body_mass),
#          abs_latitude = abs(mean_latitude))
# 
# unique(df_mammal_cv_bs$GBIF_ID)

df_lpi_fs <- merge(df_lpi_mbs, df_mammal_fs, by = "GBIF_ID")
unique(df_lpi_fs$ID)

##### Relationship between CV ~ Body Size + PC1 #####
df_mammal_cv_bs_fs <- df_lpi_fs %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value)) %>%
  rename(Binomial = Binomial.y) %>%
  group_by(Group, GBIF_ID, Binomial, cv_window) %>%
  summarise(mean_cv = mean(cv_value, na.rm = T),
            mean_log_cv = mean(log_cv, na.rm = T),
            # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
            logit_kcv = mean(logit_k_cv, na.rm = T),
            PC1 = mean(PC1, na.rm = T),
            PC2 = mean(PC2, na.rm = T),
            mean_k_cv = mean(k_cv_value, na.rm = T),
            mean_body_mass = mean(mean_body_mass, na.rm = T),
            mean_latitude = mean(Latitude, na.rm = T)) %>%
  mutate(log_mean_body_mass = log10(mean_body_mass),
         abs_latitude = abs(mean_latitude)) %>%
  ungroup()

df_mammal_cv_bs_fs_w <- df_mammal_cv_bs_fs %>% 
  filter(cv_window == 5) %>%
  mutate(sc_log_mean_body_mass = scale(log_mean_body_mass)[,1]) #dont need to scale based on global body size, only need that for priors in joint model

unique(df_mammal_cv_bs_fs_w$Binomial)


###### Model for mean differences between mammal groupings ######
priors_anova <- c(
  set_prior("normal(0, 3)", class = "b"),
  set_prior("normal(0, 3)", class = "Intercept"),
  set_prior("exponential(1)", class = "phi")
)

brms_k_cv_group <- brm(
  formula = mean_k_cv ~ Group,
  data = df_mammal_cv_bs_fs_w,
  prior = priors_anova,
  family = Beta(link = "logit"),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99)
)

loo_bs_cv_group <- brms::loo(brms_k_cv_group)
loo_bs_cv_group

pp_check(brms_k_cv_group, ndraws = 100)
bayes_R2(brms_k_cv_group, probs = c(0.055, 0.945))
summary(brms_k_cv_group, prob = 0.89)
plot(brms_k_cv_group)

#Look at mean differences
#Regrid to resonse scale
emm_resp <- emmeans(brms_k_cv_group, ~ Group, type = "response") %>%
  regrid(transform = "response")

#now pull draws from the regridded object
emm_draws <- gather_emmeans_draws(emm_resp)

median_hdci(emm_draws, .width = 0.89)

emm_draws$Group <- factor(emm_draws$Group, levels=c('Marine Mammals', 'Elephants', 'Ungulates', 'Carnivores', 'Marsupials', 'Primates', 'Glires', 'Bats', 'Insectivores'))

#Plot group posterior distributions
gg_cv_means <- ggplot(emm_draws, aes(x = .value, y = Group, fill = Group)) +
  stat_halfeye(alpha = 0.7, .width = c(0.89, 0.95)) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Estimated CV",
       y = "Group") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_cv_means

# ggsave("../R/Figures/Figure SX - Predicted CV Across Groups.jpeg", plot = gg_cv_means, width = 6, height = 6)

##### Model for CV ~ BS + Rmax #####

#Match Phylogeny up
sp  <- unique(df_mammal_cv_bs_fs_w$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_mammal_cv_bs_fs_w) <- df_mammal_cv_bs_fs_w$Binomial

#Phylogeny correlation matrix
cv_fs_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_mammal_cv_bs_fs_w$phylo <- factor(df_mammal_cv_bs_fs_w$Binomial, levels = rownames(cv_fs_corrma))

#Model Priors
priors_kcv_bs <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("exponential(1)", class = "phi")
  # set_prior("student_t(3, 0, 10)", class = "sigma")
)

priors_kcv_null <- c(
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("exponential(1)", class = "phi")
)

brms_kcv_null <- brm(
  formula = mean_k_cv ~ 1 + (1 | gr(phylo, cov = A)),
  data = df_mammal_cv_bs_fs_w,
  data2 = list(A = cv_fs_corrma),
  prior = priors_kcv_null,
  family = 
    # gaussian(),
    Beta(link = "logit"),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_kcv_null <- brms::loo(brms_kcv_null)
loo_kcv_null

brms_kcv_bs <- brm(
  formula = mean_k_cv ~ sc_log_mean_body_mass + (1 | gr(phylo, cov = A)),
  data = df_mammal_cv_bs_fs_w,
  data2 = list(A = cv_fs_corrma),
  prior = 
    priors_kcv_bs,
  family = 
    # gaussian(),
    Beta(link = "logit"),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_kcv_bs <- brms::loo(brms_kcv_bs)
loo_kcv_bs

brms_kcv_bs_pc1 <- brm(
  formula = mean_k_cv ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A)),
  data = df_mammal_cv_bs_fs_w,
  data2 = list(A = cv_fs_corrma),
  prior = 
    priors_kcv_bs,
  family = 
    # gaussian(),
    Beta(link = "logit"),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_kcv_bs_pc1 <- brms::loo(brms_kcv_bs_pc1)
loo_kcv_bs_pc1

r2_bayes(brms_kcv_null)
r2_bayes(brms_kcv_bs)
r2_bayes(brms_kcv_bs_pc1)

pp_check(brms_kcv_bs_pc1, ndraws = 100)
bayes_R2(brms_kcv_bs_pc1, probs = c(0.055, 0.945))
summary(brms_kcv_bs_pc1, prob = 0.89)
plot(brms_kcv_bs_pc1)


#Estimates of body size unscaled
draws <- as_draws_df(brms_kcv_bs_pc1) %>%
  mutate(
    beta_x  = b_sc_log_mean_body_mass / sd(df_mammal_cv_bs_fs_w$log_mean_body_mass),
    alpha_x = b_Intercept - b_sc_log_mean_body_mass * mean(df_mammal_cv_bs_fs_w$log_mean_body_mass) / sd(df_mammal_cv_bs_fs_w$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .055), u95 = quantile(beta_x, .945))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .055), u95 = quantile(alpha_x, .945))


#Prediction df
df_pred_bs_cv <- tibble(
  sc_log_mean_body_mass = seq(min(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass),
                              max(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass), length.out = 100),
  PC1 = mean(df_mammal_cv_bs_fs_w$PC1),
  log_mean_body_mass = seq(min(df_mammal_cv_bs_fs_w$log_mean_body_mass),
                              max(df_mammal_cv_bs_fs_w$log_mean_body_mass), length.out = 100)
)

#Posterior predictions
epred_draws_bs_cv <- add_epred_draws(brms_kcv_bs_pc1, newdata = df_pred_bs_cv, re_formula = NA)

gg_epred_bs_cv <- ggplot(epred_draws_bs_cv, aes(x = log_mean_body_mass, y = .epred)) +
  geom_point(data = df_mammal_cv_bs_fs_w, aes(x = log_mean_body_mass, y = mean_k_cv), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted CV (k)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_bs_cv

#Posterior predictions
linpred_draws_bs_cv <- add_linpred_draws(brms_kcv_bs_pc1, newdata = df_pred_bs_cv, re_formula = NA)

gg_linpred_bs_cv <- ggplot(linpred_draws_bs_cv, aes(x = log_mean_body_mass, y = .linpred)) +
  geom_point(data = df_mammal_cv_bs_fs_w, aes(x = log_mean_body_mass, y = qlogis(mean_k_cv)), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted CV (logit)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-4.1, 2.15), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_linpred_bs_cv

#Prediction df
df_pred_pc1_cv <- tibble(
  sc_log_mean_body_mass = mean(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass),
  PC1 = seq(min(df_mammal_cv_bs_fs_w$PC1),
               max(df_mammal_cv_bs_fs_w$PC1), length.out = 100)
)

#Get expected posterior predictions
epred_draws_pc1_cv <- add_epred_draws(brms_kcv_bs_pc1, newdata = df_pred_pc1_cv, re_formula = NA)

gg_epred_pc1_cv <- ggplot(epred_draws_pc1_cv, aes(x = PC1, y = .epred)) +
  geom_point(data = df_mammal_cv_bs_fs_w, aes(x = PC1, y = mean_k_cv), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "PC1", y = "Predicted CV (k)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_pc1_cv

#Get expected posterior predictions
linpred_draws_pc1_cv <- add_linpred_draws(brms_kcv_bs_pc1, newdata = df_pred_pc1_cv, re_formula = NA)

gg_linpred_pc1_cv <- ggplot(linpred_draws_pc1_cv, aes(x = PC1, y = .linpred)) +
  geom_point(data = df_mammal_cv_bs_fs_w, aes(x = PC1, y = qlogis(mean_k_cv)), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "PC1", y = "Predicted CV (logit)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-4.1, 2.15), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_linpred_pc1_cv

gg_epred_bs_pc1_cv <- plot_grid(gg_epred_bs_cv, gg_epred_pc1_cv, nrow = 1, align = 'hv')
gg_epred_bs_pc1_cv

gg_linpred_bs_pc1_cv <- plot_grid(gg_linpred_bs_cv, gg_linpred_pc1_cv, nrow = 1, align = 'hv', labels = c("c", "d"))

gg_linpred_bs_pc1_cv

gg_epred_fs_cv <- plot_grid(gg_epred_bs_pc1, gg_linpred_bs_pc1_cv, nrow = 2, align = "hv")
gg_epred_fs_cv

# ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 2 - Predicted Rmax & CV ~ BS + PC1.jpeg", plot = gg_epred_fs_cv, width = 8, height = 8)


##### Joint model between CV ~ rmax & rmax ~ bs + PC1 #####
df_rmax_fs_lpi <- merge(df_rmax_fs, df_lpi_mbs, by = c("GBIF_ID")) %>%
  rename(Binomial = Binomial.x)

df_rmax_fs_lpi_sub <- df_rmax_fs_lpi %>%
  group_by(Binomial, GBIF_ID, cv_window) %>%
  summarise(mean_k_cv = mean(k_cv_value),
            rmax = mean(rmax), 
            Mass_g = mean(Mass_g),
            PC1 = mean(PC1),
            PC2 = mean(PC2)) %>%
  ungroup() %>%
  mutate(log_rmax = log10(rmax),
         log_mean_body_mass = log10(Mass_g),
         logit_kcv = qlogis(mean_k_cv))

unique(df_rmax_fs_lpi_sub$Binomial)

df_rmax_fs_lpi_sub_w <- df_rmax_fs_lpi_sub %>%
  filter(cv_window == 5) %>%
  mutate(sc_log_mean_body_mass = (log_mean_body_mass - global_mean_log_body_mass)/global_sd_log_body_mass)

#Match Phylogeny up
sp  <- unique(df_rmax_fs_lpi_sub_w$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_lpi_sub_w) <- df_rmax_fs_lpi_sub_w$Binomial

#Phylogeny correlation matrix
rmax_fs_lpis_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_lpi_sub_w$phylo <- factor(df_rmax_fs_lpi_sub_w$Binomial, levels = rownames(rmax_fs_lpis_corrma))

post_fixef_bspc
post_sd_phy_bspc
post_sigma_bspc

priors_updated <- c(
  
  #Priors for rmax if using count method data only, no prior information, weakly informative priors
  # set_prior("normal(0, 1)", class = "b", resp = "logrmax"),
  # set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo", resp = "logrmax"),
  # set_prior("student_t(3, 0, 10)",    class = "sigma", resp = "logrmax"),
  
  #Priors for rmax model updating priors with information from rmax previous model that uses all rmax data without species that have temporal series data
  set_prior("normal(-0.32941983, 0.10)", class = "b", coef = "sc_log_mean_body_mass", resp = "logrmax"),
  set_prior("normal(0.16893643, 0.10)", class = "b", coef = "PC1", resp = "logrmax"),
  set_prior("normal(0, 1)", class = "Intercept", resp = "logrmax"),
  set_prior("student_t(3, 0, 0.6)", class = "sd", group = "phylo", resp = "logrmax"),
  set_prior("student_t(3, 0, 0.6)",    class = "sigma", resp = "logrmax"),
  
  #Priors for k_cv model, no prior information, weakly informative priors
  set_prior("normal(0, 1)", class = "b", resp = "meankcv"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo", resp = "meankcv"),
  set_prior("exponential(1)", class = "phi", resp = "meankcv")
  # set_prior("student_t(3, 0, 10)",    class = "sigma", resp = "logitkcv")
)


fit_joint_rmax_pc_bs <- brm(
  brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    set_rescor(F),
  data = df_rmax_fs_lpi_sub_w,
  data2 = list(A = rmax_fs_lpis_corrma),
  prior = priors_updated,
  family = list(gaussian(), 
                # gaussian()
                Beta(link = "logit")
  ),  
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)


fit_joint_rmax_rpc_bs <- brm(
    brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ log_rmax + sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
      set_rescor(F),
    data = df_rmax_fs_lpi_sub_w,
    data2 = list(A = rmax_fs_lpis_corrma),
    prior = priors_updated,
    family = list(gaussian(), 
                  # gaussian()
                  Beta(link = "logit")
    ),  
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

fit_joint_rmax_rpc <- brm(
    brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ log_rmax + PC1 + (1 | gr(phylo, cov = A))) +
      set_rescor(F),
    data = df_rmax_fs_lpi_sub_w,
    data2 = list(A = rmax_fs_lpis_corrma),
    prior = priors_updated,
    family = list(gaussian(), 
                  # gaussian()
                  Beta(link = "logit")
    ),  
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

fit_joint_rmax_rbs <- brm(
  brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ log_rmax + sc_log_mean_body_mass + (1 | gr(phylo, cov = A))) +
    set_rescor(F),
  data = df_rmax_fs_lpi_sub_w,
  data2 = list(A = rmax_fs_lpis_corrma),
  prior = priors_updated,
  family = list(gaussian(), 
                # gaussian()
                Beta(link = "logit")
  ),  
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

fit_joint_rmax <- brm(
  brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ log_rmax + (1 | gr(phylo, cov = A))) +
    set_rescor(F),
  data = df_rmax_fs_lpi_sub_w,
  data2 = list(A = rmax_fs_lpis_corrma),
  prior = priors_updated,
  family = list(gaussian(), 
                # gaussian()
                Beta(link = "logit")
  ),  
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_bs_rmax_pc_bs <- brms::loo(fit_joint_rmax_pc_bs)
loo_bs_rmax_pc_bs

loo_bs_rmax_rpc_bs <- brms::loo(fit_joint_rmax_rpc_bs)
loo_bs_rmax_rpc_bs

loo_bs_rmax_rpc <- brms::loo(fit_joint_rmax_rpc)
loo_bs_rmax_rpc

loo_bs_rmax_rbs <- brms::loo(fit_joint_rmax_rbs)
loo_bs_rmax_rbs

loo_bs_rmax <- brms::loo(fit_joint_rmax)
loo_bs_rmax

loo_compare(loo_bs_rmax_pc_bs, loo_bs_rmax_rpc_bs, loo_bs_rmax, loo_bs_rmax_rpc, loo_bs_rmax_rbs)
loo_model_weights(list(loo_bs_rmax_pc_bs, loo_bs_rmax_rpc_bs, loo_bs_rmax_rpc, loo_bs_rmax, loo_bs_rmax_rbs))

pp_check(fit_joint_rmax_pc_bs, ndraws = 100, resp = "logrmax")
pp_check(fit_joint_rmax_rpc_bs, ndraws = 100, resp = "logrmax")
pp_check(fit_joint_rmax_rpc, ndraws = 100, resp = "meankcv")
pp_check(fit_joint_rmax_rbs, ndraws = 100, resp = "meankcv")
pp_check(fit_joint_rmax, ndraws = 100, resp = "meankcv")

bayes_R2(fit_joint_rmax, resp = "logrmax", probs = c(0.055, 0.945))
bayes_R2(fit_joint_rmax, resp = "meankcv", probs = c(0.055, 0.945))
summary(fit_joint_rmax, prob = 0.89)
plot(fit_joint_rmax)

#Estimates of body size unscaled
draws <- as_draws_df(fit_joint_rmax) %>%
  mutate(
    beta_x  = b_logrmax_sc_log_mean_body_mass / sd(df_rmax_fs_lpi_sub_w$log_mean_body_mass),
    alpha_x = b_logrmax_Intercept - b_logrmax_sc_log_mean_body_mass * mean(df_rmax_fs_lpi_sub_w$log_mean_body_mass) / sd(df_rmax_fs_lpi_sub_w$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .055), u95 = quantile(beta_x, .945))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .055), u95 = quantile(alpha_x, .945))


mcmc_plot(fit_joint_rmax, 
          type = "intervals", 
          prob = 0.5, 
          prob_outer = 0.89,
          point_est = "median")

conditional_effects(fit_joint_rmax, "sc_log_mean_body_mass", resp = "logrmax")
conditional_effects(fit_joint_rmax, "PC1", resp = "logrmax")
conditional_effects(fit_joint_rmax, "log_rmax", resp = "meankcv")
# conditional_effects(fit_joint_rmax, "sc_log_mean_body_mass", resp = "meankcv")
# conditional_effects(fit_joint_rmax, "PC1", resp = "meankcv")

posterior <- as_draws_df(fit_joint_rmax)

posterior_ind <- posterior %>%
  mutate(
    indirect_bs = b_logrmax_sc_log_mean_body_mass * b_meankcv_log_rmax,
    indirect_pc1 = b_logrmax_PC1 * b_meankcv_log_rmax)

rmax_median_ci <- posterior_ind %>%
  summarise(
    lower = quantile(b_meankcv_log_rmax, 0.055),
    median = mean(b_meankcv_log_rmax),
    upper = quantile(b_meankcv_log_rmax, 0.945),
    pred = "rmax")

bs_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_bs, 0.055),
    median = mean(indirect_bs),
    upper = quantile(indirect_bs, 0.945),
    pred = "body size")

pc1_median_ci<- posterior_ind %>%
  summarise(
    lower = quantile(indirect_pc1, 0.055),
    median = mean(indirect_pc1),
    upper = quantile(indirect_pc1, 0.945),
    pred = "fast-slow continuum")

df_median_cis <- rbind(rmax_median_ci, bs_median_ci, pc1_median_ci)

gg_median_cis_rmaxbspc <- ggplot(df_median_cis, aes(x = median, y = pred)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_point(shape = 21, stroke = 0.5, fill = "aliceblue", color = "black", size = 3) +
  scale_x_continuous(limits = c(-1.3, 1.3), breaks = scales::pretty_breaks()) +
  labs(x = "Direct or Indirect Effect on CV", 
       y = "Covariate") +
  theme_classic(base_size = 14)

gg_median_cis_rmaxbspc

# ggsave("../R/Figures/Figure 3 - Bayesian CV - Rmax, BS, PC1 - Median CIs.jpeg", plot = gg_median_cis_rmaxbspc, width = 5, height = 4)

df_pred_pc1 <- tibble(
  log_mean_body_mass = mean(df_rmax_fs_lpi_sub_w$sc_log_mean_body_mass, na.rm = TRUE),
  PC1 = mean(df_rmax_fs_lpi_sub_w$PC1, na.rm = TRUE),
  log_rmax = seq(-2,
                 max(df_rmax_fs_lpi_sub_w$log_rmax, na.rm = TRUE), length.out = 100),
)

epred_pc1_cv <- add_epred_draws(fit_joint_rmax,
                                newdata = df_pred_pc1,
                                resp = "meankcv",
                                re_formula = NA
                                )

linpred_pc1_cv <- add_linpred_draws(fit_joint_rmax,
                                newdata = df_pred_pc1,
                                resp = "meankcv",
                                re_formula = NA
                                )

gg_rmax_cv <- ggplot(linpred_pc1_cv, aes(x = log_rmax, y = .linpred)) +
  geom_point(data = as.data.frame(df_rmax_fs_lpi_sub_w), aes(x = log_rmax, y = logit_kcv),
             inherit.aes = FALSE, #shape = 21, stroke = 0.5,
             size = 2.5, color = "black",alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "Growth Potential (log(rmax))", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-4.0, 1.8), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_rmax_cv

# ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 3 - Bayesian CV - Rmax.jpeg", plot = gg_rmax_cv, width = 5, height = 5)


##### Supplement Models ######
###### Relationship between rmax ~ Body Size + PC1 - Count Method Data Only  ######
#Line up Species names with Phylogeny for rmax, PC1 data frame
df_rmax_fs <- merge(df_mammal_fs, df_growth_count, by = c("GBIF_ID")) %>%
  dplyr::select(- X.x, -X.1, -X.y, - Species.y, - Binomial, Binomial = Species.x)

df_rmax_fs_sub <- df_rmax_fs %>%
  # anti_join(df_lpi_sp, by = "GBIF_ID") %>%
  group_by(Binomial, GBIF_ID) %>%
  summarise(rmax = mean(Growth_rate_per_yr), 
            Mass_g = mean(Mass_g),
            PC1 = mean(PC1),
            PC2 = mean(PC2)) %>%
  ungroup() %>%
  mutate(log_rmax = log10(rmax),
         log_mean_body_mass = log10(Mass_g),
         sc_log_rmax = scale(log_rmax)[,1],
         sc_log_mean_body_mass = scale(log_mean_body_mass)[,1])

unique(df_rmax_fs_sub$Binomial)

#Match Phylogeny up
sp  <- unique(df_rmax_fs_sub$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_sub) <- df_rmax_fs_sub$Binomial

#Phylogeny correlation matrix
rmax_fs_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_sub$phylo <- factor(df_rmax_fs_sub$Binomial, levels = rownames(rmax_fs_corrma))

#Model Priors
priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("student_t(3, 0, 10)", class = "sigma")
)

brms_bs_rmax <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_sub,
  data2 = list(A = rmax_fs_corrma),
  prior = priors,
  family = gaussian(),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_bs_rmax <- brms::loo(brms_bs_rmax)
loo_bs_rmax

brms_bs_rmax
pp_check(brms_bs_rmax, ndraws = 100)
bayes_R2(brms_bs_rmax)
summary(brms_bs_rmax)
plot(brms_bs_rmax)

#Estimates of body size unscaled
draws <- as_draws_df(brms_bs_rmax) %>%
  mutate(
    beta_x  = b_sc_log_mean_body_mass / sd(df_rmax_fs_sub$log_mean_body_mass),
    alpha_x = b_Intercept - b_sc_log_mean_body_mass * mean(df_rmax_fs_sub$log_mean_body_mass) / sd(df_rmax_fs_sub$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .05), u95 = quantile(beta_x, .95))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .05), u95 = quantile(alpha_x, .95))


#Create new data for prediction
df_pred_bs <- tibble(
  sc_log_mean_body_mass = seq(min(df_rmax_fs_sub$sc_log_mean_body_mass),
                              max(df_rmax_fs_sub$sc_log_mean_body_mass), length.out = 100),
  PC1 = mean(df_rmax_fs_sub$PC1),
  log_mean_body_mass = seq(min(df_rmax_fs_sub$log_mean_body_mass),
                           max(df_rmax_fs_sub$log_mean_body_mass), length.out = 100))


#Get expected posterior predictions
epred_draws_bs <- add_epred_draws(brms_bs_rmax, newdata = df_pred_bs, re_formula = NA)

gg_epred_bs <- ggplot(epred_draws_bs, aes(x = log_mean_body_mass, y = .epred)) +
  # geom_point(data = df_rmax_fs_sub, aes(x = log_mean_body_mass, y = log_rmax),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(aes(fill = after_stat(.width)), .width = 0.95, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted log(rmax)") +
  # scale_fill_gradient(low = "#3CA373", high = "#EB8F00") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2.2, 1.4), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_bs

#Create new data for prediction
df_pred_pc1 <- tibble(
  sc_log_mean_body_mass = mean(df_rmax_fs_sub$sc_log_mean_body_mass),
  PC1 = seq(min(df_rmax_fs_sub$PC1),
            max(df_rmax_fs_sub$PC1), length.out = 100)
)

#Get expected posterior predictions
epred_draws_pc1 <- add_epred_draws(brms_bs_rmax, newdata = df_pred_pc1, re_formula = NA)

gg_epred_pc1 <- ggplot(epred_draws_pc1, aes(x = PC1, y = .epred)) +
  # geom_point(data = df_rmax_fs_sub, aes(x = PC1, y = log_rmax),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(aes(fill = after_stat(.width)), .width = 0.95, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "PC1", y = "Predicted log(rmax)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-2.25, 1.4), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_pc1

gg_epred_bs_pc1 <- plot_grid(gg_epred_bs, gg_epred_pc1, nrow = 1, align = 'hv')
gg_epred_bs_pc1

ggsave("../R/Figures/Figure SX - Bayesian Rmax - BS + PC1 - Count Method.jpeg", plot = gg_epred_bs_pc1, width = 12, height = 4)



###### Relationship between rmax ~ Body Size + PC1 + MIMR ######
df_rmax_fs_metab <- merge(df_metab, df_rmax_fs, by = c("GBIF_ID")) %>%
  # anti_join(df_lpi_sp, by = "GBIF_ID") %>%
  group_by(Binomial, GBIF_ID) %>%
  summarise(rmax = mean(rmax), 
            Mass_g = mean(Mass_g.y),
            PC1 = mean(PC1),
            PC2 = mean(PC2),
            Metabolism_W = mean(Metabolism_W)) %>%
  ungroup() %>%
  mutate(log_rmax = log10(rmax),
         log_mean_body_mass = log10(Mass_g),
         sc_log_rmax = scale(log_rmax)[,1],
         sc_log_mean_body_mass = scale(log_mean_body_mass)[,1])



###### Adding Mass Independent Metabolic Rate ######
fit_log_smr <- lm(log10(Metabolism_W) ~ log10(Mass_g),
                  data = df_rmax_fs_metab)

df_rmax_fs_metab$MIMR <- residuals(fit_log_smr)
df_rmax_fs_metab$sc_MIMR <- scale(df_rmax_fs_metab$MIMR)[,1]


unique(df_rmax_fs_metab$Binomial)

summary(lm(data = df_rmax_fs_metab, MIMR ~ PC1))
plot(df_rmax_fs_metab$PC1, df_rmax_fs_metab$MIMR)
abline(lm(data = df_rmax_fs_metab, MIMR ~ PC1))

#Match Phylogeny up
sp  <- unique(df_rmax_fs_metab$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_metab) <- df_rmax_fs_metab$Binomial

#Phylogeny correlation matrix
rmax_fs_metab_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_metab$phylo <- factor(df_rmax_fs_metab$Binomial, levels = rownames(rmax_fs_metab_corrma))

#Model Priors
priors_1 <- c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("student_t(3, 0, 10)", class = "sigma")
)


brms_bs_rmax_metab <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + sc_MIMR + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_metab,
  data2 = list(A = rmax_fs_metab_corrma),
  prior = priors_1,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

brms_bs_rmax_pc <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_metab,
  data2 = list(A = rmax_fs_metab_corrma),
  prior = priors_1,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)

brms_bs_rmax_metab_pc <- brm(
  formula = log_rmax ~ sc_log_mean_body_mass + PC1 + sc_MIMR + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_metab,
  data2 = list(A = rmax_fs_metab_corrma),
  prior = priors_1,
  family = gaussian(),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

loo_bs_rmax_pc <- brms::loo(brms_bs_rmax_pc)
loo_bs_rmax_pc

loo_bs_rmax_metab <- brms::loo(brms_bs_rmax_metab)
loo_bs_rmax_metab

loo_bs_rmax_metab_pc <- brms::loo(brms_bs_rmax_metab_pc)
loo_bs_rmax_metab_pc

loo_compare(loo_bs_rmax_pc, loo_bs_rmax_metab, loo_bs_rmax_metab_pc)

pp_check(brms_bs_rmax_metab_pc, ndraws = 100)
bayes_R2(brms_bs_rmax_metab_pc)
summary(brms_bs_rmax_metab_pc, prob = 0.89)
posterior_summary(brms_bs_rmax_metab_pc, probs = c(0.055, 0.945))
plot(brms_bs_rmax_metab_pc)

#Create dataframe for prediction
df_pred_bs <- tibble(
  sc_log_mean_body_mass = seq(min(df_rmax_fs_metab$sc_log_mean_body_mass),
                              max(df_rmax_fs_metab$sc_log_mean_body_mass), length.out = 100),
  PC1 = mean(df_rmax_fs_metab$PC1),
  sc_MIMR = mean(df_rmax_fs_metab$sc_MIMR),
)

#Get expected posterior predictions
epred_draws_bs <- add_epred_draws(brms_bs_rmax_metab_pc, newdata = df_pred_bs, re_formula = NA)

gg_epred_bs <- ggplot(epred_draws_bs, aes(x = sc_log_mean_body_mass, y = .epred)) +
  stat_lineribbon(.width = 0.89, alpha = 0.5, linewidth = 2, color = "black", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted log(rmax)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_epred_bs

#Create dataframe for prediction
df_pred_pc1 <- tibble(
  sc_log_mean_body_mass = mean(df_rmax_fs_metab$sc_log_mean_body_mass),
  PC1 = seq(min(df_rmax_fs_metab$PC1),
               max(df_rmax_fs_metab$PC1), length.out = 100),
  sc_MIMR = mean(df_rmax_fs_metab$sc_MIMR),
)

#Get expected posterior predictions
epred_draws_pc1 <- add_epred_draws(brms_bs_rmax_metab_pc, newdata = df_pred_pc1, re_formula = NA)

gg_epred_pc1 <- ggplot(epred_draws_pc1, aes(x = PC1, y = .epred)) +
  stat_lineribbon(.width = 0.89, alpha = 0.5, linewidth = 1.5, color = "black", fill = "grey") +
  labs(x = "PC1", y = "Predicted log(rmax)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_epred_pc1

#Create dataframe for prediction
df_pred_met <- tibble(
  sc_log_mean_body_mass = mean(df_rmax_fs_metab$sc_log_mean_body_mass),
  PC1 = mean(df_rmax_fs_metab$PC1),
  sc_MIMR = seq(min(df_rmax_fs_metab$sc_MIMR),
                 max(df_rmax_fs_metab$sc_MIMR), length.out = 100)
)

#Get expected posterior predictions
epred_draws_met <- add_epred_draws(brms_bs_rmax_metab, newdata = df_pred_met, re_formula = NA)

gg_epred_met <- ggplot(epred_draws_met, aes(x = sc_MIMR, y = .epred)) +
  stat_lineribbon(.width = 0.89, alpha = 0.5, linewidth = 1.5, color = "black", fill = "grey") +
  labs(x = "Mass Independent Specific Metabolic Rate", y = "Predicted log(rmax)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

gg_epred_met

gg_epred_bs_pc1_met <- plot_grid(gg_epred_bs, gg_epred_pc1, gg_epred_met, nrow = 1, align = 'hv')
gg_epred_bs_pc1_met

ggsave("../R/Figures/Figure SX - Bayesian Rmax - Metab.jpeg", plot = gg_epred_bs_pc1_met, width = 12, height = 4)



###### Relationship between CV ~ bs + PC1 + MIMR ######
df_rmax_fs_cv_metab <- merge(df_metab, df_lpi_fs, by = c("GBIF_ID")) %>%
  # anti_join(df_lpi_sp, by = "GBIF_ID") %>%
  mutate(log_cv = log10(cv_value),
         # logit_p_cv = qlogis(p_cv_value),
         logit_k_cv = qlogis(k_cv_value)) %>%
  group_by(Binomial.y, GBIF_ID) %>%
  summarise(Mass_g = mean(Mass_g),
            PC1 = mean(PC1),
            PC2 = mean(PC2),
            Metabolism_W = mean(Metabolism_W),
            mean_cv = mean(cv_value, na.rm = T),
            mean_log_cv = mean(log_cv, na.rm = T),
            # mean_logit_p_cv = mean(logit_p_cv, na.rm = T),
            logit_kcv = mean(logit_k_cv, na.rm = T),
            PC1 = mean(PC1, na.rm = T),
            PC2 = mean(PC2, na.rm = T),
            mean_k_cv = mean(k_cv_value, na.rm = T),
            mean_body_mass = mean(mean_body_mass, na.rm = T),
            mean_latitude = mean(Latitude, na.rm = T)) %>%
  ungroup() %>%
  mutate(log_mean_body_mass = log10(Mass_g),
         sc_log_mean_body_mass = scale(log_mean_body_mass)[,1]) %>%
  rename(Binomial = Binomial.y)

###### Adding Mass Independent Metabolic Rate ######
fit_log_smr <- lm(log10(Metabolism_W) ~ log10(Mass_g),
                  data = df_rmax_fs_cv_metab)

df_rmax_fs_cv_metab$MIMR <- residuals(fit_log_smr)
df_rmax_fs_cv_metab$sc_MIMR <- scale(df_rmax_fs_cv_metab$MIMR)[,1]

unique(df_rmax_fs_cv_metab$Binomial)

summary(lm(data = df_rmax_fs_cv_metab, MIMR ~ PC1))
plot(df_rmax_fs_cv_metab$PC1, df_rmax_fs_cv_metab$MIMR)
abline(lm(data = df_rmax_fs_cv_metab, MIMR ~ PC1))

unique(df_rmax_fs_cv_metab$Binomial.y)


#Match Phylogeny up
sp  <- unique(df_rmax_fs_cv_metab$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_cv_metab) <- df_rmax_fs_cv_metab$Binomial

#Phylogeny correlation matrix
cv_fs_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_cv_metab$phylo <- factor(df_rmax_fs_cv_metab$Binomial, levels = rownames(cv_fs_corrma))

#Model Priors
priors_kcv_bs <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0, 1)", class = "Intercept"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo"),
  set_prior("exponential(1)", class = "phi")
  # set_prior("student_t(3, 0, 10)", class = "sigma")
)

brms_k_cv_bs_met <- brm(
  formula = mean_k_cv ~ sc_log_mean_body_mass + PC1 + sc_MIMR + (1 | gr(phylo, cov = A)),
  data = df_rmax_fs_cv_metab,
  data2 = list(A = cv_fs_corrma),
  prior = priors_kcv_bs,
  family = 
    # gaussian(),
    Beta(link = "logit"),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.999),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_bs_metab <- brms::loo(brms_k_cv_bs_met)
loo_bs_metab

pp_check(brms_k_cv_bs_met, ndraws = 100)
bayes_R2(brms_k_cv_bs_met, probs = c(0.055, 0.945))
summary(brms_k_cv_bs_met, prob = 0.89)
plot(brms_k_cv_bs_met)


#Estimates of body size unscaled
draws <- as_draws_df(brms_k_cv_bs) %>%
  mutate(
    beta_x  = b_sc_log_mean_body_mass / sd(df_mammal_cv_bs_fs_w$log_mean_body_mass),
    alpha_x = b_Intercept - b_sc_log_mean_body_mass * mean(df_mammal_cv_bs_fs_w$log_mean_body_mass) / sd(df_mammal_cv_bs_fs_w$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .055), u95 = quantile(beta_x, .945))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .055), u95 = quantile(alpha_x, .945))


#Prediction df
df_pred_bs_cv <- tibble(
  sc_log_mean_body_mass = seq(min(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass),
                              max(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass), length.out = 100),
  PC1 = mean(df_mammal_cv_bs_fs_w$PC1),
  log_mean_body_mass = seq(min(df_mammal_cv_bs_fs_w$log_mean_body_mass),
                           max(df_mammal_cv_bs_fs_w$log_mean_body_mass), length.out = 100)
)

#Posterior predictions
epred_draws_bs_cv <- add_epred_draws(brms_k_cv_bs, newdata = df_pred_bs_cv, re_formula = NA)

gg_epred_bs_cv <- ggplot(epred_draws_bs_cv, aes(x = log_mean_body_mass, y = .epred)) +
  # geom_point(data = df_mammal_cv_bs_fs_w, aes(x = log_mean_body_mass, y = mean_k_cv),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(aes(fill = after_stat(.width)), .width = 0.95, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted CV (k)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_bs_cv

#Posterior predictions
linpred_draws_bs_cv <- add_linpred_draws(brms_k_cv_bs, newdata = df_pred_bs_cv, re_formula = NA)

gg_linpred_bs_cv <- ggplot(linpred_draws_bs_cv, aes(x = log_mean_body_mass, y = .linpred)) +
  # geom_point(data = df_mammal_cv_bs_fs_w, aes(x = log_mean_body_mass, y = mean_k_cv),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(aes(fill = after_stat(.width)), .width = 0.89, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "log(Body Mass)", y = "Predicted CV (logit)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-3.75, 2.15), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_linpred_bs_cv

#Prediction df
df_pred_pc1_cv <- tibble(
  sc_log_mean_body_mass = mean(df_mammal_cv_bs_fs_w$sc_log_mean_body_mass),
  PC1 = seq(min(df_mammal_cv_bs_fs_w$PC1),
            max(df_mammal_cv_bs_fs_w$PC1), length.out = 100)
)

#Get expected posterior predictions
epred_draws_pc1_cv <- add_epred_draws(brms_k_cv_bs, newdata = df_pred_pc1_cv, re_formula = NA)

gg_epred_pc1_cv <- ggplot(epred_draws_pc1_cv, aes(x = PC1, y = .epred)) +
  # geom_point(data = df_mammal_cv_bs_fs_w, aes(x = PC1, y = mean_k_cv),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "PC1", y = "Predicted CV (k)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_epred_pc1_cv

#Get expected posterior predictions
linpred_draws_pc1_cv <- add_linpred_draws(brms_k_cv_bs, newdata = df_pred_pc1_cv, re_formula = NA)

gg_linpred_pc1_cv <- ggplot(linpred_draws_pc1_cv, aes(x = PC1, y = .linpred)) +
  # geom_point(data = df_mammal_cv_bs_fs_w, aes(x = PC1, y = mean_k_cv),
  #            shape = 21, size = 2.5, stroke = 0.5, color = "black", fill = "aliceblue") +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0, color = "black", fill = "grey") +
  labs(x = "PC1", y = "Predicted CV (logit)") +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(-3.75, 2.15), breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_linpred_pc1_cv


gg_epred_bs_pc1_cv <- plot_grid(gg_epred_bs_cv, gg_epred_pc1_cv, nrow = 1, align = 'hv')
gg_epred_bs_pc1_cv

gg_linpred_bs_pc1_cv <- plot_grid(gg_linpred_bs_cv, gg_linpred_pc1_cv, nrow = 1, align = 'hv')
gg_linpred_bs_pc1_cv

gg_epred_fs_cv <- plot_grid(gg_epred_bs_pc1, gg_linpred_bs_pc1_cv, nrow = 2, align = "hv")
gg_epred_fs_cv

# ggsave("../R/Figures/Figure 2 - Predicted Rmax & CV ~ BS + PC1.jpeg", plot = gg_epred_fs_cv, width = 8, height = 8)


###### Joint model between CV ~ rmax + MISMR & rmax ~ bs + PC1 ######
df_rmax_fs_lpi_met <- merge(df_rmax_fs_lpi, df_metab, by = c("GBIF_ID")) %>%
  filter(cv_window == 5) %>%
  group_by(Binomial, GBIF_ID) %>%
  summarise(rmax = mean(rmax), 
            Mass_g = mean(Mass_g.y),
            mean_k_cv = mean(k_cv_value),
            PC1 = mean(PC1),
            PC2 = mean(PC2),
            Mass_Spec_Meta_Watt = mean(Mass_Spec_Meta_Watt)) %>%
  ungroup() %>%
  mutate(log_rmax = log10(rmax),
         log_mean_body_mass = log10(Mass_g),
         logit_kcv = qlogis(mean_k_cv))

fit_log_smr <- lm(log10(Mass_Spec_Meta_Watt) ~ log10(Mass_g),
                  data = df_rmax_fs_lpi_met)

df_rmax_fs_lpi_met$MIMR <- residuals(fit_log_smr)
df_rmax_fs_lpi_met$sc_MIMR <- scale(df_rmax_fs_lpi_met$MIMR)[,1]

unique(df_rmax_fs_lpi_met$Binomial)

df_rmax_fs_lpi_sub_w <- df_rmax_fs_lpi_met %>%
  mutate(sc_log_mean_body_mass = scale(log_mean_body_mass)[,1])

#Match Phylogeny up
sp  <- unique(df_rmax_fs_lpi_sub_w$Binomial) 

#Prune tree
setdiff(sp, df_tree$tip.label)
keep <- intersect(df_tree$tip.label, sp)
tr_pruned <- drop.tip(df_tree, setdiff(df_tree$tip.label, keep))

rownames(df_rmax_fs_lpi_sub_w) <- df_rmax_fs_lpi_sub_w$Binomial

#Phylogeny correlation matrix
rmax_fs_lpis_corrma <- vcv(tr_pruned, corr = TRUE)

#Make grouping factor phylo
df_rmax_fs_lpi_sub_w$phylo <- factor(df_rmax_fs_lpi_sub_w$Binomial, levels = rownames(rmax_fs_lpis_corrma))

post_fixef_bspc
post_sd_phy_bspc
post_sigma_bspc

priors_updated <- c(

  # #Priors for rmax model from previous
  # set_prior("normal(-0.32941983, 0.10)", class = "b", coef = "sc_log_mean_body_mass", resp = "logrmax"),
  # set_prior("normal(0.16893643, 0.10)", class = "b", coef = "PC1", resp = "logrmax"),
  # set_prior("normal(0, 1)", class = "Intercept", resp = "logrmax"),
  # set_prior("student_t(3, 0, 0.6)", class = "sd", group = "phylo", resp = "logrmax"),
  # set_prior("student_t(3, 0, 0.6)",    class = "sigma", resp = "logrmax"),
  
  set_prior("normal(0, 1)", class = "b", resp = "logrmax"),
  set_prior("normal(0, 1)", class = "Intercept", resp = "logrmax"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo", resp = "logrmax"),
  set_prior("student_t(3, 0, 10)",    class = "sigma", resp = "logrmax"),

  #Priors for k_cv model, no prior information, weakly informative priors
  set_prior("normal(0, 1)", class = "b", resp = "meankcv"),
  set_prior("student_t(3, 0, 10)", class = "sd", group = "phylo", resp = "meankcv"),
  set_prior("exponential(1)", class = "phi", resp = "meankcv")
  # set_prior("student_t(3, 0, 0.6)",    class = "sigma", resp = "logitkcv")
)

# 
# fit_joint_rmax_rpc <- brm(
#     brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
#     brms::bf(logit_kcv ~ log_rmax + PC1 + (1 | gr(phylo, cov = A))) +
#     set_rescor(TRUE),
#   data = df_rmax_fs_lpi_sub_w,
#   data2 = list(A = rmax_fs_lpis_corrma),
#   prior = priors_updated,
#   family = list(gaussian(), gaussian()
#   ),
#   chains = 4, cores = 4,
#   warmup = 2000, iter = 4000,
#   control = list(adapt_delta = 0.95),
#   sample_prior = "yes",
#   save_pars = save_pars(all = TRUE)
# )
# 
# fit_joint_rmax_rbs <- brm(
#   brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
#     brms::bf(logit_kcv ~ log_rmax + log_mean_body_mass + (1 | gr(phylo, cov = A))) +
#     set_rescor(TRUE),
#   data = df_rmax_fs_lpi_sub_w,
#   data2 = list(A = rmax_fs_lpis_corrma),
#   prior = priors_updated,
#   family = list(gaussian(), gaussian()
#   ),  
#   chains = 4, cores = 4, 
#   warmup = 2000, iter = 4000,
#   control = list(adapt_delta = 0.95),
#   sample_prior = "yes",
#   save_pars = save_pars(all = TRUE)
# )

fit_joint_rmax <- brm(
  brms::bf(log_rmax ~ sc_log_mean_body_mass + PC1 + (1 | gr(phylo, cov = A))) +
    brms::bf(mean_k_cv ~ log_rmax + sc_MIMR + (1 | gr(phylo, cov = A))) +
    set_rescor(F),
  data = df_rmax_fs_lpi_sub_w,
  data2 = list(A = rmax_fs_lpis_corrma),
  prior = priors_updated,
  family = #Beta(link = "logit"),
    list(
    gaussian(),
                # gaussian(),
                Beta(link = "logit")
  ),
  chains = 4, cores = 4, 
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.99),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

# loo_bs_rmax_rpc <- brms::loo(fit_joint_rmax_rpc)
# 
# loo_bs_rmax_rbs <- brms::loo(fit_joint_rmax_rbs)

loo_bs_rmax <- brms::loo(fit_joint_rmax)

# loo_compare(loo_bs_rmax_rpc, loo_bs_rmax_rbs, loo_bs_rmax)

loo_bs_rmax

pp_check(fit_joint_rmax, ndraws = 100, resp = "logrmax")
pp_check(fit_joint_rmax, ndraws = 100, resp = "meankcv")
bayes_R2(fit_joint_rmax, resp = "logrmax")
bayes_R2(fit_joint_rmax, resp = "meankcv")
summary(fit_joint_rmax, prob = 0.89)
plot(fit_joint_rmax)

#Estimates of body size unscaled
draws <- as_draws_df(fit_joint_rmax) %>%
  mutate(
    beta_x  = b_logrmax_sc_log_mean_body_mass / sd(df_rmax_fs_lpi_sub_w$log_mean_body_mass),
    alpha_x = b_logrmax_Intercept - b_logrmax_sc_log_mean_body_mass * mean(df_rmax_fs_lpi_sub_w$log_mean_body_mass) / sd(df_rmax_fs_lpi_sub_w$log_mean_body_mass)
  )

draws %>% 
  summarise(Coefficient = mean(beta_x), l95 = quantile(beta_x, .05), u95 = quantile(beta_x, .95))

draws %>% 
  summarise(Intercept = mean(alpha_x), l95 = quantile(alpha_x, .05), u95 = quantile(alpha_x, .95))


mcmc_plot(fit_joint_rmax, 
          type = "intervals", 
          prob = 0.5, 
          prob_outer = 0.89,
          point_est = "median")

conditional_effects(fit_joint_rmax, "sc_log_mean_body_mass", resp = "logrmax")
conditional_effects(fit_joint_rmax, "PC1", resp = "logrmax")
conditional_effects(fit_joint_rmax, "log_rmax", resp = "meankcv")
conditional_effects(fit_joint_rmax, "sc_MIMR", resp = "meankcv")

posterior <- as_draws_df(fit_joint_rmax)

posterior_ind <- posterior %>%
  mutate(
    indirect_bs = b_logrmax_sc_log_mean_body_mass * b_meankcv_log_rmax,
    indirect_pc1 = b_logrmax_PC1 * b_meankcv_log_rmax)

rmax_median_ci <- posterior_ind %>%
  summarise(
    lower = quantile(b_meankcv_log_rmax, 0.025),
    median = median(b_meankcv_log_rmax),
    upper = quantile(b_meankcv_log_rmax, 0.975),
    pred = "rmax")

MIMR_median_ci <- posterior_ind %>%
  summarise(
    lower = quantile(b_meankcv_sc_MIMR, 0.025),
    median = median(b_meankcv_sc_MIMR),
    upper = quantile(b_meankcv_sc_MIMR, 0.975),
    pred = "MIBMR")

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

df_median_cis <- rbind(rmax_median_ci, MIMR_median_ci, bs_median_ci, pc1_median_ci)

ggplot(df_median_cis, aes(x = median, y = pred)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0) +
  geom_point(shape = 21, stroke = 0.5, fill = "aliceblue", color = "black", size = 3) +
  # xlim(-1.25, 1.25) +
  labs(x = "Direct or Indirect Effect on CV", 
       y = "Covariate") +
  theme_classic(base_size = 14)

df_pred_pc1 <- tibble(
  log_mean_body_mass = mean(df_rmax_fs_lpi_sub_w$sc_log_mean_body_mass, na.rm = TRUE),
  PC1 = mean(df_rmax_fs_lpi_sub_w$PC1, na.rm = TRUE),
  sc_MIMR = mean(df_rmax_fs_lpi_sub_w$sc_MIMR, na.rm = TRUE),
  log_rmax = seq(min(df_rmax_fs_lpi_sub_w$log_rmax, na.rm = TRUE),
                 max(df_rmax_fs_lpi_sub_w$log_rmax, na.rm = TRUE), length.out = 100),
)

epred_pc1_cv <- add_epred_draws(fit_joint_rmax,
                                newdata = df_pred_pc1,
                                resp = "meankcv",
                                re_formula = NA)

gg_rmax_cv <- ggplot(epred_pc1_cv, aes(x = log_rmax, y = .epred)) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "black", fill = "grey60") +
  labs(x = "Growth Potential (log(rmax))", y = "Predicted CV") +
  theme_classic(base_size = 14) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme(legend.position = "none")

gg_rmax_cv










