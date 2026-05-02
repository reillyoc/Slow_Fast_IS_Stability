# Linking Growth Potential to Interaction Strength in Predator-Prey Experiements

# Author(s): Reilly O'Connor
# Version: 2026-04-15

# Load Pkgs
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(loo)
library(lmodel2)

# load data
forage <- read.csv("../Slow_Fast_IS_Stability/R/Data/coblentz_etal_2025_forage_modified.csv", header = T) %>% 
  filter(Fittted.h..day. > 1e-6)

colnames(forage)[c(43,46)] <- c('Obs_a', 'Obs_h')

##### Looking at consumer 'rmax' and interaction strengths (Calculated as Gilbert/Nilsson) #####
#assume a constant conversion efficiency...
forage_IS <- forage %>%
  mutate(handling = Obs_h/Mass_g,
         clearance = Obs_a/Mass_g,
         clearance_yi = Obs_a/PredMass_g,
         handling_yi = Obs_h*PredMass_g/Mass_g,
         jmax_obs = 1/handling_yi,
         clearance_yi_pred = X.Type2a_median./PredMass_g,
         handling_yi_pred = X.Type2h_median.*PredMass_g/Mass_g,
         jmax_pred = 1/handling_yi_pred,
         amax_obs = 1/handling,
         a_h_obs = clearance * handling,
         b_obs = 1/a_h_obs,
         eff = 0.80,
         E = 14,
         Mn = Mass_g,
         Mp = PredMass_g,
         m = (PredMetabolism/E),
         X = m/PredMass_g,
         a_obs = clearance,
         h_obs = handling,
         a_yi = clearance_yi,
         h_yi = handling_yi,
         b_yi = 1/(clearance_yi*handling_yi),
         b_yi_pred = 1/(clearance_yi_pred*handling_yi_pred),
         K = PreyAbundance_90*Mass_g) %>%
  mutate(# IS_yi = K * a_yi * (((eff)/(m/PredMass_g)) - h_yi),
    # IS_obs = K * a_obs * (((eff)/(m)) - h_obs),
    IS_obs = K * ((jmax_obs * eff) - X) / X * b_yi,
    IS_pred = K * ((jmax_pred * eff) - X) / X * b_yi,
    max_consump = jmax_obs * K * eff,
    IS_pred_2 = (Mp/Mn) * ((eff * jmax_obs / X) - X/X),
    # growth_yi = ((eff*jmax_obs) - (m)),
    # rmax_yi = (growth_yi/Mp),
    eco_scope = b_yi/K,
    rmax_met = eff*jmax_obs/X,
    growth_obs = ((eff*jmax_obs) - (X)),
    rmax_obs = (growth_obs)) %>%
  filter(IS_obs > 0 & rmax_obs > 0 & IS_pred > 0
         # & IS_yi > 0 & rmax_yi > 0
  ) %>%
  # filter(! Dim == 2.5) %>%
  mutate(log_pred_body_size = log10(Predator.mass..mg.),
         log_prey_body_size = log10(Prey.mass..mg.),
         log_size_ratio = log10(Predator.mass..mg./Prey.mass..mg.),
         log_K = log10(K),
         log_b_yi = log10(b_yi),
         log_a_obs = log10(a_obs),
         log_h_obs = log10(h_obs),
         log_amax_obs = log10(amax_obs),
         log_jmax_obs = log10(jmax_obs),
         log_X = log10(X),
         log_max_consump = log10(max_consump),
         log_IS_obs = log10(IS_obs),
         log_IS_pred = log10(IS_pred),
         log_IS_obs_mass = log10(IS_obs/Mass_g),
         log_IS_pred_2 = log10(IS_pred_2),
         log_eco_scope = log10(eco_scope),
         log_rmax_met = log10(rmax_met),
         # log_IS_yi = log10(IS_yi),
         # log_IS_Nils_obs = log10(IS_Nils_obs),
         # log_growth_yi = log10(growth_yi),
         # log_rmax_yi = log10(rmax_yi),
         log_growth_obs = log10(growth_obs),
         log_rmax_obs = log10(rmax_obs)) %>%
  filter(log_b_yi < 4)

forage_IS_vertebrates <- forage_IS %>%
  filter(Dim == 2) %>%
  mutate(Dim = as.factor(Dim),
         Habitat, as.factor(Habitat),
         Vert.invert = as.factor(Vert.invert),
         Vert.invert.1 = as.factor(Vert.invert.1),
         VertPair = interaction(Vert.invert, Vert.invert.1, sep = "_")) %>%
  filter(!(Data.set == 1954 | Data.set == 1957))

unique(forage_IS_vertebrates$Data.set)

lm_test<- lm(data = forage_IS_vertebrates, log_rmax_met ~ log_size_ratio)
summary(lm_test)
plot(data = forage_IS_vertebrates, log_rmax_met ~ log_size_ratio)
abline(lm_test)

rma_K_b <- lmodel2(log_b_yi ~ log_K, 
                   data = forage_IS_vertebrates, 
                   range.y = "interval", 
                   range.x = "interval")

rma_K_b$regression.results
rma_K_b$rsquare
rma_K_b$confidence.intervals

##### BRMS regression models #####
##### Predator Body Size - Rmax #####
priors <- c(
  set_prior("normal(0, 1)", class = "b"),          
  set_prior("normal(0, 2)", class = "Intercept"),  
  set_prior("student_t(3, 0, 1)", class = "sigma") 
)

rmx_pred_bs <- brm(
  formula = log_rmax_obs ~ log_pred_body_size,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_rmx_pred_bs <- loo(rmx_pred_bs)
loo_rmx_pred_bs

r2_bayes(rmx_pred_bs)

rmx_pred_bs
pp_check(rmx_pred_bs, ndraws = 100)
bayes_R2(rmx_pred_bs, probs = c(0.055, 0.945))
summary(rmx_pred_bs, prob = 0.89)
plot(rmx_pred_bs)

#Plot predicted rmax
df_pred_rmx <- tibble(
  log_pred_body_size = seq(min(forage_IS_vertebrates$log_pred_body_size),
                           max(forage_IS_vertebrates$log_pred_body_size),
                           length.out = 200)
)

epred_draws_rmx <- add_epred_draws(rmx_pred_bs, newdata = df_pred_rmx)

gg_rmx_bs <- ggplot(epred_draws_rmx, aes(x = log_pred_body_size, y = .epred)) +
  geom_point(data = forage_IS_vertebrates, 
             aes(x = log_pred_body_size, y = log_rmax_obs), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "Log Predator Body Size (mg)", 
       y = "Log Growth Potential (rmax)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_rmx_bs

##### Interaction Strength - Rmax #####
IS_pred_rmax <- brm(
  formula = log_IS_obs ~ log_rmax_obs,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_IS_pred_rmax <- loo(IS_pred_rmax)
loo_IS_pred_rmax

r2_bayes(IS_pred_rmax)
pp_check(IS_pred_rmax, ndraws = 100)
bayes_R2(IS_pred_rmax, probs = c(0.055, 0.945))
summary(IS_pred_rmax, prob = 0.89)
plot(IS_pred_rmax)

# Plot predicted IS ~ rmax
df_pred_IS <- tibble(
  log_rmax_obs = seq(min(forage_IS_vertebrates$log_rmax_obs),
                     max(forage_IS_vertebrates$log_rmax_obs),
                     length.out = 200)
)

epred_draws_IS <- add_epred_draws(IS_pred_rmax, newdata = df_pred_IS)

gg_IS_rmax <- ggplot(epred_draws_IS, aes(x = log_rmax_obs, y = .epred)) +
  geom_point(data = forage_IS_vertebrates,
             aes(x = log_rmax_obs, y = log_IS_obs), size = 2.5, 
             # stroke = 0.5, shape = 21,
             color = "black", #fill = "aliceblue",
             alpha = 0.15) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "blue3", fill = "grey") +
  labs(x = "Log Growth Potential (rmax)",
       y = "Log Relative Energy Flux (IS)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_IS_rmax


#Figure 5
gg_rmax_IS_grid <- plot_grid(gg_rmx_bs, gg_IS_rmax, 
                             align = "hv", nrow = 1,
                             labels = c("a", "b"))
gg_rmax_IS_grid


# ggsave("../Slow_Fast_IS_Stability/R/Figures/Figure 5 - Bayesian Rmax - Body Size - IS.jpeg", plot = gg_rmax_IS_grid, width = 8, height = 4)

##### Interaction Strength - Size Ratio #####
IS_size_ratio <- brm(
  formula = log_IS_obs ~ log_size_ratio,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_IS_size_ratio <- loo(IS_size_ratio)
loo_IS_size_ratio

r2_bayes(IS_size_ratio)
pp_check(IS_size_ratio, ndraws = 100)
bayes_R2(IS_size_ratio, probs = c(0.055, 0.945))
summary(IS_size_ratio, prob = 0.89)
plot(IS_size_ratio)

df_pred_IS_sr <- tibble(
  log_size_ratio = seq(min(forage_IS_vertebrates$log_size_ratio),
                       max(forage_IS_vertebrates$log_size_ratio),
                       length.out = 200)
)

epred_draws_IS_sr <- add_epred_draws(IS_size_ratio, newdata = df_pred_IS_sr)

gg_IS_sr <- ggplot(epred_draws_IS_sr, aes(x = log_size_ratio, y = .epred)) +
  geom_point(data = forage_IS_vertebrates,
             aes(x = log_size_ratio, y = log_IS_obs),
             shape = 21, size = 2.5, stroke = 0.5,
             color = "black", fill = "grey15", alpha = 0.25) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "black", fill = "grey") +
  labs(x = "Log Predator:Prey Size Ratio",
       y = "Log Relative Energy Flux (IS)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_IS_sr

##### Interaction Strength - Predator Body Size ##### 
IS_pred_bs <- brm(
  formula = log_IS_obs ~ log_pred_body_size,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

df_pred_IS_pbs <- tibble(
  log_pred_body_size = seq(min(forage_IS_vertebrates$log_pred_body_size),
                           max(forage_IS_vertebrates$log_pred_body_size),
                           length.out = 200)
)

epred_draws_IS_pbs <- add_epred_draws(IS_pred_bs, newdata = df_pred_IS_pbs)

gg_IS_pbs <- ggplot(epred_draws_IS_pbs, aes(x = log_pred_body_size, y = .epred)) +
  geom_point(data = forage_IS_vertebrates,
             aes(x = log_pred_body_size, y = log_IS_obs),
             shape = 21, size = 2.5, stroke = 0.5,
             color = "black", fill = "grey15", alpha = 0.25) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "black", fill = "grey") +
  labs(x = "Log Predator Body Size (mg)",
       y = "Log Relative Energy Flux (IS)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")


##### Interaction Strength - Prey Body Size ##### 
IS_prey_bs <- brm(
  formula = log_IS_obs ~ log_prey_body_size,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

df_pred_IS_prbs <- tibble(
  log_prey_body_size = seq(min(forage_IS_vertebrates$log_prey_body_size),
                           max(forage_IS_vertebrates$log_prey_body_size),
                           length.out = 200)
)

epred_draws_IS_prbs <- add_epred_draws(IS_prey_bs, newdata = df_pred_IS_prbs)

gg_IS_prbs <- ggplot(epred_draws_IS_prbs, aes(x = log_prey_body_size, y = .epred)) +
  geom_point(data = forage_IS_vertebrates,
             aes(x = log_prey_body_size, y = log_IS_obs),
             shape = 21, size = 2.5, stroke = 0.5,
             color = "black", fill = "grey15", alpha = 0.25) +
  stat_lineribbon(.width = 0.89, alpha = 0.65, linewidth = 2.0,
                  color = "black", fill = "grey") +
  labs(x = "Log Prey Body Size (mg)",
       y = "Log Relative Energy Flux (IS)") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")


##### Combined grid ######
gg_IS_grid <- plot_grid(gg_IS_sr, gg_IS_pbs, gg_IS_prbs,
                        align = "hv", nrow = 1,
                        labels = c("a", "b", "c"))
gg_IS_grid


##### rmax ~ Size Ratio Model #####
rmx_size_ratio <- brm(
  formula = log_rmax_obs ~ log_size_ratio,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

loo_rmx_size_ratio <- loo(rmx_size_ratio)
loo_rmx_size_ratio

r2_bayes(rmx_size_ratio)
pp_check(rmx_size_ratio, ndraws = 100)
bayes_R2(rmx_size_ratio, probs = c(0.055, 0.945))
summary(rmx_size_ratio, prob = 0.89)
plot(rmx_size_ratio)


IS_prey_bs <- brm(
  formula = log_IS_obs ~ log_prey_body_size,
  data = forage_IS_vertebrates,
  prior = priors,
  backend = 'cmdstanr',
  family = gaussian(),
  chains = 4, cores = 4,
  warmup = 2000, iter = 4000,
  control = list(adapt_delta = 0.95),
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)


##### Coefficient Forest Plots #####
# Standardise predictors for comparable coefficients across forest plots
forage_IS_vertebrates <- forage_IS_vertebrates %>%
  mutate(
    sc_log_pred_body_size = scale(log_pred_body_size)[,1],
    sc_log_prey_body_size = scale(log_prey_body_size)[,1],
    sc_log_size_ratio     = scale(log_size_ratio)[,1],
    sc_log_rmax_obs       = scale(log_rmax_obs)[,1]
  )

# Refit all models on standardised predictors
rmx_pred_bs_sc   <- update(rmx_pred_bs, newdata = forage_IS_vertebrates, formula = log_rmax_obs ~ sc_log_pred_body_size)
rmx_size_ratio_sc <- update(rmx_size_ratio, newdata = forage_IS_vertebrates, formula = log_rmax_obs ~ sc_log_size_ratio)

IS_pred_rmax_sc  <- update(IS_pred_rmax, newdata = forage_IS_vertebrates, formula = log_IS_obs ~ sc_log_rmax_obs)
IS_pred_bs_sc    <- update(IS_pred_bs, newdata = forage_IS_vertebrates, formula = log_IS_obs ~ sc_log_pred_body_size)
IS_prey_bs_sc    <- update(IS_prey_bs, newdata = forage_IS_vertebrates, formula = log_IS_obs ~ sc_log_prey_body_size)
IS_size_ratio_sc <- update(IS_size_ratio, newdata = forage_IS_vertebrates, formula = log_IS_obs ~ sc_log_size_ratio)

# Extract slope posteriors — rmax models
draws_rmax_coef <- bind_rows(
  as_draws_df(rmx_pred_bs_sc) %>%
    transmute(estimate = b_sc_log_pred_body_size,
              predictor = "Predator Body Size"),
  as_draws_df(rmx_size_ratio_sc) %>%
    transmute(estimate = b_sc_log_size_ratio,
              predictor = "Pred:Prey Size Ratio")
) %>% mutate(response = "log(rmax)")

# Extract slope posteriors — IS models
draws_IS_coef <- bind_rows(
  as_draws_df(IS_pred_rmax_sc) %>%
    transmute(estimate = b_sc_log_rmax_obs,
              predictor = "rmax"),
  as_draws_df(IS_pred_bs_sc) %>%
    transmute(estimate = b_sc_log_pred_body_size,
              predictor = "Predator Body Size"),
  as_draws_df(IS_prey_bs_sc) %>%
    transmute(estimate = b_sc_log_prey_body_size,
              predictor = "Prey Body Size"),
  as_draws_df(IS_size_ratio_sc) %>%
    transmute(estimate = b_sc_log_size_ratio,
              predictor = "Pred:Prey Size Ratio")
) %>% mutate(response = "log(IS)")

# Forest plot — rmax
gg_coef_rmax <- ggplot(draws_rmax_coef,
                        aes(x = estimate, y = predictor, fill = predictor)) +
  stat_halfeye(alpha = 0.8, .width = 0.89) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.75) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Standardised Coefficient", y = NULL,
       title = "log(rmax)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_coef_rmax

# Forest plot — IS
gg_coef_IS <- ggplot(draws_IS_coef,
                      aes(x = estimate, y = predictor, fill = predictor)) +
  stat_halfeye(alpha = 0.8, .width = 0.89) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.75) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Standardised Coefficient", y = NULL,
       title = "log(IS)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

gg_coef_IS

# Combined forest plot grid
gg_coef_grid <- plot_grid(gg_coef_rmax, gg_coef_IS,
                           nrow = 1, align = "hv",
                           labels = c("a", "b"),
                           rel_widths = c(1, 1.5))
gg_coef_grid

# ggsave("../R/Figures/Figure SX - Coefficient Forest Plots.jpeg", plot = gg_coef_grid, width = 10, height = 5)








##### OLS and Type 2 Regressions #####
lm_rmax_bs_obs <- lm(log_rmax_obs ~ log_pred_body_size, data = forage_IS_vertebrates)
summary(lm_rmax_bs_obs)

rma_rmax_bs_obs <- lmodel2(log_rmax_obs ~ log_pred_body_size, data = forage_IS_vertebrates, range.y = "interval", range.x = "interval")
rma_rmax_bs_obs$regression.results
rma_rmax_bs_obs$rsquare

gg_rmax_bs_obs <- ggplot(forage_IS_vertebrates, aes(x = log_pred_body_size, y = log_rmax_obs)) +
  geom_point(shape = 21, color = "black", fill = "aliceblue", stroke = 0.25, size = 3) +
  # geom_point() +
  geom_smooth(method = "lm", se = F, linetype = "solid", linewidth = 2, color = "black", alpha = 0.75) +
  # geom_abline(slope = -0.2556172, intercept = 0.4883376, linewidth = 2, color = "black", alpha = 0.75) +
  # facet_wrap(~ VertPair) +
  labs(x = "Log Predator Body Size",
       y = "Log Growth Potential (rmax)") +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  theme_bw(base_size = 16)

gg_rmax_bs_obs

lm_is_rmax_obs <- lm(log_IS_obs ~ log_rmax_obs, data = forage_IS_vertebrates)
summary(lm_is_rmax_obs)

rma_rmax_IS <- lmodel2(log_IS_obs ~ log_rmax_obs, data = forage_IS_vertebrates, range.y = "interval", range.x = "interval")
rma_rmax_IS$regression.results
rma_rmax_IS$rsquare

gg_IS_rmax_obs <- ggplot(forage_IS_vertebrates, aes(x = log_rmax_obs, y = log_IS_obs)) +
  geom_point(shape = 21, color = "black", fill = "aliceblue", stroke = 0.25, size = 3) +
  # geom_point() +
  geom_smooth(method = "lm", se = F, linetype = "solid", linewidth = 2, color = "black", alpha = 0.75) +
  # geom_abline(slope = 2.047541, intercept = 3.496617, linewidth = 2, color = "black", alpha = 0.75) +
  # facet_wrap(~ VertPair) +
  labs(x = "Log Growth Potential (rmax)",
       y = "Log Relative Energy Flux (IS max)") +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  theme_bw(base_size = 16)

gg_IS_rmax_obs

gg_rmax_IS_grid <- plot_grid(gg_rmax_bs_obs, gg_IS_rmax_obs, align = "hv", nrow = 1)
gg_rmax_IS_grid

# ggsave("../KyleCoblentz-FR_Prediction-3c8e63c/Outputs/Figure SX, rmax vs body size, rmax vs Relative Energy Flux.jpeg", plot = gg_rmax_IS_grid, width = 10, height = 5)

rma_IS_bs_obs <- lmodel2(log_IS_obs ~ log_size_ratio, data = forage_IS_vertebrates, range.y = "interval", range.x = "interval")
rma_IS_bs_obs$regression.results
rma_IS_bs_obs$rsquare

gg_IS_bs_obs <- ggplot(forage_IS_vertebrates, aes(x = log_size_ratio, y = log_IS_obs)) +
  geom_point(shape = 21, color = "black", fill = "aliceblue", stroke = 0.25, size = 3) +
  # geom_point() +
  geom_smooth(method = "lm", se = F, linetype = "solid", linewidth = 2, color = "black", alpha = 0.75) +
  # geom_abline(slope = -0.4582887, intercept = 3.977773, linewidth = 2, color = "black", alpha = 0.75) +
  # facet_wrap(~ VertPair) +
  labs(x = "Log Predator:Prey Body Size Ratio",
       y = "Log Relative Energy Flux (ISmax)") +
  scale_y_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(), expand = expansion(mult = 0.10)) +
  theme_bw(base_size = 16)

gg_IS_bs_obs

# ggsave("../KyleCoblentz-FR_Prediction-3c8e63c/Outputs/Figure SX, Relative Energy Flux vs Pred Prey Body Size Ratio.jpeg", plot = gg_IS_bs_obs, width = 6, height = 6)

