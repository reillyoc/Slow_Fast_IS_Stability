# Living Planet Index -  Calculating CV over multiple moving windows for each .csv

#Author(s): Reilly O'Connor
#Version: 2024-01-15

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(zoo)
library(RcppRoll)
library(easystats)
library(reshape2)
library(tseries)
library(urca)
library(wql)
library(trend)

source("../R/Functions.R")

#load data
#All TS
#df_lpi <- read.csv("../R/Data/Living Planet Index/TS_lpi_all.csv", header = T)

#All TS without significant trends
#df_lpi <- read.csv("../R/Data/Living Planet Index/TS_lpi_notrend.csv", header = T)

#All TS without trends + detrended
df_lpi <- read.csv("../R/Data/TS_lpi_detrended resolved taxonomy.csv", header = T)


##### CV Calculation #####
#Filter by Maximum number of Years in TS
#df_lpi <- df_lpi %>% filter(unique_years > 19)

unique(df_lpi$Binomial)
unique(df_lpi$ID)

#List of Main IDs
main_ids <- unique(df_lpi$ID)


#Initialize the dataframe to store the results
df_lpi_cv <- data.frame(ID = integer(), cv_window = integer(), cv_value = numeric(), p_cv_value = numeric(), k_cv_value = numeric())

for (id in main_ids) {
  
  df_cv_id <- df_lpi %>%
    filter(ID == id)
  
  df_id <- df_cv_id %>% group_by(Year) %>%
    reframe(Population = sum(Value)) %>% 
    arrange(Year) %>%
    dplyr::select(Population)
  
  time_series_length <- nrow(df_id)
  
  #Calculate and store CVs for window sizes from 3 to the minimum of 5 or the time series length
  max_window_size <- min(time_series_length)
  
  for (window_size in 3:max_window_size) {
    #Calculate rolling CV for the current window size
    rolling_cv <- rollapply(data = df_id$Population, width = window_size, FUN = cv, by = 1, align = 'center', partial = F)
    rolling_p_cv <- rollapply(data = df_id$Population, width = window_size, FUN = p_cv, by = 1, align = 'center', partial = F)
    rolling_k_cv <- rollapply(data = df_id$Population, width = window_size, FUN = k_cv, by = 1, align = 'center', partial = F)
    
    mean_cv_for_window <- mean(rolling_cv, na.rm = TRUE)
    mean_p_cv_for_window <- mean(rolling_p_cv, na.rm = TRUE)
    mean_k_cv_for_window <- mean(rolling_k_cv, na.rm = TRUE)
    
    #Store the results in the dataframe
    df_lpi_cv <- rbind(df_lpi_cv, data.frame(ID = id, cv_window = window_size, cv_value = mean_cv_for_window, p_cv_value = mean_p_cv_for_window, k_cv_value = mean_k_cv_for_window))
  }
}

df_lpi_cv_wind <- df_lpi_cv %>% filter(cv_window > 4) %>%
  mutate(ID = as.character(ID)) #%>%
  #filter(cv_window < 21)

df_lpi_cv_mean <- df_lpi_cv_wind %>% group_by(cv_window) %>%
  reframe(mean_cv = mean(cv_value, na.rm = T),
          se_cv = standard_error(cv_value),
          mean_p_cv = mean(p_cv_value, na.rm = T),
          mean_k_cv = mean(k_cv_value, na.rm = T),
          se_k_cv = standard_error(k_cv_value)) #%>%
  #filter(cv_window < 21)

df_lpi_cv_window_count <- df_lpi_cv_wind %>%
  group_by(cv_window) %>%
  reframe(count = n())

summary(lm(data = df_lpi_cv_mean, mean_k_cv ~ mean_cv))

ggplot(df_lpi_cv_mean, aes(y = (mean_k_cv), x = mean_cv)) +
  geom_abline(linetype = "dashed") +
  #geom_line(alpha = 0.1) +
  geom_point(alpha = 0.5, shape = 21, color = "black") +
  ylim(0.1, 0.9) +
  xlim(0.1, 0.9) +
  geom_smooth(method = "lm")
  

ggplot(df_lpi_cv_wind, aes(y = (k_cv_value), x = as.factor(cv_window), fill = ID, group = ID)) +
  geom_line(alpha = 0.1) +
  geom_point(alpha = 0.5, shape = 21, color = "black") +
  geom_line(data = df_lpi_cv_mean, aes(y = (mean_p_cv), x = as.factor(cv_window), fill = NA, group = NA)) +
  geom_point(data = df_lpi_cv_mean, aes(y = (mean_p_cv), x = as.factor(cv_window), fill = NA, group = NA), shape = 21, size = 5, fill = "aliceblue", color = "black") +
  theme_classic() +
  #ylim(0,3) +
  #scale_y_log10(limits = c(0.001, 10)) +
  ylab("Population Varaibility (CV - log 10)") +
  xlab("Timing Window Size (Years)") +
  theme(axis.line = element_line(linetype = "solid"),
      axis.ticks = element_line(linetype = "solid"),
      panel.background = element_rect(fill = NA),
      legend.key = element_rect(fill = NA),
      legend.background = element_rect(fill = NA),
      legend.position = "none",
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 14),
      axis.title.y =element_text(size = 14), 
      axis.title.x = element_text(size = 8), 
      text = element_text(family = "Arial"))

ggplot(df_lpi_cv_mean, aes(y = (mean_k_cv), x = as.numeric(cv_window))) +
  geom_line(alpha = 0.25) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin = mean_k_cv - se_k_cv, ymax = mean_k_cv + se_k_cv), width = 0.2)

#Visually inspect all data...
main_ids <- unique(df_lpi_cv_wind$ID)
hist(qlogis(df_lpi_cv_wind$p_cv_value))

for (id in main_ids) {

  df_id <- df_lpi_cv_wind %>% filter(ID == id)
  
  gg_id <- ggplot(df_id, aes(y = (p_cv_value), x = as.numeric(cv_window), group = ID)) +
    geom_point(alpha = 0.5) +
    geom_line(alpha = 0.25) +
    theme_classic() +
    theme(legend.position = "none") +
    #scale_y_log10(limits = c(0.001, 10)) +
    ylim(0, 1) +
    xlim(0, 25)
  
  print(gg_id)
  print(unique(df_id$ID))
  readline(prompt = "Press [Enter] to continue")
}


df_lpi_info <- df_lpi %>% dplyr::select(-X, -X.1, -Year, -Value)
df_lpi_info_unique <- as.data.frame(unique(df_lpi_info))

df_lpi_cv_fin <- merge(df_lpi_cv, df_lpi_info_unique, by = "ID")

#write.csv(df_lpi_cv_fin, "../R/Data/lpi mammal population cv all.csv")
#write.csv(df_lpi_cv_fin, "../R/Data/lpi mammal population cv no trends.csv")
#write.csv(df_lpi_cv_fin, "../R/Data/lpi mammal population cv detrended resolved taxonomy.csv")


##### Calculate CV based on Lifespan ##### 
#Beccari et al 2024
df_mammal_long <- read.csv("../R/Data/Amniote database resolved taxonomy.csv", header = T)

df_mammal_fs_sub <- df_mammal_long %>% 
  dplyr::select(GBIF_ID, Species, female_maturity_d) %>%
  filter(! (female_maturity_d == -999.0000000)) %>%
  group_by(GBIF_ID) %>%
  reframe(mean_generation_yr = mean(female_maturity_d, na.rm = T)/365) %>%
  #filter(mean_generation_yr > 0.90) %>%
  mutate(generation_window = mean_generation_yr*3) %>%
  #filter(generation_window >= 3) %>%
  mutate(Window_Size = round(generation_window)) %>%
  mutate(Window_Size = ifelse(Window_Size < 3, 3, Window_Size))

df_fs_long_lpi <- merge(df_lpi, df_mammal_fs_sub, by = c("GBIF_ID"))

unique_species_lpi <- unique(df_lpi$Binomial)
unique_species_lifespan <- unique(df_fs_long_lpi$Binomial)

# Values in df1$column1 but not in df2$column2
df_sp_missing_lifespan <- as.data.frame(setdiff(unique_species_lpi, unique_species_lifespan))

unique(df_fs_long_lpi$Binomial)
unique(df_fs_long_lpi$ID)

#List of Main IDs
main_ids <- unique(df_fs_long_lpi$ID)

# Initialize the dataframe to store results
df_lpi_cv_lifespan <- data.frame(
  ID = integer(),
  cv_window = integer(),
  cv_value = numeric(),
  k_cv_value = numeric(),
  full_length = logical() # TRUE if CV was calculated over the full time series
)

# Loop over each unique ID
for (id in main_ids) {
  
  # Filter data for the current ID
  df_cv_id <- df_fs_long_lpi %>%
    filter(ID == id)
  
  # Summarize Population by Year and arrange by Year
  df_id <- df_cv_id %>%
    group_by(Year) %>%
    summarise(Population = sum(Value, na.rm = TRUE)) %>%
    arrange(Year) %>%
    dplyr::select(Population)
  
  # Get time series length
  time_series_length <- nrow(df_id)
  
  # Get the specified Window_Size from the column
  window_size <- unique(df_cv_id$Window_Size)
  
  # Check if the time series is shorter than the Window_Size
  if (time_series_length < window_size) {
    # Calculate CV for the full time series
    full_cv <- cv(df_id$Population)
    full_kcv <- k_cv(df_id$Population)
    
    # Store the result
    df_lpi_cv_lifespan <- rbind(
      df_lpi_cv_lifespan,
      data.frame(
        ID = id,
        cv_window = time_series_length, # Use the full length as the "window"
        cv_value = full_cv,
        k_cv_value = full_kcv,
        full_length = TRUE # Mark as full time series
      )
    )
  } else {
    # Calculate rolling CV using the specified Window_Size
    rolling_cv <- rollapply(
      data = df_id$Population,
      width = window_size,
      FUN = cv,
      by = 1,
      align = 'center',
      partial = TRUE
    )
    mean_cv_for_window <- mean(rolling_cv, na.rm = TRUE)
    
    rolling_kcv <- rollapply(
      data = df_id$Population,
      width = window_size,
      FUN = k_cv,
      by = 1,
      align = 'center',
      partial = TRUE
    )
    mean_kcv_for_window <- mean(rolling_kcv, na.rm = TRUE)
    
    # Store the result
    df_lpi_cv_lifespan <- rbind(
      df_lpi_cv_lifespan,
      data.frame(
        ID = id,
        cv_window = window_size,
        cv_value = mean_cv_for_window,
        k_cv_value = mean_kcv_for_window,
        full_length = FALSE # Mark as rolling window
      )
    )
  }
}

ggplot(df_lpi_cv_lifespan, aes(y = cv_value, x = as.factor(cv_window))) +
  geom_boxplot() +
  theme_classic()

ggplot(df_lpi_cv_lifespan, aes(y = k_cv_value, x = as.factor(cv_window))) +
  geom_boxplot() +
  theme_classic()

df_lpi_cv_lifespan_t <- df_lpi_cv_lifespan %>% filter(! (full_length == "TRUE"))

ggplot(df_lpi_cv_lifespan_t, aes(y = cv_value, x = as.factor(cv_window))) +
  geom_boxplot() +
  theme_classic()

ggplot(df_lpi_cv_lifespan_t, aes(y = k_cv_value, x = as.factor(cv_window))) +
  geom_boxplot() +
  theme_classic()


df_lpi_cv_lifespan_fin <- merge(df_lpi_cv_lifespan, df_lpi_info_unique, by = "ID")
unique(df_lpi_cv_lifespan_fin$ID)
unique(df_lpi_cv_lifespan_fin$Binomial)

df_lpi_cv_lifespan_t_fin <- merge(df_lpi_cv_lifespan_t, df_lpi_info_unique, by = "ID")
unique(df_lpi_cv_lifespan_t_fin$ID)
unique(df_lpi_cv_lifespan_t_fin$Binomial)

# write.csv(df_lpi_cv_lifespan_t_fin, "../R/Data/lpi mammal population cv detrended resolved taxonomy lifespan windows.csv")

