#Living Planet Index - Checking for suitable time series and detrending

#Author(s): Reilly O'Connor
#Version: 2024-12-15

#Pkgs
library(NSM3)
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

#load data
df_lpi <- read.csv("../Slow_Fast_IS_Stability/R/Data/Living Planet Index/LPD2022_public/LPD2022_public.csv", header = T)

##### Code #####
##### Organize Data and Subset for Mammals #####
#Filter LPI for only mammal data
df_lpi_mammals <- df_lpi %>% filter(Class == "Mammalia")

#Transpose dataframe
df_lpi_melt <- melt(df_lpi_mammals, 
                  id.vars = c("ID", "Binomial", "Replicate", "Excluded_LPR_2022", "Citation", 
                              "Class", "Order", "Family", "Genus", "Species", "Subspecies", 
                              "Authority", "Common_name", "Location", "Country", "All_countries", 
                              "Region", "IPBES_region", "IPBES_subregion", "Latitude", "Longitude", 
                              "Specific_location", "System", "T_realm", "T_biome", "FW_realm", 
                              "FW_biome", "M_realm", "M_ocean", "M_biome", "Migratory_fw_fish", 
                              "Units", "Method"),
                  variable.name = "Year",
                  value.name = "Value")

#remove X and null values, set year as numerical
df_lpi_melt$Year <- sub("^X", "", df_lpi_melt$Year)
df_lpi_mammal_melt <- df_lpi_melt %>% filter(! Value == "NULL")
df_lpi_mammal_melt$Year <- as.numeric(df_lpi_mammal_melt$Year)

##### Assess TS for Appropriate number of continuous years (minimum 10 or more) #####
#Summarize study years
df_unique_years <- df_lpi_mammal_melt %>% 
  group_by(ID) %>%
  reframe(start_year = min(Year),
            end_year = max(Year),
            unique_years = n_distinct(Year),
            is_continuous = (max(Year) - min(Year) + 1) == n_distinct(Year))

unique(df_unique_years$ID)

#Initial filter based on years
df_lpi_min <- df_unique_years %>% 
  filter(unique_years > 4, is_continuous == TRUE) %>%
  group_by(ID) %>%
  arrange(ID)

#Studies that meet more than 10 years but are not continuous to potentially add back in
df_lpi_potential <- df_unique_years %>% 
  filter(unique_years > 4, is_continuous == FALSE) %>%
  group_by(ID) %>%
  arrange(ID)

df_sub_lpi_f <- merge(df_lpi_mammal_melt, df_lpi_potential, by = "ID")

df_sub_lpi_f <- df_sub_lpi_f %>%
  arrange(ID, Year)

#Function to find longest continuous section of TS
find_longest_continuous <- function(sub_df) {
  sub_df <- sub_df %>% arrange(Year)
  sub_df$diff <- c(1, diff(sub_df$Year))
  sub_df$group <- cumsum(sub_df$diff != 1)
  longest_group <- sub_df %>% group_by(group) %>% 
    summarise(length = n()) %>% 
    arrange(desc(length)) %>% 
    slice(1)
  sub_df %>% filter(group == longest_group$group) %>% 
    dplyr::select(-diff, -group)
}

#Apply function to each non-continuous time series and combine results
result <- df_sub_lpi_f %>% group_by(ID) %>% 
  do(find_longest_continuous(.)) %>%
  dplyr::select(-start_year, -end_year, -unique_years, -is_continuous)

#Add a column with the length of each series and filter for time series that have 5 years or more of data
result_filtered <- result  %>% 
  group_by(ID) %>%
  reframe(start_year = min(Year),
          end_year = max(Year),
          unique_years = n_distinct(Year),
          is_continuous = (max(Year) - min(Year) + 1) == n_distinct(Year)) %>%
  filter(unique_years > 4)

df_lpi_recovered <- merge(result, result_filtered, by = "ID")

df_lpi_recovered <- df_lpi_recovered %>%
  arrange(ID, Year)

df_sub_lpi_true <- merge(df_lpi_mammal_melt, df_lpi_min, by = "ID")

df_sub_lpi_true <- df_sub_lpi_true %>% group_by(ID) %>%
  arrange(Year)

df_sub_lpi_tf <- rbind(df_sub_lpi_true, df_lpi_recovered)

df_sub_lpi_tf <- df_sub_lpi_tf %>%
  arrange(ID, Year)

#Change value to numeric
df_sub_lpi_tf$Value <- as.numeric(df_sub_lpi_tf$Value)

unique(df_sub_lpi_tf$Binomial)
unique(df_sub_lpi_tf$ID)

##### Data to remove based on invalid methods or unknowns in/about data #####
df_lpi_true_rm <- df_sub_lpi_tf %>% filter(
  
  #All TS with 10 or more years of Data to be removed
        ! (ID == 262 | ID == 506 | ID == 562 | ID == 563 | ID == 566 | ID == 567 | ID == 568 |
           ID == 1962 | ID == 1967 | ID == 1974 | ID == 1975 | ID == 1980 | ID == 1983 |ID == 1985 |
           ID == 1993 | ID == 1997 | ID == 2008 | ID == 2197 | ID == 2234 | ID == 2356 | ID == 2359 | ID == 3416 | 
           ID == 3417 | ID == 3419 | ID == 3420 | ID == 3435 | ID == 3436 | ID == 3437 | ID == 3438 | ID == 3439 | 
           ID == 3440 | ID == 3441 | ID == 3442 | ID == 3443 | ID == 3456 | ID == 3457 | ID == 3458 | ID == 3459 | 
           ID == 3460 | ID == 3461 | ID == 3462 | ID == 3463 | ID == 3464 | ID == 3484 | ID == 3486 | ID == 3487 | 
           ID == 3488 | ID == 3489 | ID == 3490 | ID == 3491 | ID == 3492 | ID == 3493 | ID == 3494 | ID == 4517 | 
           ID == 4703 | ID == 4835 | ID == 4836 | ID == 4837 | ID == 4838 | ID == 4839 | ID == 4840 | ID == 4841 | 
           ID == 4842 | ID == 4843 | ID == 4963 | ID == 5689 | ID == 5817 | ID == 5818 | ID == 5827 | ID == 5828 | 
           ID == 5839 | ID == 5852 | ID == 5876 | ID == 6227 | ID == 6339 | ID == 8103 | ID == 8221 | ID == 8226 | 
           ID == 8321 | ID == 8322 | ID == 8323 | ID == 8421 | ID == 9104 | ID == 9134 | ID == 10694 | ID == 10695 | 
           ID == 10696 | ID == 10740 | ID == 10776 | ID == 10778 | ID == 10871 | ID == 10872 | ID == 10873 | ID == 11103 |
           ID == 11104 | ID == 11105 | ID == 11175 | ID == 11473 | ID == 11502 | ID == 12438 | ID == 12445 | ID == 12685 | 
           ID == 13287 | ID == 14553 | ID == 17787 | ID == 17815 | ID == 17816 | ID == 19088 | ID == 19722 | ID == 19802 | 
           ID == 23279 | ID == 27413 | ID == 27585 | ID == 27590 | ID == 27741 | ID == 27747 | ID == 27789 | ID == 27860 | 
           ID == 108092 | ID == 202588 | ID == 202589 | ID == 117 | ID == 1952 | ID == 1953 | ID == 1954 | ID == 1956 | 
           ID == 1957 | ID == 1958 | ID == 1959 | ID == 1960 | ID == 1963 | ID == 1964 | ID == 1968 | ID == 1969 |
           ID == 1970 | ID == 1972 | ID == 1973 | ID == 1976 | ID == 1977 | ID == 1978 | ID == 1981 | ID == 1982 | 
           ID == 1984 | ID == 1986 | ID == 1987 | ID == 1989 | ID == 1990 | ID == 1994 | ID == 1995 | ID == 1996 | 
           ID == 1998 | ID == 1999 | ID == 2000 | ID == 2001 | ID == 2002 | ID == 2003 | ID == 2004 | ID == 2006 | 
           ID == 2200 | ID == 3428 | ID == 3480 | ID == 3488 | ID == 4625 | ID == 4633 | ID == 4903 | ID == 5085 | 
           ID == 5086 | ID == 5087 | ID == 5088 | ID == 6394 | ID == 6395 | ID == 7353 | ID == 8147 | ID == 8148 |
           ID == 9044 | ID == 9126 | ID == 9234 | ID == 10692 | ID == 11651 | ID == 11886 | ID == 11887 | ID == 18290 | 
           ID == 27629 | ID == 27632 | ID == 27634 | ID == 27635 | ID == 27637 | ID == 27654 | ID == 27673 | ID == 27678 | 
           ID == 106552 | ID == 107163 | ID == 107165 | ID == 107169 | ID == 107170 | ID == 107172 | ID == 107177 | ID == 107178 | 
           ID == 107179 | ID == 108101 | ID == 108112 | ID == 108113 | ID == 108117 | ID == 108118 | ID == 108119 | ID == 108121 |    
           ID == 108122 | ID == 500004 | ID == 19075 | ID == 19076| ID == 19081| ID == 19084| ID == 19086| ID == 19091| ID == 21296| 
           ID == 21297| ID == 21298| ID == 21299 | ID == 21300| ID == 21302| ID == 21303| ID == 21304| ID == 21305| ID == 107175 | 
           ID == 107162| ID == 27645| ID == 11147 | ID == 27880 | ID == 27879 |
             
           #log/ln TS
           ID == 5823 | ID == 5855 | ID == 6261 | ID == 6291 | ID == 6292 | ID == 6451 | ID == 14544 | ID == 5840 | ID == 11747 | ID == 11748 |
            
    #All TS with 5 or More years of Data to be removed
           ID == 507 | ID == 510 | ID == 511 | ID == 1355 | ID == 1966 | ID == 1979 | ID == 2226 | ID == 2357 | ID == 2364 |
           ID == 2367 | ID == 3418 | ID == 4710 | ID == 5269 | ID == 5270 | ID == 5494 | ID == 5754 | ID == 5836 | ID == 5837 |
           ID == 5854 | ID == 5865 | ID == 5866 | ID == 5868 | ID == 5869 | ID == 6073 | ID == 6074 | ID == 6353 | ID == 6468 |
           ID == 6550 | ID == 9051 | ID == 9052 | ID == 9114 | ID == 9309 | ID == 9813 | ID == 10468 | ID == 10473 | ID == 10474 |
           ID == 10555 | ID == 10706 | ID == 10723 | ID == 10730 | ID == 10731 | ID == 10732 | ID == 11166 | ID == 11276 | ID == 11278 |
           ID == 11284 | ID == 11297 | ID == 11901 | ID == 12534 | ID == 13483 | ID == 13490 | ID == 17826 | ID == 18004 | ID == 18256 | ID == 18311 |
           ID == 18877 | ID == 18878 | ID == 18880 | ID == 18882 | ID == 18885 | ID == 19514 | ID == 19519 | ID == 19520 | ID == 19522 |
           ID == 19523 | ID == 19550 | ID == 19726)) %>%
  
  #Remove endangered spp
      filter(
        !(Binomial == "Burramys_parvus" |
           Binomial == "Bettongia_penicillata" |	
           Binomial == "Gymnobelideus_leadbeateri" |
           Binomial == "Moschus_moschiferus" |
           Binomial == "Dasyurus_geoffroii" |
           Binomial == "Pan_troglodytes" |
           Binomial == "Panthera_tigris" |
           Binomial == "Parantechinus_apicalis"|
           Binomial == "Platanista_minor" |
           Binomial == "Pongo_pygmaeus" |
      str_detect(Citation, "^Threatened Species Index data portal:")))

#Rename Binomial column after removing unsuitable time series
df_lpi_final_sp <- df_lpi_true_rm %>% mutate(new_species = str_replace(Binomial, "_", " ")) %>%
  arrange(ID, Year)

unique(df_lpi_final_sp$Binomial)
unique(df_lpi_final_sp$ID)

##### Evaluate TS for significant trends over time #####
#List of Unique TS IDs
unique_ids <- unique(df_lpi_final_sp$ID)

#Empty dataframe for results
results <- data.frame(ID = integer(), PValue_lm = numeric(), PValue_mk = numeric(), DetrendNeeded_lm = logical(), DetrendNeeded_mk = logical())

#Loop through each TS ID
for(id in unique_ids) {
  #Extract the time series for the current ID
  time_series_data <- subset(df_lpi_final_sp, ID == id)
  print(unique(time_series_data$ID))
  
  #Linear & Mann Kendall Tests for trends
  lm_model <- lm(Value ~ Year, data = time_series_data)
  p_value_lm <- summary(lm_model)$coefficients[2,4]
  
  mk_model <- mannKen(time_series_data$Value)
  p_value_mk <- mk_model$p.value
  
  
  #Determine if significant trend in TS
  detrend_needed_lm <- p_value_lm < 0.05
  detrend_needed_mk <- p_value_mk < 0.05
  
  #Append the results
  results <- rbind(results, data.frame(ID = id, PValue_lm = p_value_lm, PValue_mk = p_value_mk, DetrendNeeded_lm = detrend_needed_lm, DetrendNeeded_mk = detrend_needed_mk))
}

###### Visually Inspect Each Time Series to Assess Validity #####
#List of Main IDs
df_lpi_trend_all <- results %>% dplyr::select(ID, DetrendNeeded_mk)

df_lpi_ts_vis <- merge(df_lpi_final_sp, df_lpi_trend_all, by = "ID")

df_lpi_ts_vis <- df_lpi_ts_vis %>% arrange(ID, Year)

main_ids <- unique(df_lpi_ts_vis$ID)

for (id in main_ids) {
  #Filter the data for the current ID
  id_data <- subset(df_lpi_ts_vis, ID == id)
  Common_name <- unique(id_data$Common_name)
  Trend <- unique(id_data$DetrendNeeded_mk)
  
  #Plotting the data
  plot(id_data$Year, id_data$Value, main = paste("Time Series for ID:", id, "-", Common_name, "-", Trend), 
       xlab = "Year", ylab = "Value", pch = 19, col = "blue")
  lines(id_data$Year, id_data$Value, col = "black", type = "l")
  
  print(unique(id_data$Binomial))
  readline(prompt = "Press [Enter] to continue")
}

df_lpi_ts_ids <- merge(df_lpi_final_sp, results, by = "ID")

##### Remove TS IDs that do not pass visual inspection #####
#ID 9863 - many 0s and all repetitive values after 2003
#IDs 9896 + 9899 - data is identical and source cannot be located for correct values
#ID 27831, 27832, 27833 All values nearly identical
#ID 23155 - TS is mostly 0s
#ID 1991, 1992, 4860 perfect fits - modelled data remove

df_lpi_ts_rm <- df_lpi_ts_ids %>%
  filter(! (ID == 1991 | ID == 1992 | ID == 4860 | ID == 9863| ID == 9896| ID == 9899| ID == 23155 | ID == 27832 | ID == 27831))

##### TS that need years shaved off for relative stationarity (populations with initial re-introduction phase from around 0, steep decline to near 0 at end of TS or with more than three repetitive values at start or end) #####
#list of IDs that need the years reduced
list_ts_spp <- list(
  "28" = 1975:1994,
  "3431" = 1974:1988,
  "9055" = 1980:1992,
  "9864" = 1996:2007,
  "9883" = 1990:2001,
  "10734" = 1971:2008,
  "11143" = 1977:1996,
  "11148" = 1971:1993,
  "27571" = 1972:1993,
  "27714" = 1992:2008,
  "27718" = 1990:2012
)

df_lpi_ts_rm_st <- df_lpi_ts_rm %>%
  filter(case_when(
    ID == 28 ~ Year %in% list_ts_spp[["28"]],
    ID == 3431 ~ Year %in% list_ts_spp[["3431"]],
    ID == 9055 ~ Year %in% list_ts_spp[["9055"]],
    ID == 9864 ~ Year %in% list_ts_spp[["9864"]],
    ID == 9883 ~ Year %in% list_ts_spp[["9883"]],
    ID == 10734 ~ Year %in% list_ts_spp[["10734"]],
    ID == 11143 ~ Year %in% list_ts_spp[["11143"]],
    ID == 11148 ~ Year %in% list_ts_spp[["11148"]],
    ID == 27571 ~ Year %in% list_ts_spp[["27571"]],
    ID == 27714 ~ Year %in% list_ts_spp[["27714"]],
    ID == 27718 ~ Year %in% list_ts_spp[["27718"]],
    TRUE ~ TRUE))

unique(df_lpi_ts_rm_st$new_species)
unique(df_lpi_ts_rm_st$ID)

#write.csv(df_lpi_ts_rm_st, "../R/Data/Living Planet Index/TS_lpi_all_trended.csv")

#CSV for data with no sig trends
df_lpi_notrend <- df_lpi_ts_rm_st %>% filter(DetrendNeeded_mk == "FALSE")

df_lpi_notrend <- df_lpi_notrend %>% arrange(ID, Year)

unique(df_lpi_notrend$new_species)
unique(df_lpi_notrend$ID)

#write.csv(df_lpi_notrend, "../R/Data/Living Planet Index/TS_lpi_notrend.csv")

##### Detrend trended data based on Theil-Sen Slope #####
df_lpi_ts_detrend <- df_lpi_ts_rm_st %>% filter(DetrendNeeded_mk == TRUE)

#Initialize an empty data frame to store the results
detrend_results <- data.frame(ID = integer(), Year = numeric(), Value = numeric())

main_ids <- unique(df_lpi_ts_detrend$ID)

#Loop through each TS ID
for(id in main_ids) {
  #Extract the time series for the current ID
  time_series_data <- subset(df_lpi_ts_detrend, ID == id)
  print(unique(time_series_data$ID))
  
  #Fit Theil-Sen to calculate Residuals and Detrend data
  model_fit <- mblm(Value ~ Year, repeated = F, data = time_series_data)
  summary(model_fit)
  
  mean_value <- mean(time_series_data$Value)
  residual_value <- (mean_value + resid(model_fit))
  years <- time_series_data$Year
  
  #Append the results
  detrend_results <- rbind(detrend_results, data.frame(ID = id, Year = years, Value = residual_value))
}

#Visually inspect all data...
main_ids <- unique(detrend_results$ID)

for (id in main_ids) {
  #Filter the data for the current ID
  id_data <- subset(detrend_results, ID == id)
  Common_name <- unique(id_data$Common_name)
  Trend <- unique(id_data$DetrendNeeded_mk)
  
  #Plotting the data
  plot(id_data$Year, id_data$Value, main = paste("Time Series for ID:", id, "-", Common_name, "-", Trend), 
       xlab = "Year", ylab = "Value", pch = 19, col = "blue")
  lines(id_data$Year, id_data$Value, col = "black", type = "l")
  
  print(unique(id_data$Binomial))
  readline(prompt = "Press [Enter] to continue")
}

#Data to remove after secondary visual inspection after detrending
df_lpi_ts_detrend_rm <- df_lpi_ts_detrend %>%
  filter(! (ID == 4860 | ID == 7673 | ID == 10096 |
              ID == 19783))

df_lpi_ts_detrended_rm <- detrend_results %>%
  filter(! (ID == 4860 | ID == 7673 | ID == 10096 |
            ID == 19783))

#Add back taxonomic Data
df_lpi_ts_detrend_sub <- df_lpi_ts_detrend_rm %>% dplyr::select(-Value, -Year, -ID)

df_lpi_detrended <- cbind(df_lpi_ts_detrend_sub, df_lpi_ts_detrended_rm)

df_lpi_detrended_notrend <- rbind(df_lpi_notrend, df_lpi_detrended)

unique(df_lpi_detrended_notrend$ID)
unique(df_lpi_detrended_notrend$Binomial)

write.csv(df_lpi_detrended_notrend, "../Slow_Fast_IS_Stability/R/Outputs/TS_lpi_detrended.csv")

