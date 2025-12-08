#Checking LPI taxonomic names

#Author(s): Reilly O'Connor
#Version: 2024-03-03

#Bold and taxize have been discontinued and need to be installed manually
#install.packages("remotes")
#remotes::install_github("ropensci/bold")
#remotes::install_github("ropensci/taxize")

#Pkgs
library(tidyverse)
library(RColorBrewer)
library(beepr)
library(taxize)

#load data
df_lpi <- read.csv("../Slow_Fast_IS_Stability/R/Outputs/TS_lpi_detrended.csv")

##### Code #####
unique(df_lpi$ID)
unique(df_lpi$Binomial)

#Taxonomic Name Corrections
df_lpi_sub <- df_lpi %>%
  mutate(Binomial = case_when(
    Binomial == "Arctocephalus_townsendi" ~ "Arctocephalus_philippii",
    Binomial == "Spermophilus_franklinii" ~ "Poliocitellus_franklinii",
    Binomial == "Spermophilus_tridecemlineatus" ~ "Ictidomys_tridecemlineatus",
    Binomial == "Erethizon_dorsata" ~ "Erethizon_dorsatus",
    Binomial == "Equus_burchellii" ~ "Equus_quagga",
    Binomial == "Tragelaphus_oryx" ~ "Taurotragus_oryx",
    Binomial == "Spermophilus_parryii" ~ "Urocitellus_parryii",
    Binomial == "Dicrostonyx_vinogradovi" ~ "Dicrostonyx_torquatus",
    Binomial == "Monachus_schauinslandi" ~ "Neomonachus_schauinslandi",
    Binomial == "Clethrionomys_gapperi" ~ "Myodes_gapperi",
    Binomial == "Spermophilus_columbianus" ~ "Urocitellus_columbianus",
    Binomial == "Neovison_vison" ~ "Mustela_vison",
    Binomial == "Alexandromys_oeconomus" ~ "Microtus_oeconomus",
    Binomial == "Panthera_uncia" ~ "Uncia_uncia",
    Binomial == "Lasiurus_cinereus" ~ "Aeorestes_cinereus",
    Binomial == "Otaria_flavescens" ~ "Otaria_byronia",
    Binomial == "Lama_guanicoe" ~ "Lama_glama",
    Binomial == "Procyon_gymnocercus" ~ "Lycalopex_gymnocercus",
    Binomial == "Martes_pennanti" ~ "Pekania_pennanti",
    Binomial == "Colobus_tephrosceles" ~ "Piliocolobus_tephrosceles",
    Binomial == "Physeter_catodon" ~ "Physeter_macrocephalus",
    Binomial == "Sorex_bendirri" ~ "Sorex_bendirii",
    Binomial == "Delphinus_capensis" ~ "Delphinus_delphis",
    TRUE ~ Binomial))


df_lpi_sub$species_names <- gsub("_", " ", df_lpi_sub$Binomial)

species_names <- unique(df_lpi_sub$species_names)

lpi_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

gbif_ids <- lpi_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(gbif_ids, "match")
names <- attr(gbif_ids, "names")

# Combine into a data frame
df_lpi_gbif_ids <- data.frame(
  gbif_id = gbif_ids,
  species_names = names,
  stringsAsFactors = FALSE)

df_lpi_gbif_ids <- df_lpi_gbif_ids %>% select(gbif_id.ids, species_names) %>%
  rename(GBIF_ID = gbif_id.ids)

df_lpi_ids_final <- merge(df_lpi_sub, df_lpi_gbif_ids, by = "species_names")
unique(df_lpi_ids_final$Binomial)

write.csv(df_lpi_ids_final, "../Slow_Fast_IS_Stability/R/Outputs/TS_lpi_detrended resolved taxonomy.csv")

##### Hatton Dataframes #####
df_metab <- read.csv("../Slow_Fast_IS_Stability/R/Data/Hatton et al 2019/slow_fast_metabolism.csv", header = T)
df_pop_metab <- read.csv("../Slow_Fast_IS_Stability/R/Data/Hatton et al 2019/pop_metab.csv", header = T)
df_growth <- read.csv("../Slow_Fast_IS_Stability/R/Data/Hatton et al 2019/slow_fast_growth.csv", header = T)

#metabolism
df_metab_sub <- df_metab %>%
  filter(Major_taxa == "Mammal") %>%
  mutate(Species = case_when(
    Species == "Arctocephalus townsendi" ~ "Arctocephalus philippii",
    Species == "Spermophilus franklinii" ~ "Poliocitellus franklinii",
    Species == "Spermophilus tridecemlineatus" ~ "Ictidomys tridecemlineatus",
    Species == "Erethizon dorsata" ~ "Erethizon dorsatus",
    Species == "Equus burchellii" ~ "Equus quagga",
    Species == "Tragelaphus oryx" ~ "Taurotragus oryx",
    Species == "Spermophilus parryii" ~ "Urocitellus parryii",
    Species == "Dicrostonyx vinogradovi" ~ "Dicrostonyx torquatus",
    Species == "Monachus schauinslandi" ~ "Neomonachus schauinslandi",
    Species == "Clethrionomys gapperi" ~ "Myodes gapperi",
    Species == "Spermophilus columbianus" ~ "Urocitellus columbianus",
    Species == "Neovison vison" ~ "Mustela vison",
    Species == "Alexandromys oeconomus" ~ "Microtus oeconomus",
    Species == "Panthera uncia" ~ "Uncia uncia",
    Species == "Lasiurus cinereus" ~ "Aeorestes cinereus",
    Species == "Otaria flavescens" ~ "Otaria byronia",
    Species == "Lama guanicoe" ~ "Lama glama",
    Species == "Procyon gymnocercus" ~ "Lycalopex gymnocercus",
    Species == "Martes pennanti" ~ "Pekania pennanti",
    Species == "Colobus tephrosceles" ~ "Piliocolobus tephrosceles",
    Species == "Physeter catodon" ~ "Physeter macrocephalus",
    Species == "Sorex bendirri" ~ "Sorex bendirii",
    Species == "Delphinus capensis" ~ "Delphinus delphis",
    TRUE ~ Species))


species_names <- unique(df_metab_sub$Species)
metab_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

metab_gbif_ids <- metab_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(metab_gbif_ids, "match")
names <- attr(metab_gbif_ids, "names")

# Combine into a data frame
df_metab_gbif_ids <- data.frame(
  gbif_id = metab_gbif_ids,
  Species = names,
  stringsAsFactors = FALSE)

df_metab_gbif_ids <- df_metab_gbif_ids %>% select(gbif_id.ids, Species) %>%
  rename(GBIF_ID = gbif_id.ids)

df_metab_gbif_final <- merge(df_metab_sub, df_metab_gbif_ids, by = "Species")

write.csv(df_metab_gbif_final, "../Slow_Fast_IS_Stability/R/Outputs/slow_fast_metabolism resolved taxonomy.csv")

#rmax
df_growth_sub <- df_growth %>%
  filter(Major_taxa == "Mammal")%>%
  mutate(Species = case_when(
    Species == "Arctocephalus townsendi" ~ "Arctocephalus philippii",
    Species == "Spermophilus franklinii" ~ "Poliocitellus franklinii",
    Species == "Spermophilus tridecemlineatus" ~ "Ictidomys tridecemlineatus",
    Species == "Erethizon dorsata" ~ "Erethizon dorsatus",
    Species == "Equus burchellii" ~ "Equus quagga",
    Species == "Tragelaphus oryx" ~ "Taurotragus oryx",
    Species == "Spermophilus parryii" ~ "Urocitellus parryii",
    Species == "Dicrostonyx vinogradovi" ~ "Dicrostonyx torquatus",
    Species == "Monachus schauinslandi" ~ "Neomonachus schauinslandi",
    Species == "Clethrionomys gapperi" ~ "Myodes gapperi",
    Species == "Spermophilus columbianus" ~ "Urocitellus columbianus",
    Species == "Neovison vison" ~ "Mustela vison",
    Species == "Alexandromys oeconomus" ~ "Microtus oeconomus",
    Species == "Panthera uncia" ~ "Uncia uncia",
    Species == "Lasiurus cinereus" ~ "Aeorestes cinereus",
    Species == "Otaria flavescens" ~ "Otaria byronia",
    Species == "Lama guanicoe" ~ "Lama glama",
    Species == "Procyon gymnocercus" ~ "Lycalopex gymnocercus",
    Species == "Martes pennanti" ~ "Pekania pennanti",
    Species == "Colobus tephrosceles" ~ "Piliocolobus tephrosceles",
    Species == "Physeter catodon" ~ "Physeter macrocephalus",
    Species == "Sorex bendirri" ~ "Sorex bendirii",
    Species == "Delphinus capensis" ~ "Delphinus delphis",
    TRUE ~ Species))

species_names <- unique(df_growth_sub$Species)
grw_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

grw_gbif_ids <- grw_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(grw_gbif_ids, "match")
names <- attr(grw_gbif_ids, "names")

# Combine into a data frame
df_grw_gbif_ids <- data.frame(
  gbif_id = grw_gbif_ids,
  Species = names,
  stringsAsFactors = FALSE)

df_grw_gbif_ids <- df_grw_gbif_ids %>% select(gbif_id.ids, Species) %>%
  rename(GBIF_ID = gbif_id.ids)

df_grw_gbif_final <- merge(df_growth_sub, df_grw_gbif_ids, by = "Species")

write.csv(df_grw_gbif_final, "../Slow_Fast_IS_Stability/R/Outputs/slow_fast_growth resolved taxonomy.csv")

#rmax subset - count method rm only
Growth <- read.csv("../Slow_Fast_IS_Stability/R/Data/Hatton et al 2019/Growth.csv", header = T)

mammal_growth <- Growth %>%
  filter(Major_taxa == "Mammal")

mammal_growth_filtered <- mammal_growth %>%
  filter(grepl("method: count", Reference) | 
           grepl("Global Population Dynamics Database", Reference))
#rmax
df_mammal_growth_filtered_sub <- mammal_growth_filtered %>%
  filter(Major_taxa == "Mammal")%>%
  mutate(Species = case_when(
    Species == "Arctocephalus townsendi" ~ "Arctocephalus philippii",
    Species == "Spermophilus franklinii" ~ "Poliocitellus franklinii",
    Species == "Spermophilus tridecemlineatus" ~ "Ictidomys tridecemlineatus",
    Species == "Erethizon dorsata" ~ "Erethizon dorsatus",
    Species == "Equus burchellii" ~ "Equus quagga",
    Species == "Tragelaphus oryx" ~ "Taurotragus oryx",
    Species == "Spermophilus parryii" ~ "Urocitellus parryii",
    Species == "Dicrostonyx vinogradovi" ~ "Dicrostonyx torquatus",
    Species == "Monachus schauinslandi" ~ "Neomonachus schauinslandi",
    Species == "Clethrionomys gapperi" ~ "Myodes gapperi",
    Species == "Spermophilus columbianus" ~ "Urocitellus columbianus",
    Species == "Neovison vison" ~ "Mustela vison",
    Species == "Alexandromys oeconomus" ~ "Microtus oeconomus",
    Species == "Panthera uncia" ~ "Uncia uncia",
    Species == "Lasiurus cinereus" ~ "Aeorestes cinereus",
    Species == "Otaria flavescens" ~ "Otaria byronia",
    Species == "Lama guanicoe" ~ "Lama glama",
    Species == "Procyon gymnocercus" ~ "Lycalopex gymnocercus",
    Species == "Martes pennanti" ~ "Pekania pennanti",
    Species == "Colobus tephrosceles" ~ "Piliocolobus tephrosceles",
    Species == "Physeter catodon" ~ "Physeter macrocephalus",
    Species == "Sorex bendirri" ~ "Sorex bendirii",
    Species == "Delphinus capensis" ~ "Delphinus delphis",
    TRUE ~ Species))

species_names <- unique(df_mammal_growth_filtered_sub$Species)
grw_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

grw_gbif_ids <- grw_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(grw_gbif_ids, "match")
names <- attr(grw_gbif_ids, "names")

# Combine into a data frame
df_grw_gbif_ids <- data.frame(
  gbif_id = grw_gbif_ids,
  Species = names,
  stringsAsFactors = FALSE)

df_grw_gbif_ids <- df_grw_gbif_ids %>% select(gbif_id.ids, Species) %>%
  rename(GBIF_ID = gbif_id.ids)

df_grw_gbif_final <- merge(df_mammal_growth_filtered_sub, df_grw_gbif_ids, by = "Species")

write.csv(df_grw_gbif_final, "../Slow_Fast_IS_Stability/R/Outputs/slow_fast_growth resolved taxonomy - count method.csv")


##### Mhyrvold Database #####
df_body_sizes <- read.csv("../Slow_Fast_IS_Stability/R/Data/Amniote Life History Database/Data_Files/Amniote_Database_Aug_2015.csv")

df_body_sizes_bs <- df_body_sizes %>% mutate(Species = paste(genus, species)) %>%
  rename(Mass_g = adult_body_mass_g) %>%
  filter(class == "Mammalia") %>%
  select(Species, Mass_g, litter_or_clutch_size_n, litters_or_clutches_per_y, longevity_y, female_maturity_d, gestation_d, weaning_d) %>%
  filter(Mass_g > 0)

df_body_sizes_bs$Binomial <- gsub(" ", "_", df_body_sizes_bs$Species)

df_bs_sub <- df_body_sizes_bs%>%
  mutate(Species = case_when(
    Species == "Arctocephalus townsendi" ~ "Arctocephalus philippii",
    Species == "Spermophilus franklinii" ~ "Poliocitellus franklinii",
    Species == "Spermophilus tridecemlineatus" ~ "Ictidomys tridecemlineatus",
    Species == "Erethizon dorsata" ~ "Erethizon dorsatus",
    Species == "Equus burchellii" ~ "Equus quagga",
    Species == "Tragelaphus oryx" ~ "Taurotragus oryx",
    Species == "Spermophilus parryii" ~ "Urocitellus parryii",
    Species == "Dicrostonyx vinogradovi" ~ "Dicrostonyx torquatus",
    Species == "Monachus schauinslandi" ~ "Neomonachus schauinslandi",
    Species == "Clethrionomys gapperi" ~ "Myodes gapperi",
    Species == "Spermophilus columbianus" ~ "Urocitellus columbianus",
    Species == "Neovison vison" ~ "Mustela vison",
    Species == "Alexandromys oeconomus" ~ "Microtus oeconomus",
    Species == "Panthera uncia" ~ "Uncia uncia",
    Species == "Lasiurus cinereus" ~ "Aeorestes cinereus",
    Species == "Otaria flavescens" ~ "Otaria byronia",
    Species == "Lama guanicoe" ~ "Lama glama",
    Species == "Procyon gymnocercus" ~ "Lycalopex gymnocercus",
    Species == "Martes pennanti" ~ "Pekania pennanti",
    Species == "Colobus tephrosceles" ~ "Piliocolobus tephrosceles",
    Species == "Physeter catodon" ~ "Physeter macrocephalus",
    Species == "Sorex bendirri" ~ "Sorex bendirii",
    Species == "Delphinus capensis" ~ "Delphinus delphis",
    TRUE ~ Species))

species_names <- unique(df_bs_sub$Species)
bs_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

bs_gbif_ids <- bs_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(bs_gbif_ids, "match")
names <- attr(bs_gbif_ids, "names")

# Combine into a data frame
df_bs_gbif_ids <- data.frame(
  gbif_id = bs_gbif_ids,
  Species = names,
  stringsAsFactors = FALSE)

df_bs_gbif_ids <- df_bs_gbif_ids %>% select(gbif_id.ids, Species) %>%
  rename(GBIF_ID = gbif_id.ids)

df_bs_gbif_final <- merge(df_bs_sub, df_bs_gbif_ids, by = "Species")

write.csv(df_bs_gbif_final, "../Slow_Fast_IS_Stability/R/Outputs/Amniote database resolved taxonomy.csv")

##### Beccari et al 2024 Data #####
#Beccari et al 2024
df_fast_slow <- read.table("../Slow_Fast_IS_Stability/R/Data/Beccari et al 2024/Mammal_BMcor_Def copy.txt", header = T, sep = "", dec = ".")

df_fast_slow <- df_fast_slow %>%
  rownames_to_column(var = "Species")

df_fast_slow$Binomial <- gsub(" ", "_", df_fast_slow$Species)

df_fs_sub <- df_fast_slow %>%
  mutate(Species = case_when(
    Species == "Arctocephalus townsendi" ~ "Arctocephalus philippii",
    Species == "Spermophilus franklinii" ~ "Poliocitellus franklinii",
    Species == "Spermophilus tridecemlineatus" ~ "Ictidomys tridecemlineatus",
    Species == "Erethizon dorsata" ~ "Erethizon dorsatus",
    Species == "Equus burchellii" ~ "Equus quagga",
    Species == "Tragelaphus oryx" ~ "Taurotragus oryx",
    Species == "Spermophilus parryii" ~ "Urocitellus parryii",
    Species == "Dicrostonyx vinogradovi" ~ "Dicrostonyx torquatus",
    Species == "Monachus schauinslandi" ~ "Neomonachus schauinslandi",
    Species == "Clethrionomys gapperi" ~ "Myodes gapperi",
    Species == "Spermophilus columbianus" ~ "Urocitellus columbianus",
    Species == "Neovison vison" ~ "Mustela vison",
    Species == "Alexandromys oeconomus" ~ "Microtus oeconomus",
    Species == "Panthera uncia" ~ "Uncia uncia",
    Species == "Lasiurus cinereus" ~ "Aeorestes cinereus",
    Species == "Otaria flavescens" ~ "Otaria byronia",
    Species == "Lama guanicoe" ~ "Lama glama",
    Species == "Procyon gymnocercus" ~ "Lycalopex gymnocercus",
    Species == "Martes pennanti" ~ "Pekania pennanti",
    Species == "Colobus tephrosceles" ~ "Piliocolobus tephrosceles",
    Species == "Physeter catodon" ~ "Physeter macrocephalus",
    Species == "Sorex bendirri" ~ "Sorex bendirii",
    Species == "Delphinus capensis" ~ "Delphinus delphis",
    TRUE ~ Species))

species_names <- unique(df_fs_sub$Species)
fs_gbif_ids <- get_ids(species_names, db = "gbif", rows=1, messages = TRUE, searchtype = "scientific", accepted = TRUE)

fs_gbif_ids <- fs_gbif_ids$gbif

#Extract the species names and ID #s
matches <- attr(fs_gbif_ids, "match")
names <- attr(fs_gbif_ids, "names")

#Combine into a data frame
df_fs_gbif_ids <- data.frame(
  gbif_id = fs_gbif_ids,
  Species = names,
  stringsAsFactors = FALSE)

df_fs_gbif_ids <- df_fs_gbif_ids %>% select(gbif_id.ids, Species) %>%
  rename(GBIF_ID = gbif_id.ids)

df_fs_gbif_final <- merge(df_fs_sub, df_fs_gbif_ids, by = "Species")


write.csv(df_fs_gbif_final, "../Slow_Fast_IS_Stability/R/Outputs/Beccari Fast-Slow resolved taxonomy.csv")



