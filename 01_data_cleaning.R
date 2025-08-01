###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2010-2025 CMR data

# script objective: clean up raw data for analysis

###############################################################################

rm(list = ls())

library(tidyverse) # for data manipulation
library(stringr) # for fixing capture histories
library(RMark) # for running rmark

## Read in raw CMR data -------------------------------------------------------

dat <- read.csv("data/raw/LHI_FFSH_clean_live_recap.csv", skip = 1) %>% 
  mutate(across(where(is.integer), ~ ifelse(is.na(.), 0, 1))) %>% # fill NAs with 0 and non-NAs with 1
  rename_with(~ gsub("^X", "", .), starts_with("X")) %>% 
  mutate(`2019` = ifelse(is.na(`2019`), 0, `2019`)) %>% # no recaptures for 2019 breeding season, fill with 0
  glimpse()

# remove 2009-banded birds because there is no env data
dat <- dplyr::select(dat, -`2009`)

# add a 'capture history' (ch) column
dat <- dat %>%
  mutate(ch = apply(select(., `2010`:`2024`), 1, paste0, collapse = "")) %>%
  glimpse()

# remove transients (birds banded as chicks that were only seen once- adults who are only seen once are not considered as transients)
dat <- dat[dat$Age != "P" | str_count(dat$ch, "1") > 1, ]


# format chick data so that you're culling chick years
dat <- dat %>%
  mutate(ch = ifelse(Age == "P", 
                     sub("^.*?1.*?(1.*)$", "\\1", ch), 
                     ch)) %>%
  glimpse()

# fill the culled capture histories so everything's the right format
max_length <- max(nchar(dat$ch))  # Find the longest string length

dat <- dat %>%
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>%
  glimpse()

# select relevant columns
dat <- dat %>%
  dplyr::select(Band, ch, Age) %>% 
  glimpse()

# age class '?' are '1+' birds, so technically adults
dat$Age <- ifelse(dat$Age == "1+", "A", dat$Age)

# split the capture histories again (needed for GOF tests)
dat <- dat %>%
  mutate(ch_split = strsplit(ch, "")) %>%
  unnest_wider(ch_split, names_sep = "_") %>%
  rename_with(~ as.character(2010:2024), starts_with("ch_split_")) %>% 
  glimpse(dat)

# save the cleaned data
saveRDS(dat, "data/tidy/LHI_FFSH_capture_histories.rds")


## Read in raw CMR data (known fate, archive) ---------------------------------

dat <- read.csv("data/raw/LHI_FFSH_clean_live_recap.csv") %>% 
  mutate(across(where(is.integer), ~ ifelse(is.na(.), 0, 1))) %>% # fill NAs with 0 and non-NAs with 1
  rename_with(~ gsub("^X", "", .), starts_with("X")) %>% 
  mutate(`2019` = ifelse(is.na(`2019`), 0, `2019`)) %>% # no recaptures for 2019 breeding season, fill with 0
  glimpse()

# remove 2009-banded birds because there is no env data
dat <- dplyr::select(dat, -`2009`)

# remove transients
date_cols <- as.character(2010:2024)

dat <- dat %>%
  filter(rowSums(select(., all_of(date_cols))) > 1)
# dat <- dat %>%
#   filter(str_count(ch, "1") > 1)

# add in the known-deceased birds 
# dead_lookup <- read_csv("data/raw/LHI_FFSH_known_deceased.csv") %>%
#   transmute(Band = band,
#             year = as.character(breeding_season_recovered),
#             dead = 1)
# 
# dead_matrix <- expand.grid(Band = dat$Band,
#                            year = names(dat)[-(1:2)],
#                            stringsAsFactors = FALSE) %>%
#   left_join(dead_lookup, by = c("Band", "year")) %>%
#   mutate(dead = replace_na(dead, 0)) %>%
#   pivot_wider(names_from = year,
#               values_from = dead,
#               names_prefix = "",
#               names_sep = "_dead") %>%
#   rename_with(~ paste0(gsub("_dead", "", .), "_dead"), -Band)
# dat <- bind_cols(dat, dead_matrix %>% select(-Band))
# 
# # reorder the dataframe columns
# years <- grep("^\\d{4}$", names(dat), value = TRUE)
# interleaved_cols <- unlist(lapply(years, function(y) c(y, paste0(y, "_dead"))))
# dat_reordered <- dat %>% 
#   select(Age, Band, all_of(interleaved_cols))

# # add a 'capture history' (ch) column
# dat <- dat_reordered %>%
#   mutate(ch = apply(select(., `2010`:`2024_dead`), 1, paste0, collapse = "")) %>%
#   glimpse()



# format chick data so that you're culling chick years
dat <- dat %>%
  mutate(ch = ifelse(Age == "P", 
                     sub("^.*?1.*?(1.*)$", "\\1", ch), 
                     ch)) %>%
  glimpse()

# fill the culled capture histories so everything's the right format
max_length <- max(nchar(dat$ch))  # Find the longest string length

dat <- dat %>%
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>%
  glimpse()

# select relevant columns
dat <- dat %>%
  dplyr::select(Band, ch, Age) %>% 
  glimpse()


# save the cleaned data
saveRDS(dat, "data/tidy/LHI_FFSH_capture_histories_known_fate.rds")


## Environmental data ---------------------------------------------------------

enso <- read.csv("data/raw/ENSO values.csv", skip = 1) %>%  # skip the first row (data source info)
  rename(ENSO_sep_oct_nov = SON, # rename for clarity
         ENSO_may_jun_jul = MJJ
         ) %>% 
  glimpse()

# 1-year lag
enso$ENSO_prev_yr_sep_oct_nov <- dplyr::lag(enso$ENSO_sep_oct_nov, 1)
enso$ENSO_prev_yr_may_jun_jul <- dplyr::lag(enso$ENSO_may_jun_jul, 1)

# Select the desired columns
enso <- enso[, c("time", "ENSO_may_jun_jul", "ENSO_sep_oct_nov", "ENSO_prev_yr_may_jun_jul", "ENSO_prev_yr_sep_oct_nov")]
enso <- enso[enso$time >= 2009 & enso$time <= 2024,]
# View result
print(enso)

# save the tidy environmental dat
saveRDS(enso, "data/tidy/ENSO_tidy.rds")


# Pacific Decadal Oscillation -------------------------------------------------

pdo <- readLines("data/raw/PDO values.txt")
pdo_data <- pdo[-c(1,2)] # Remove the first two header lines

pdo_split <- strsplit(pdo_data, "\\s+") # Split each line by whitespace
pdo_matrix <- do.call(rbind, lapply(pdo_split, function(x) x[nzchar(x)]))
pdo <- as.data.frame(pdo_matrix, stringsAsFactors = FALSE) # Convert to data frame 
colnames(pdo) <- c("Year", month.abb) # name columns
pdo[] <- lapply(pdo, function(x) ifelse(suppressWarnings(!is.na(as.numeric(x))), as.numeric(x), x)) # Convert numeric columns

pdo_clean <- pdo %>%
  mutate(
    time = Year,
    PDO_sep_oct_nov = rowMeans(select(., Sep, Oct, Nov), na.rm = TRUE),
    PDO_may_jun_jul = rowMeans(select(., May, Jun, Jul), na.rm = TRUE),
    PDO_prev_yr_sep_oct_nov = lag(PDO_sep_oct_nov),
    PDO_prev_yr_may_jun_jul = lag(PDO_may_jun_jul)
  ) %>%
  dplyr::select(time, PDO_sep_oct_nov, PDO_may_jun_jul, PDO_prev_yr_sep_oct_nov, PDO_prev_yr_may_jun_jul) %>% 
  dplyr::filter(time >= 2009 & time <= 2024) %>% 
  glimpse()

print(pdo_clean) # resulting data has breeding season PDO values, non-breeding season, as well as 1 year lag of each

# save the tidy environmental dat
saveRDS(pdo_clean, "data/tidy/PDO_tidy.rds")


# Antarctic Oscillation -------------------------------------------------------

aao <- readLines("data/raw/AAO values.txt")
aao_data <- aao[-c(1, 2)] # Remove the first two header lines

aao_split <- strsplit(aao_data, "\\s+") # Split each line by whitespace
aao_matrix <- do.call(rbind, lapply(aao_split, function(x) x[nzchar(x)]))
aao <- as.data.frame(aao_matrix, stringsAsFactors = FALSE) # Convert to data frame 
colnames(aao) <- c("Year", month.abb) # name columns
aao[] <- lapply(aao, function(x) ifelse(suppressWarnings(!is.na(as.numeric(x))), as.numeric(x), x)) # Convert numeric columns

aao_clean <- aao %>%
  mutate(
    time = Year,
    AAO_sep_oct_nov = rowMeans(select(., Sep, Oct, Nov), na.rm = TRUE),
    AAO_may_jun_jul = rowMeans(select(., May, Jun, Jul), na.rm = TRUE),
    AAO_prev_yr_sep_oct_nov = lag(AAO_sep_oct_nov),
    AAO_prev_yr_may_jun_jul = lag(AAO_may_jun_jul)
  ) %>%
  dplyr::select(time, AAO_sep_oct_nov, AAO_may_jun_jul, AAO_prev_yr_sep_oct_nov, AAO_prev_yr_may_jun_jul) %>% 
  dplyr::filter(time >= 2009 & time <= 2024) %>% 
  glimpse()

print(aao_clean) # resulting data has breeding season PDO values, non-breeding season, as well as 1 year lag of each

# save the tidy environmental dat
saveRDS(aao_clean, "data/tidy/AAO_tidy.rds")


# SST data --------------------------------------------------------------------

sst <- read.csv("data/raw/IMOS_-_Australian_National_Mooring_Network_(ANMN)_-_CTD_Profiles.csv", skip = 29) %>%
    mutate(datetime = ymd_hms(TIME)) %>% 
  dplyr::select(c(datetime, 
                  LATITUDE, LONGITUDE, 
                  DEPTH,
                  TEMP, TEMP_quality_control)) %>% # select relevant columns
  dplyr::filter(DEPTH > 10 & DEPTH < 15) %>% 
  dplyr::filter(TEMP_quality_control == 1 | TEMP_quality_control == 2) %>%  # quality control codes 1 = good data, 2 = prbably good data, 4 is bad data
  group_by(datetime) %>% 
  mutate(temp_10_15m = mean(TEMP)) %>% 
  glimpse()

# Extract year and month
sst_temp <- sst %>%
  mutate(
    time = year(datetime), # has to be 'time' for RMark
    month = month(datetime)
  ) %>%
  filter(month %in% c(5, 6, 7, 9, 10, 11)) %>% # filter for only months we care about
  mutate(
    season = case_when(
      month %in% c(5, 6, 7) ~ "may_jun_jul",
      month %in% c(9, 10, 11) ~ "sep_oct_nov"
    )
  ) %>%
  group_by(time, season) %>%
  summarise(
    mean_temp = mean(temp_10_15m, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  glimpse()

# Pivot to wide format
sst_wide <- sst_temp %>%
  pivot_wider(
    names_from = season,
    values_from = mean_temp,
    names_prefix = "temp_"
  ) %>% 
  glimpse()

# Create lagged (previous year) temperature columns
sst_final <- sst_wide %>%
  arrange(time) %>%
  mutate(
    temp_prev_may_jun_jul = lag(temp_may_jun_jul),
    temp_prev_sep_oct_nov = lag(temp_sep_oct_nov)
  ) %>% 
  glimpse()

# Remove rows with NA if needed (i.e., first year)
# sst_final <- sst_final %>%
#   filter(!is.na(temp_prev_may_jun_jul) & !is.na(temp_prev_sep_oct_nov))
print(sst_final)

saveRDS(sst_final, "data/tidy/temp_tidy.rds")

### END ###