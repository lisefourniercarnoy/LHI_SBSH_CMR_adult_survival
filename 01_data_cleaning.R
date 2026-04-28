###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2010-2025 CMR data

# script objective: clean up raw data for analysis.

###############################################################################

rm(list = ls())

library(tidyverse) # for data manipulation
library(stringr) # for fixing capture histories
library(RMark) # for running rmark
library(zoo)

## Read in raw CMR data -------------------------------------------------------

dat <- read.csv("data/raw/LHI_FFSH_clean_live_recap.csv", skip = 1) %>% 
  mutate(across(where(is.integer), ~ ifelse(is.na(.), 0, 1))) %>% # fill NAs with 0 and non-NAs with 1
  rename_with(~ gsub("^X", "", .), starts_with("X")) %>% 
  mutate(`2019` = ifelse(is.na(`2019`), 0, `2019`)) %>% # no recaptures for 2019 breeding season, fill with 0
  glimpse()

dat %>%
  select(`2009`:`2024`) %>%
  # for each bird, find the first year they were marked (first column with a 1)
  mutate(first_year = names(.)[max.col(. != 0, ties.method = "first")]) %>%
  count(first_year) %>%
  mutate(first_year = as.integer(first_year)) %>%
  ggplot(aes(x = first_year, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = 2009:2024) +
  labs(
    y = "Number of birds marked"
  ) +
  theme_classic()

# remove birds from early years. (there are few birds in these years, so estimating survival in these years is difficult numerically and not very informative)
dat <- dplyr::select(dat, -c(`2009`:`2012`))

# add a 'capture history' (ch) column
dat <- dat %>%
  mutate(ch = apply(select(., first(starts_with("20")):`2024`), 1, paste0, collapse = "")) %>%
  glimpse()


# format chick data so that you're culling chick years
dat <- dat %>%
  mutate(ch = ifelse(Age == "P", 
                     sub("^.*?1.*?(1.*)$", "\\1", ch), 
                     ch)) %>%
  glimpse()

# fill the culled capture histories so everything's the right format
max_length <- max(nchar(dat$ch))  # Find the longest string length

dat <- dat %>%
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>% # pad the character string with zeroes on the left
  glimpse()

# age class '?' are '1+' birds, so technically adults
dat$Age <- ifelse(dat$Age == "1+", "A", dat$Age)

# remove transients (birds banded as chicks that were only seen once- adults who are only seen once are not considered as transients)
dat <- dat[dat$Age != "P" | str_count(dat$ch, "1") > 1, ]
dat <- dat[str_count(dat$ch, "1") != 0, ]


glimpse(dat)
colSums(dat[3:14]) # make sure the first few years actually have birds in them, otherwise remove those years (line 26)

# save the cleaned data
saveRDS(dat, "data/tidy/LHI_FFSH_capture_histories.rds")


## Environmental data ---------------------------------------------------------

# we want to understand whether/what environmental index has an influence over survival.
# the wintering period (May-September) is what we will test. We'll use a 3-month average of each variable (July-September)


### El Nino Southern Oscillation (ENSO) ---------------------------------------

enso <- read.csv("data/raw/ENSO values.csv", skip = 1) %>%  # skip the first row (data source info)
  rename(# rename for clarity
         ENSO_jun_jul_aug = JJA,
         ENSO_jul_aug_sep = JAS,
         #ENSO_aug_sep_oct = ASO
         ) %>% 
  glimpse()

# Select the desired columns
enso <- enso[, c("time", 
                 "ENSO_jun_jul_aug", 
                 "ENSO_jul_aug_sep" 
                 #"ENSO_aug_sep_oct"
                 )]
enso <- enso[enso$time >= 2013 & enso$time <= 2024,]

print(enso) # check the format - years have to be called 'time' for rmark to recognise it.

# save the tidy environmental dat
saveRDS(enso, "data/tidy/ENSO_tidy.rds")


### Pacific Decadal Oscillation -----------------------------------------------

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
    PDO_jun_jul_aug = rowMeans(select(., Jun, Jul, Aug), na.rm = TRUE),
    PDO_jul_aug_sep = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE)
  ) %>%
  dplyr::select(time,
                PDO_jun_jul_aug, 
                PDO_jul_aug_sep
                ) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

print(pdo_clean) # check the format - years have to be called 'time' for rmark to recognise it. 

# save the tidy environmental dat
saveRDS(pdo_clean, "data/tidy/PDO_tidy.rds")


### Antarctic Oscillation -----------------------------------------------------

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
    AAO_jun_jul_aug = rowMeans(select(., Jun, Jul, Aug), na.rm = TRUE),
    AAO_jul_aug_sep = rowMeans(select(., Jul, Aug, Sep), na.rm = TRUE)
  ) %>%
  dplyr::select(time, 
                AAO_jun_jul_aug, 
                AAO_jul_aug_sep
                ) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

print(aao_clean) # resulting data has breeding season PDO values, non-breeding season, as well as 1 year lag of each

# save the tidy environmental data
saveRDS(aao_clean, "data/tidy/AAO_tidy.rds")


### Arctic Oscillation --------------------------------------------------------

AO <- read.csv("data/raw/AO values.csv", skip = 1) %>%  # skip the first row (data source info)
  glimpse()

# Select the desired columns
ao_clean <- AO %>%
  mutate(
    AAO_jun_jul_aug = rowMeans(select(., June, July, August), na.rm = TRUE),
    AO_jul_aug_sep = rowMeans(select(., July, August, September), na.rm = TRUE)
  ) %>%
  dplyr::select(time, 
                AAO_jun_jul_aug, 
                AO_jul_aug_sep
  ) %>% 
  dplyr::filter(time >= 2013 & time <= 2024) %>% 
  glimpse()

# save the tidy environmental dat
saveRDS(ao_clean, "data/tidy/AO_tidy.rds")


# SST data --------------------------------------------------------------------

## sst is the only environmental thing that is missing data to make an 'in-between breeding seasons' metric (july-sept.) so we are lumping june-oct for this one only.

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
unique(sst$datetime)


# Extract year and month
sst_temp <- sst %>%
  mutate(
    time = year(datetime),
    month = month(datetime)
  ) %>%
  filter(month %in% 6:10) %>%
  group_by(time) %>%
  summarise(
    temp_jun_oct = mean(temp_10_15m[month %in% c(6,7,8,9,10)], na.rm = TRUE),
    
    #temp_jun_jul_aug = mean(temp_10_15m[month %in% c(6,7,8)], na.rm = TRUE),
    #temp_jul_aug_sep = mean(temp_10_15m[month %in% c(7,8,9)], na.rm = TRUE),
    #temp_aug_sep_oct = mean(temp_10_15m[month %in% c(8,9,10)], na.rm = TRUE),
    .groups = "drop"
  )

# Remove rows with NA if needed (i.e., first year)
# sst_final <- sst_final %>%
#   filter(!is.na(temp_prev_may_jun_jul) & !is.na(temp_prev_sep_oct_nov))
sst_final <- sst_temp %>% dplyr::filter(time >= 2013)
print(sst_final)

saveRDS(sst_final, "data/tidy/temp_tidy.rds")


# Wintering temperature -------------------------------------------------------

sst_w <- read.csv("data/raw/SST japan values.csv", skip = 1) %>%
  glimpse()

sst_w <- data.frame(time = sst_w[,1], sst_w_mean = rowMeans(sst_w[,-1]))

plot(sst_w$sst_w_mean ~ sst_w$time)

saveRDS(sst_w, "data/tidy/temp_w_tidy.rds")



### END ###