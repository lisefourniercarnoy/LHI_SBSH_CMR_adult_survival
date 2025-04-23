###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2011-2024 CMR data

# script objective: clean up raw data for analysis

###############################################################################

library(tidyverse)
library(stringr)
library(RMark)


## Read in raw CMR data -------------------------------------------------------

dat <- read.csv("data/raw/FFSH encounters 20250311.csv") %>% glimpse()
summary(as.factor(dat$Age))

# add a 'capture history' (ch) column
dat <- dat %>%
  mutate(ch = apply(select(., X2010:X2024), 1, paste0, collapse = "")) %>% 
  glimpse()

# remove transients
dat <- dat %>% 
  dplyr::filter(Include.in.analysis == "TRUE") %>% 
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
  mutate(ch = str_pad(ch, width = max_length, side = "left", pad = "0")) %>%
  glimpse()


# select relevant columns
dat <- dat %>%
  dplyr::select(Band, ch, Age) %>% 
  glimpse()


# save the cleaned data
saveRDS(dat, "data/tidy/capture_histories.rds")


## Environmental data ---------------------------------------------------------

env <- read.csv("data/raw/ENSO values.csv", skip = 1) %>% glimpse() # skip the first row (data source info)

# (any transformations needed)

# save the tidy environmental dat
saveRDS(env, "data/tidy/climate_indices.rds")

### END ###