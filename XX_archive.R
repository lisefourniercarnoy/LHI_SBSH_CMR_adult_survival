###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2011-2024 CMR data

# script objective: keep old code that was once useful, but no more.

###############################################################################

rm(list = ls())

library(tidyverse)
library(stringr)
library(RMark)
library(R2ucare)
library(corrplot)
library(jagsUI)


## Principal Component (PC) in case covariates are collinear ------------------


# if covariate values are very collinear with each other, we'll group them by running a Principal Component (PC) analysis with them and using this combined axis as our variable.
cov <- clim[, c("cov1", "cov2")]
cov <- scale(cov)
pca_result <- prcomp(cov, center = TRUE, scale. = TRUE)
summary(pca_result) # check that the axis of the analysis describes the variation in the covariate values well
plot(pca_result, type = "l")
cov_pc <- pca_result$x[,1] # we will now use PC1 as our new combined variable

# rebuild the climate variable dataset with clean, combined variables
clean_clim <- data.frame(
  time = clim$time,
  
  # the combined variables
  cov = cov_pc,

  # the unchanged variables
  AAO_sep_oct_nov = clim$AAO_sep_oct_nov,
  AAO_jun_jul_aug = clim$AAO_jun_jul_aug,
  temp_sep_oct_nov = clim$temp_sep_oct_nov
)

m <- cor(clean_clim)
test <-  cor.mtest(clean_clim, conf.level = 0.95)
corrplot(m, p.mat = test$p, sig.level = 0.05, order = 'hclust')
# here fewer covariates should be collinear


## END ##
