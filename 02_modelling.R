###############################################################################

# authors: Jenn Lavers, Lise Fournier-Carnoy, Alex Bond
# project: Sable Shearwater, adult CMR survival
# data: LHI SBSH 2011-2024 CMR data

# script objective: run RMark analysis

###############################################################################
rm(list = ls())

library(tidyverse)
library(stringr)
library(RMark)

## Read in data ---------------------------------------------------------------

sbsh <- readRDS("data/tidy/capture_histories.rds")
clim <- readRDS("data/tidy/climate_indices.rds")

## Process the data -----------------------------------------------------------

sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2010, groups = "Age")

sbsh_ddl <- make.design.data(sbsh_process)

# add ENSO values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim) # ENSO is assumed to influence survival Phi, not recapture probability


## Define parameters to test --------------------------------------------------

# choose predictors to test
sbsh.models <- function()
{
  # select a few variables that could influence survival Phi
  Phi.dot           <- list(formula = ~1) # survival doesn't change
  Phi.ENSO          <- list(formula = ~ENSO_sep_oct_nov) # survival depends on ENSO
  Phi.cohort        <- list(formula = ~cohort) # survival depends on cohort
  Phi.age_cohort    <- list(formula = ~age + cohort) # survival depends on cohort and age
  Phi.time          <- list(formula = ~time) # survival depends on time
  
  # select what could influence recapture
  p.dot <- list(formula = ~1) # recapture probability doesn't change
  p.time <- list(formula = ~group) # recapture probability changes over time 
  
  cml <- create.model.list("CJS")
  results <- mark.wrapper(cml, 
                          data = sbsh_process, 
                          ddl = sbsh_ddl)
  return(results)
}

# run that function
sbsh.results <- sbsh.models()
sbsh.results

# parameter averages
sbsh.mod.avg <- model.average(sbsh.results, vcv = TRUE)
sbsh.mod.avg.Phi <- model.average(sbsh.results, "Phi", vcv = TRUE)

sbsh.mod.avg$estimates
sbsh.mod.avg.Phi$estimates

# Create a plot
plot(sbsh.mod.avg.Phi$estimates$time, sbsh.mod.avg.Phi$estimates$estimate,
     ylim = c(0, 1),
     pch = 19,
     xlab = "Time",
     ylab = "Estimate",
     main = "Survival Estimates over Time with CI",
     type = "b")


### END ###