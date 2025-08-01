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

sbsh  <- readRDS("data/tidy/LHI_FFSH_capture_histories.rds") %>% 
  rename(age_capture = Age) %>% 
  glimpse()
enso  <- readRDS("data/tidy/ENSO_tidy.rds") %>% glimpse()
pdo   <- readRDS("data/tidy/PDO_tidy.rds") %>% glimpse()
aao   <- readRDS("data/tidy/AAO_tidy.rds") %>% glimpse()
temp  <- readRDS("data/tidy/temp_tidy.rds") %>% glimpse()

clim  <- merge(enso, pdo, by = "time")
clim  <- merge(clim, aao, by = "time")
clim  <- merge(clim, temp, by = "time")

glimpse(clim)
  
## Process the data -----------------------------------------------------------

sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2010)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim) # climate indices are assumed to influence survival Phi, not recapture probability


# check CJS assumptions by running a goodness-of-fit test
library("R2ucare")

ch <- sbsh %>% dplyr::select(`2010`:`2024`) %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
n <- rep(1, nrow(sbsh))
overall_CJS(ch, n) # the GOF test for the whole model is not significant, we can go ahead with the CJS model (Gimenez et al. 2017)


## Define parameters to test --------------------------------------------------


sbsh.models <- function() {
  
  # Step 1: Define variable mapping for cleaner names
  var_map <- c(
    ENSO_sep_oct_nov = "ENSO_son",
    ENSO_may_jun_jul = "ENSO_mjj",
    ENSO_prev_yr_sep_oct_nov = "ENSO_son_prev",
    ENSO_prev_yr_may_jun_jul = "ENSO_mjj_prev",
    # PDO_sep_oct_nov = "PDO_son",
    # PDO_may_jun_jul = "PDO_mjj",
    # PDO_prev_yr_sep_oct_nov = "PDO_son_prev",
    # PDO_prev_yr_may_jun_jul = "PDO_mjj_prev",
    # AAO_sep_oct_nov = "AAO_son",
    # AAO_may_jun_jul = "AAO_mjj",
    # AAO_prev_yr_sep_oct_nov = "AAO_son_prev",
    # AAO_prev_yr_may_jun_jul = "AAO_mjj_prev",
    temp_sep_oct_nov = "temp_son",
    temp_may_jun_jul = "temp_mjj",
    temp_prev_sep_oct_nov = "temp_son_prev",
    temp_prev_may_jun_jul = "temp_mjj_prev",
    cohort = "cohort",
    age = "age",
    time = "time"
  )
  phi_vars <- names(var_map)
  
  # ------------------------------
  # Step 2: Generate Phi combinations and assign to function env
  # ------------------------------
  for (n in 1:2) {
    combos <- combn(phi_vars, n, simplify = FALSE)
    for (vars in combos) {
      formula_text <- paste("~", paste(vars, collapse = " + "))
      formula_obj <- list(formula = as.formula(formula_text))
      
      # Create name from var_map
      model_name <- paste0("Phi.", paste(var_map[vars], collapse = "_"))
      
      # Assign into the local environment (this function’s env)
      assign(model_name, formula_obj, envir = environment())
    }
  }
  
  # ------------------------------
  # Step 3: Define p models
  # ------------------------------
  p.dot <- list(formula = ~1, fixed = list(time = c(2010:2019), value = 0)) # here, the recapture rate 2010-2019 is set to zero, because no recapture effort.

  # ------------------------------
  # Step 4: Create model list and run
  # ------------------------------
  cml <- create.model.list("CJS")
  results <- mark.wrapper(cml, data = sbsh_process, ddl = sbsh_ddl)
  return(results)
}

summary(test)

# run that function
sbsh.results <- sbsh.models()
sbsh.results

phi.ensoMJJ.tempprMJJ <- list(formula = ~ENSO_may_jun_jul + temp_prev_may_jun_jul)
p.dot <- list(formula = ~1)
test <- mark(sbsh_process, sbsh_ddl, model.parameters = list(Phi = phi.ensoMJJ.tempprMJJ, p = p.dot))


# parameter averages
sbsh.mod.avg <- model.average(sbsh.results, vcv = TRUE)
sbsh.mod.avg.Phi <- model.average(sbsh.results, "Phi", vcv = TRUE)
sbsh.mod.avg.p <- model.average(sbsh.results, "p", vcv = TRUE)


sbsh.mod.avg$estimates
sbsh.mod.avg.Phi$estimates
sbsh.mod.avg.p$estimates

# survival estimate averaged over all models
plot(sbsh.mod.avg.Phi$estimates$time, sbsh.mod.avg.Phi$estimates$estimate,
     ylim = c(0, 1),
     pch = 19,
     xlab = "Time",
     ylab = "Estimate",
     main = "yikes",
     type = "b")
abline(h = 0.94)


summary_table <- sbsh.mod.avg.Phi$estimates %>%
  group_by(Time) %>%
  summarise(Average_Survival = mean(estimate, na.rm = TRUE)) %>%
  arrange(Time)

print(summary_table)

# extract best model
Phi.AAOp.Temp <- list(formula = ~AAO_prev_yr_sep_oct_nov + temp_sep_oct_nov) 
p.group <- list(formula = ~group) 

sbsh.best <- mark(sbsh_process, sbsh_ddl,
                  model.parameters = list(Phi = Phi.AAOp.Temp,
                                          p = p.group),
                  filename = "data/mark_outputs/best")
summary(sbsh.best)
### END ###