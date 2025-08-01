
dat <- readRDS("data/tidy/LHI_FFSH_capture_histories_known_fate.rds") %>% glimpse()

sbsh_process <- process.data(dat, model = "Known", begin.time = 2009)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$S <- merge_design.covariates(sbsh_ddl$S, clim) # climate indices are assumed to influence survival Phi, not recapture probability

mark(sbsh_process, sbsh_ddl, model.parameters = list(S = list(formula = ~ Time)))

climate_vars <- c("ENSO_may_jun_jul", "ENSO_sep_oct_nov", 
                  "ENSO_prev_yr_may_jun_jul", "ENSO_prev_yr_sep_oct_nov", 
                  "PDO_sep_oct_nov", "PDO_may_jun_jul", 
                  "PDO_prev_yr_sep_oct_nov", "PDO_prev_yr_may_jun_jul", 
                  "AAO_sep_oct_nov", "AAO_may_jun_jul", 
                  "AAO_prev_yr_sep_oct_nov", "AAO_prev_yr_may_jun_jul", 
                  "temp_may_jun_jul", "temp_sep_oct_nov", 
                  "temp_prev_may_jun_jul", "temp_prev_sep_oct_nov")

# 2. Generate all 2-variable combinations
combos <- combn(climate_vars, 2, simplify = FALSE)

# 3. Loop over combinations and fit models
results <- list()

for (i in seq_along(combos)) {
  vars <- combos[[i]]
  form <- as.formula(paste("~", paste(vars, collapse = " + ")))
  
  model_name <- paste(vars, collapse = "_")
  
  cat("Fitting model:", model_name, "\n")
  
  results[[model_name]] <- mark(sbsh_process, sbsh_ddl,
                                model.parameters = list(S = list(formula = form)))
}
glimpse(results)
aic_df <- data.frame(
  model = names(results),
  AICc = sapply(results, function(x) x$results$AICc)
)

# Sort by lowest AICc
aic_df <- aic_df %>% arrange(AICc)
print(aic_df)
results$ENSO_prev_yr_sep_oct_nov_AAO_prev_yr_sep_oct_nov
names(results)
