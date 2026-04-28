
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


## Bayesian test graveyard


ch_matrix <- do.call(rbind, strsplit(sbsh$ch, split = "")) |> 
  apply(2, as.numeric)

enso <- sbsh_ddl$Phi$enso_cy_pc[!duplicated(sbsh_ddl$Phi$Time)]


# FIRST TEST
jags_data <- list(
  N = nrow(ch_matrix),          # number of individuals
  T = ncol(ch_matrix),          # number of time occasions
  y = ch_matrix,                # capture history
  first = first_capture,        # first capture
  enso = enso                   # covariate for survival (length T-1)
)

model {
  # Priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  alpha ~ dnorm(0, 0.001)  # constant recapture prob (on logit scale)
  
  # Individual loop
  for (i in 1:N) { # for every individual..
    # Latent alive/dead state
    for (t in first[i]:(T-1)) { # from their first capture onwards
      logit(phi[i, t]) <- beta0 + beta1 * enso[t] # fit a survival model to enso values
      z[i, t + 1] ~ dbern(z[i, t] * phi[i, t])
    }
    
    # Observation model
    for (t in (first[i] + 1):T) { # from their first capture onwards
      logit(p[i, t]) <- alpha # recapture is constant
      y[i, t] ~ dbern(z[i, t] * p[i, t])
    }
    
    # Initial state is alive at first capture
    z[i, first[i]] <- 1
  }
}

inits <- function() {
  z_init <- matrix(NA, nrow = jags_data$N, ncol = jags_data$T)
  
  for (i in 1:jags_data$N) {
    first <- jags_data$first[i]
    
    # Don't assign z[i,first], because it's deterministic in model
    for (t in (first + 1):jags_data$T) {
      if (jags_data$y[i, t] == 1) {
        z_init[i, t] <- 1  # Observed = alive
      } else {
        z_init[i, t] <- 1  # or 0, or a guess — but only if t > first
      }
    }
    
    # t < first[i] is not defined in model, so leave as NA
  }
  
  list(
    beta0 = 0,
    beta1 = 0,
    alpha = 0,
    z = z_init
  )
}
params <- c("beta0", "beta1", "alpha")


# Define first capture for each individual
get_first <- function(ch) min(which(ch == 1))
first_capture <- apply(ch_matrix, 1, get_first)
valid <- first_capture < jags_data$T
jags_data$first <- first_capture[valid]
jags_data$y <- jags_data$y[valid, ]
jags_data$N <- sum(valid)

fit <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "cjs_test.txt",
            n.chains = 3,
            n.iter = 5000,
            n.burnin = 1000,
            n.thin = 2)

plot(fit)