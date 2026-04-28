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
library(R2ucare)
library(corrplot)
library(jagsUI)


## Read in data ---------------------------------------------------------------

sbsh  <- readRDS("data/tidy/LHI_FFSH_capture_histories.rds") %>% 
  rename(age_capture = Age) %>% 
  glimpse()

enso    <- readRDS("data/tidy/ENSO_tidy.rds") %>% glimpse()
pdo     <- readRDS("data/tidy/PDO_tidy.rds") %>% glimpse()
aao     <- readRDS("data/tidy/AAO_tidy.rds") %>% glimpse()
ao      <- readRDS("data/tidy/AO_tidy.rds") %>% glimpse()
temp    <- readRDS("data/tidy/temp_tidy.rds") %>% glimpse()
temp_w  <- readRDS("data/tidy/temp_w_tidy.rds") %>% dplyr::filter(time>=2013) %>% glimpse()

clim_list = list(enso, 
                 #pdo, 
                 #aao, 
                 ao, 
                 #temp, #this is temperature near Sydney which is irrelevant because birds arent there Jul-Sep
                 temp_w)
clim <- clim_list %>% 
  reduce(~ inner_join(.x, .y, by = "time"))

clim <- clim %>% dplyr::filter(time >= 2013)
glimpse(clim)

# make two time groups for recapture
recap_group <- data.frame(time = 2013:2024) %>% 
  mutate(recap_group = ifelse(time > 2019, "A", "B"))


## Check for collinearity -----------------------------------------------------

# following method by Grosbois et al. 2008, check the collinearity of variables before testing them.
m <- cor(clim)
test <-  cor.mtest(clim, conf.level = 0.95)
corrplot(m, p.mat = test$p, sig.level = 0.1, order = 'hclust')
# any circle not crossed is a pair of collinear variables.
# time and SST_w are collinear, but we will ignore this because we

# see XX_archive script to fix collinear covariates

# we can now move on !


## Process the data -----------------------------------------------------------

# we are analysing an open-population system, and the model we fit has assumptions that need to be checked.

# check CJS assumptions by running a goodness-of-fit test
ch <- sbsh %>% dplyr::select(`2013`:`2024`) %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
n <- rep(1, nrow(sbsh))

# I'm following https://jamesepaterson.github.io/jamespatersonblog/2020-05-20_gof_for_CJS

# test 1: 'is there evidence of animals having equal capture probabilities and equal survival'.
overall_CJS(ch, n) # test is not significant, first assumption checked. We can go ahead with the CJS model (Gimenez et al. 2017)
# if this was significant, our data would not be suited for a CJS model, and results would be very wrong.

# if test 1 is significant, run test 2-4 to understand whether recapture depends on when animal was marked or whether marking affects survival.
# test 2: 'does recapture deend on when an animal was first marked?' 
test2ct(ch, n) # here p-val is not significant. in $details, some tests cannot be conducted because of low sample sizes.
test2cl(ch, n) # again, p-val is not significant. in $details, some tests cannot be conducted because of low sample sizes.

# test 3: 'does marking affect survival?'
test3sr(ch, n) # p-val not significant, same thing for $details

#test 4: 'for animals seen again, does when they are recaptured depend on whether they were marked on or before time t?'
test3sm(ch, n)

# We're all good to go. Gimenez et al. 2017 confirms this.


## Checking Variance Inflation Factor (c-hat) ---------------------------------

# c-hat can help adjust model output based on extra-binomial noise (more variance in capture history data) or overdispersion (larger-than-normal discrepancies between observed and predicted values for the binomial model).
# here, an estimation of c-hat is Chi2 / df, and is valid for time-dependent models. (which i initially tested here)
# which is essentially the Chi2 and df from the overall_CJS call just above.

test <- overall_CJS(ch, n)
test$chi2/test$degree_of_freedom # simple calculation from https://jamesepaterson.github.io/jamespatersonblog/2020-06-30_chat_for_CJS.html

# another estimation from gentle intro to mark, chapter 5, section 5.5.1
t2_chi2 <- test2ct(ch, n)$test2ct[1]
t3_chi2 <- test3sr(ch, n)$test3sr[1]
t2_df <- test2ct(ch, n)$test2ct[2]
t3_df <- test3sr(ch, n)$test3sr[2]

(t2_chi2 + t3_chi2)/(t2_df + t3_df) # chi2 of (test2 + test3) divided by degrees of freedom of (test2 + test3).

# a value ~1 indicates good variance given the model, <1 is underdispersed, and >1 overdispersed
# here we have c-hat << 1, which is likely because of small sample size.
# underdispersion is rare, and there doesn't seem to be either 1. a biological reason for it or 2. a fix for it.
# for this reason, we will set c-hat to 1 and just move on.


## Test reference models ------------------------------------------------------

# again, following Grosbois et al. 2008, we are first assessing the general shape of survival. is it constant, it is logarithmic etc.
# in addition, we will add recapture group (pre-2019 and post-2019) in the test formulas. if the formula is not good, we can ignore recapture group for the rest of the analysis.
# with the best shape, we can then add the environmental variables to test.

# process the data
sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2013)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim) # climate indices are assumed to influence survival Phi, not recapture probability

# add recapture group to the design data
sbsh_ddl$p <- merge_design.covariates(sbsh_ddl$p, recap_group)


# we will test 4 survival shapes, with 4 recapture shapes
Phi.dot <- list(formula = ~1) # constant survival
Phi.linear.time <- list(formula = ~time) # linear time survival
Phi.time <- list(formula = ~factor(time)) # fully time-dependent survival
Phi.quad.time <- list(formula = ~time + I(time^2)) # quadratic time survival

p.dot <- list(formula = ~1) # constant recapture
p.time <- list(formula = ~factor(time)) # time-dependent recapture
p.linear.time <- list(formula = ~time) # linear time survival
p.rg <- list(formula = ~recap_group) # recapture-group-dependent recapture -- this is not in Grosbois et al.
p.t.rg <- list(formula = ~recap_group + time) # recapture dependent on recapture group and linear time -- this is not in Grosbois et al.


# now test the reference models
sbsh_ddl$Phi$time <- as.numeric(as.character(sbsh_ddl$Phi$time))
sbsh_ddl$p$time   <- as.numeric(as.character(sbsh_ddl$p$time))

# define formulas
phi.formulas <- list(
  const     = list(formula = ~1),
  time      = list(formula = ~factor(time)),
  linear    = list(formula = ~time),
  quadratic = list(formula = ~time + I(time^2))
)

p.formulas <- list(
  const       = list(formula = ~1),
  group       = list(formula = ~recap_group),
  const_time  = list(formula = ~time),
  group_time  = list(formula = ~recap_group + time),
  full_time   = list(formula = ~factor(time))
)

# cross all combinations as per Table 5 in Grosbois et al. 
model.combos <- list(
  list(phi = phi.formulas$time,      p = p.formulas$const),  # 1
  list(phi = phi.formulas$linear,    p = p.formulas$const),  # 2
  list(phi = phi.formulas$const,     p = p.formulas$const),  # 3
  list(phi = phi.formulas$quadratic, p = p.formulas$const),  # 4
  
  list(phi = phi.formulas$time,      p = p.formulas$group),  # 5
  list(phi = phi.formulas$linear,    p = p.formulas$group),  # 6
  list(phi = phi.formulas$const,     p = p.formulas$group),  # 7
  list(phi = phi.formulas$quadratic, p = p.formulas$group),  # 8
  
  list(phi = phi.formulas$time,      p = p.formulas$const_time),   # 9
  list(phi = phi.formulas$linear,    p = p.formulas$const_time),   # 10
  list(phi = phi.formulas$const,     p = p.formulas$const_time),   # 11
  list(phi = phi.formulas$quadratic, p = p.formulas$const_time),   # 12
  
  list(phi = phi.formulas$time,      p = p.formulas$group_time),   # 13
  list(phi = phi.formulas$linear,    p = p.formulas$group_time),   # 14
  list(phi = phi.formulas$const,     p = p.formulas$group_time),   # 15
  list(phi = phi.formulas$quadratic, p = p.formulas$group_time),   # 16
  
  list(phi = phi.formulas$time,      p = p.formulas$full_time),   # 17
  list(phi = phi.formulas$linear,    p = p.formulas$full_time),   # 18
  list(phi = phi.formulas$const,     p = p.formulas$full_time),   # 19
  list(phi = phi.formulas$quadratic, p = p.formulas$full_time)    # 20
)

# Initialize a list to hold models
models <- list()

# create a folder to store all outputs
dir_name <- paste0("data/rmark_outputs/reference_models/", today())
dir.create(path = dir_name)

# Loop through model combinations
for (i in seq_along(model.combos)) {
  models[[paste0("model", i)]] <- mark(
    data = sbsh_process,
    ddl = sbsh_ddl,
    model.parameters = list(
      Phi = model.combos[[i]]$phi,
      p   = model.combos[[i]]$p
    ),
    
    filename = paste0(dir_name, "/linear_covariate_model_", i),
    model.name = paste0("model", i),
    delete = FALSE, # keep the model in memory so collect.models() can find it
    chat = 1
  )
}


# compare the reference models with AICc
reference_model_set <- model.table(models, use.AIC = FALSE) # use.AIC = FALSE, to use corrected AIC (AICc) which is more appropriate for small sample sizes
reference_model_set # this ranks the reference models by AICc, telling us the Phi~1 and p~factor(time) are the best to use
saveRDS(reference_model_set, paste0(dir_name, "/reference_models_table.rds"))

reference_model_set <- readRDS(paste0(dir_name, "/reference_models_table.rds"))


## Test covariate models ------------------------------------------------------

# we will use the top reference model to test covariates: Phi~covariate and p~factor(time)

# standardise climate covariates
clim_std <- clim
clim_std[ , -1] <- scale(clim[ , -1])

clim_long <- clim_std %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

ggplot(clim_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "year",
       y = "standardised value",
       color = "covariate") +
  theme_minimal()

# process the data
sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2013)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim_std) # climate indices are assumed to influence survival Phi, not recapture probability

# add recapture group to the design data
sbsh_ddl$p <- merge_design.covariates(sbsh_ddl$p, recap_group)


# create sets of linear and quadratic survivals to test
glimpse(sbsh_ddl)
names(clim)

# current covariates
covariates <- names(clim[,-1])

# make quadratic covariates
for (cov in covariates) {
  sq_name <- paste0(cov, "2")
  sbsh_ddl$Phi[[sq_name]] <- sbsh_ddl$Phi[[cov]]^2
  
}

# create a folder to store all outputs
dir_name <- paste0("data/rmark_outputs/covariate_models/", today())
dir.create(path = dir_name)

# make a loop to test the linear and quadratic models of each of the variables
results <- list()
for (cov in covariates) {
  # define formulas
  f_linear        <- as.formula(paste0("~ ", cov))
  f_quadratic     <- as.formula(paste0("~ ", cov, " + ", cov, "2"))
  
  # fit models (using the p formula we found in the reference model)
  mod_linear      <- mark(data = sbsh_process, ddl = sbsh_ddl,
                          model.parameters = list(Phi = list(formula = f_linear),
                                                  p = list(formula = ~ factor(time))),
                          silent = TRUE,
                          filename = paste0(dir_name, "/linear_covariate_model_", cov),
                          model.name = paste0("model_", cov),
                          chat = 1
                          )
  
  mod_quadratic  <- mark(data = sbsh_process, ddl = sbsh_ddl,
                         model.parameters = list(Phi = list(formula = f_quadratic),
                                                 p = list(formula = ~ factor(time))),
                         silent = TRUE,
                         filename = paste0(dir_name, "/quadratic_covariate_model_", cov),
                         model.name = paste0("model_", cov),
                         chat = 1
                         )
  
  # store results
  results[[paste0(cov, "_linear")]] <- mod_linear
  results[[paste0(cov, "_quadratic")]] <- mod_quadratic
}
test_model_set <- model.table(results, use.AIC = FALSE) # as above, use.AIC = FALSE to use AICc (better for small sample sizes)
test_model_set # these are our tests models
saveRDS(test_model_set, "data/rmark_outputs/reference_models/test_models_table.rds")

test_model_set <- readRDS("data/rmark_outputs/reference_models/test_models_table.rds")


## Hypothesis testing ---------------------------------------------------------

# now we need to test whether the top model explains a significant amount of variation in survival
# Grosbois suggests 3 tests based on likelihood-ratio.

# set up what we need first: 
# deviance values
Dev_cst <- reference_model_set[reference_model_set$Phi == "~1" & reference_model_set$p == "~1",]$Deviance   # Fcst
Dev_co  <- test_model_set[1,]$Deviance        # Fco (row 1 = lowest AICc) - or QAICc/QDeviance if c-hat has been adjusted
Dev_t   <- reference_model_set[reference_model_set$Phi == "~factor(time)" & reference_model_set$p == "~factor(time)",]$Deviance # Ft

# Number of parameters
J <- test_model_set[1,]$npar # number of parameters in covariate model
n <- reference_model_set[reference_model_set$Phi == "~factor(time)" & reference_model_set$p == "~factor(time)",]$npar # number of survival estimates in Ft (i.e number of parameters because it's the fully-time-dependent model)



# 1. we'll test the following null hypothesis: " the climate covariate model fits the data as well as the time-dependent model. "
lrt <- Dev_co - Dev_t # Fco - Ft in Grosbois.
lrt
df <- n - J
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value 
# p>0.05, we can't reject the null hypothesis. this means we can move to the second test


# 2. we'll now test a second null hypothesis: " the climate covariate has no effect on survival. "
lrt <- Dev_cst - Dev_co # Fcst - Fco in Grosbois
lrt
df <- J - 1
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value
# p<0.05, we reject this null hypothesis, the p-value is highly significant, so the climate covariate has an effect on survival. this means we can do the third test.


# 3. this is an alternate null hypothesis to 2. " the climate covariate has no effect on survival. ", in case your 1. was significant.
# c_hat <- (Dev_co - Dev_t) / (n - J)
# F_stat <- ((Dev_cst - Dev_co) / (J - 1)) / c_hat
# df1 <- J - 1
# df2 <- n - J
# p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)
# 
# cat("F-statistic =", F_stat, "\n")
# cat("Degrees of freedom =", df1, "and", df2, "\n")
# cat("p-value =", p_value, "\n")
# # the output is corrected for multiple testing... if p>0.05, so we reject this null hypothesis, the climate covariate has an effect on survival.



# select based on information theory, AIC values...

# reference models AIC
AICc_ft   <- reference_model_set[reference_model_set$Phi == "~factor(time)" & reference_model_set$p == "~factor(time)", "AICc"]
AICc_fcst <- reference_model_set[reference_model_set$Phi == "~1" & reference_model_set$p == "~1", "AICc"]

# below compares every model's AIC with the constant model (fcst) and the time-dependent model (ft)
# if the delta_AICc > 2, it means the model is substantially better than the fcst/ft models
info_theory <- test_model_set %>% 
  dplyr::select(Phi, AICc) %>% 
  mutate( # delta AIC
    delta_AICc_co_cst = AICc - AICc_fcst,
    delta_AICc_co_t = AICc - AICc_ft
         ) %>% 
  # we'll test linear terms first
  mutate(
    support_co_cst = ifelse(delta_AICc_co_cst < -2, "linear covariate support", "no support"),
    support_co_t = ifelse(delta_AICc_co_t <= 2, "linear covariate support", "no support"),
  ) %>% 
  glimpse()


## obtain survival estimates --------------------------------------------------

test_model_set[[1]][1] # best model
best_model <- results$AO_jun_jul_aug_linear # extract the best model

# extract survival estimates
phi_estimates <- get.real(best_model, "Phi", vcv = TRUE) # extract phi
phi_estimates$estimates  # includes estimate, se, lcl, ucl

# now we need to obtain survival estimates
phi_unique <- phi_estimates$estimates %>%
  distinct(par.index, .keep_all = TRUE) %>%
  arrange(par.index) %>%
  mutate(year = 2013:(2013 + n() - 1))  # add readable year labels

phi_unique

ggplot(phi_unique, aes(x = year, y = estimate)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.5) +
  scale_x_continuous(breaks = c(2013:2023)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year", y = "Annual survival probability") +
  theme_light()

mean(phi_unique$estimate) # 0.8368935. yoikes


# extract recapture estimates
p_estimates <- get.real(best_model, "p", vcv = TRUE) # extract phi
p_estimates$estimates  # includes estimate, se, lcl, ucl

# now we need to obtain survival estimates
p_unique <- p_estimates$estimates %>%
  distinct(par.index, .keep_all = TRUE) %>%
  arrange(par.index) %>%
  mutate(year = 2013:(2013 + n() - 1))  # add readable year labels

p_unique

ggplot(p_unique, aes(x = year, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.3) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year", y = "Annual recapture probability") +
  theme_minimal()

mean(p_unique$estimate) # 0.09477007



## Bayesian approach ----------------------------------------------------------

# okay so. i use the crashcourse on Bayes' theorem to understand this: https://www.youtube.com/watch?v=9TDjifpGj-k&t=67s
# Bayesian statistics is: updating your beliefs about something when new info comes along. PRIOR beliefs update to your POSTERIOR belief
# we update what we think survival is, based on the new climate information we have.

# to test whether a climate covariate is actually useful, the Markov Chain Monte Carlo is going to test many models with different parameters of covariates.
# for example it'll test A. survival ~ 0.1*temp, and B. survival ~ 0.2*temp, and will tell you how likely A and B are, given your data.
# the MCMC 'samples from the POSTERIOR distribution' of parameter values. If the parameter value has a higher likelihood with our data, those values get higher POSTERIOR probability.
# so we test whether with a small effect of temp (in model A) is likely given the recapture data we have. same with a slighly high effect of temp (like in B) etc.
# this gives us a curve: x being the temp parameter, y being the likelihood.
# if the parameter value of 0.1 is given high posterior probability (like in A.), it means that the data highly supports this value.
# say it tests values -1 to 1, it'll give us a 95% confidence interval (we're 95% sure the real parameter value is between _ and _).
# if the confidence interval contains 0, then we consider the covariate doesn't have statistical support, because there's a good chance there is no effect of the covariate on survival.

glimpse(sbsh)
glimpse(sbsh_ddl)

# we're going to test all the models that were within 2 AIC points of each other, to see the amount of support for each
top_covariates <- test_model_set$Phi[test_model_set$DeltaAICc < 2]
top_covariates <- gsub("~", "", top_covariates)

top_covariates

covariate_list <- list(
  
  AO_jun_jul_aug = as.matrix(sbsh_ddl$Phi[,"AO_jun_jul_aug"])
  # temp_quad   = cbind(
  #   sbsh_ddl$Phi[,"temp_jun_oct"],
  #   sbsh_ddl$Phi[,"temp_jun_oct"]^2
  # )
)

model_names <- names(covariate_list); model_names

# obtain the first year of capture for each individual
get_first <- function(ch) min(which(ch == 1))
first_capture <- apply(ch, 1, get_first)
ch_matrix <- as.matrix(ch)

# okay here we go.

# below is the function that we're using. it uses JAGS ('just another gibbs sampler') to fit the MCMC.
# important thing to understand here is: 
# the model will use the matrix y (the observed recaptures) to estimate whether the individual is alive or not when we don't recapture it (in matrix z, the 'latent state')
# for example, in capture history '1010', it will estimate how likely it is that the individual survived years 2 and 4, given it survived year 3. 


writeLines("model {
  
  beta0 ~ dnorm(0, 1)  # Prior for intercept on survival
  alpha ~ dnorm(0, 1)  # contant recapture probability
  # for (t in 1:T) {        # Time-specific recapture probabilities
  #   alpha[t] ~ dnorm(0, 1)
  # }
  
  for (k in 1:K) {        # Priors for covariate effects on survival
    beta[k] ~ dnorm(0, 1)
  }
  
  for (i in 1:N) {
    z[i, first[i]] <- 1  # Known to be alive at first capture
    
    for (t in (first[i] + 1):T) {
      logit(phi[i, t - 1]) <- beta0 + inprod(beta[1:K], X[t - 1, 1:K])
      z[i, t] ~ dbern(z[i, t - 1] * phi[i, t - 1])
    }
    
    for (t in first[i]:T) {
      # logit(p[i, t]) <- alpha[t] # time dependent recapture probability
      logit(p[i, t]) <- alpha # constant recapture probability

      y[i, t] ~ dbern(z[i, t] * p[i, t])
    }
  }
}", "data/JAGS_models/cjs_model.txt")


results <- list() # initiate an object to store all results.

set.seed(123)

# THIS LOOP TAKES A WHILE TO RUN
for (i in seq_along(covariate_list)) { # for each covariate model... 

  X <- covariate_list[[i]] # ...select the covariate in the list... 
  if (i %in% c(2, 4)) {X[,2] <- X[,2] - mean(X[,2])} # centre the temp quad
  
  if (is.vector(X)) X <- matrix(X, ncol = 1) # convert to a matrix format, if it's a vector it doesn't work 
  
  jags_data_i <- list( # this is essentially all the data that the model is going to use 
    N     = nrow(ch_matrix), # number of individuals 
    T     = ncol(ch_matrix), # number of years in the data 
    y     = ch_matrix, # the observed recapture data 
    first = first_capture, # the year of first capture 
    K     = ncol(X), # number of covariates tested 
    X     = X # the covariate values of the current model 
    ) # Update initial values function for beta vector of length K
  
  inits <- function() { # here we set initial values. these are assumption about the survival of the individuals, and will serve in the
    z_init <- matrix(NA, nrow = jags_data_i$N, ncol = jags_data_i$T) # initial values matrix for the latent state.
    for (j in 1:jags_data_i$N) { # for every individual
      first_j <- jags_data_i$first[j] # obtain the first capture,
      if (first_j < jags_data_i$T) { # then for each year after that,
        for (t in (first_j + 1):jags_data_i$T) { # we assume that the bird is alive
          z_init[j, t] <- 1
        }
      }
    }
    
    list(beta0 = 0, 
         beta = rep(0, jags_data_i$K), 
         z = z_init
         ) 
    }
  
  # Parameters to save 
  params <- c("beta0", "beta", "alpha", "deviance") 
  
  # Run JAGS 
  fit <- jags(data = jags_data_i, 
              inits = inits, 
              parameters.to.save = params, 
              model.file = "cjs_model.txt", 
              n.chains = 3, 
              n.iter = 100000, 
              n.burnin = 20000, 
              n.thin = 1) 
  results[[model_names[i]]] <- fit 
  } # this loop will test the likelihood of the observed recapture data, given different values of beta0 (intercept for survival), beta (effect of covariate) and alpha (recapture probability).

names(results)

# each element in the 'results' list is a CJS model, the performance of which we're assessing below.
plot(results$temp_jun_oct) # all good
plot(results$temp_quad) # all good

plot(results$PDO_jul_aug_sep) # all good
plot(results$pdo_quad) # all good

plot(results$AAO_jul_aug_sep) # all good
plot(results$AAO_quad) # all good

plot(results$ENSO_jul_aug_sep) # all good
plot(results$enso_quad) # all good


# we are looking for each chain's values (red, green and blue lines) to overlap well. here, they're all looking good (yay!)
# the R-hat value should also be <1.1

## obtain the converged data --------------------------------------------------

dat <- results[c(1:7)] # remove the enso_quad which doesnt want to converge.

# see which covariates are highly likely to affect survival
extract_beta_ci <- function(model_name, model_result) {
  summary_df <- as.data.frame(model_result$summary)
  beta_rows <- rownames(summary_df)[grepl("^beta(\\[1\\])?$", rownames(summary_df))]
  
  summary_df[beta_rows, ] %>%
    rownames_to_column("parameter") %>%
    mutate(model = model_name)
}

# Combine summaries from all models
beta_summaries <- purrr::map2_dfr(names(dat), dat, extract_beta_ci)

beta_summaries <- beta_summaries %>%
  mutate(
    model = factor(model, levels = names(results)),
    sig = ifelse(`2.5%` > 0 | `97.5%` < 0, "yes", "no"),
    parameter = ifelse(parameter == "beta", "beta[1]", parameter)  # for consistency
  )

ggplot(beta_summaries, aes(x = mean, y = model, color = sig)) +
  geom_point() +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~parameter, scales = "free_x") +
  scale_color_manual(values = c("yes" = "blue", "no" = "red")) +
  labs(
    title = "Posterior estimates of beta (with 95% Credible Intervals)",
    x = "Estimate (mean)", y = "Model",
    color = "Excludes 0?"
  ) +
  theme_minimal()


# see the survival estimate for all models 

extract_model_parameters <- function(model_result) {
  summary_df <- model_result$summary %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter") %>%
    filter(grepl("^beta0$|^beta$|^beta\\[\\d+\\]$|^alpha$|^alpha\\[\\d+\\]$", parameter)) %>%
    select(parameter, mean, sd, `2.5%`, `97.5%`)
  
  return(summary_df)
}

model_exclude <- c("enso_quad")
dat <- dat[!names(dat) %in% model_exclude]
subset_results <- lapply(dat, extract_model_parameters)


parameter_summary <- bind_rows(subset_results, .id = "model")

print(parameter_summary)

beta0_all <- parameter_summary %>% 
  dplyr::filter(parameter == "beta0") %>% 
  mutate(beta0_t = 1 / (1 + exp(-mean)) )


# fit models again to get survival
get_phi <- function(model) {
  phi_df <- model$results$real
  phi_df <- phi_df[grepl("^Phi", rownames(phi_df)), ]
  return(phi_df)
}

match <- data.frame(
  model = sort(names(results)),
  covariate = sort(top_covariates[!top_covariates %in% c("ENSO_jul_aug_sep + ENSO_jul_aug_sep2")])
)
best_mod <- list()
for (COV in names(results)) {
  i <- match(COV, names(results))
  cov <- match[i, "covariate"]
  
  best_mod[[COV]] <- mark(data = sbsh_process, ddl  = sbsh_ddl,
                          model.parameters = list(Phi = list(formula = as.formula(paste("~", cov))),
                                                  p   = list(formula = ~ factor(time))),
                          silent = TRUE,
                          filename   = paste0("data/rmark_outputs/best_models/", COV),
                          model.name = paste0("model_", COV)
  )
}





# below is archive

# Get only the Phi rows from each model's real parameter estimates
get_phi <- function(model) {
  phi_df <- model$results$real
  phi_df <- phi_df[grepl("^Phi", rownames(phi_df)), ]
  return(phi_df)
}
summary(best_mod1)

real_df <- best_mod3$results$real
phi_df <- real_df[grepl("^Phi", rownames(real_df)), ]
phi_df$year <- as.numeric(sub(".*t", "", rownames(phi_df)))
phi_plot <- phi_df[!duplicated(phi_df$year), ]
ggplot(phi_plot, aes(x = year, y = estimate)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl),
              fill = "steelblue", alpha = 0.3) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Year",
    y = expression(" Annual survival probability ("*Phi*")"),
    title = "SST (SON, quad), Estimated annual survival with 95% confidence intervals"
  ) +
  theme_classic()


phi1 <- get_phi(best_mod1)
phi2 <- get_phi(best_mod2)
phi3 <- get_phi(best_mod3)
phi4 <- get_phi(best_mod4)
phi5 <- get_phi(best_mod5)
phi6 <- get_phi(best_mod6)


n_years <- nrow(phi1)
years <- 2013:2023

phi1$year <- years
phi2$year <- years
phi3$year <- years
phi4$year <- years
phi5$year <- years
phi6$year <- years


surv_df <- data.frame(
  year = years,
  model1 = phi1$estimate,
  model2 = phi2$estimate,
  model3 = phi3$estimate,
  model4 = phi4$estimate,
  model5 = phi5$estimate,
  model6 = phi6$estimate
)

# Calculate average survival per year across models
surv_df <- surv_df %>%
  mutate(mean_survival = rowMeans(select(., model1, model2)))


ggplot(surv_df, aes(x = year)) +
  geom_point(aes(y = model1), color = "lightblue") +
  geom_line(aes(y = model1), color = "lightblue") +
  
  geom_point(aes(y = model2), color = "lightblue2") +
  geom_line(aes(y = model2), color = "lightblue2") +
  
  
  geom_point(aes(y = model3), color = "lightblue3") +
  geom_line(aes(y = model3), color = "lightblue3") +

  geom_line(aes(y = mean_survival, color = "Average"), color = "red", size = 1.2) +
  ylim(0, 1) +
  labs(
    title = paste0("average survival = ", 
                  round(mean(mean(surv_df$model1), mean(surv_df$model2)
    ), 2)
    ),
    x = "Year",
    y = "Survival (Φ)"
  ) +
  theme_minimal()



# for later :: parameter averages
sbsh.mod.avg <- model.average(best_mod, vcv = TRUE)
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