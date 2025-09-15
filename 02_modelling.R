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
clim  <- merge(clim, temp, by = "time") %>% filter(time != 2009)

glimpse(clim)

# make two time groups for recapture
recap_group <- data.frame(time = 2009:2024) %>% 
  mutate(recap_group = ifelse(time > 2019, "A", "B"))
  

## Check for collinearity -----------------------------------------------------

# following method by Grosbois et al. 2008, check the collinearity of variables before testing them.
m <- cor(clim)
test <-  cor.mtest(clim, conf.level = 0.95)
corrplot(m, p.mat = test$p, sig.level = 0.10, order = 'hclust')
# here several variables are collinear (any circle not crossed is a pair of collinear variables), and need to be combined.


# PDO values are all collinear with each other, we'll group them by running a Principal Component (PC) analysis with them and using this combined axis as our PDO variable.
pdo <- clim[, c("PDO_sep_oct_nov", "PDO_may_jun_jul", "PDO_prev_yr_sep_oct_nov", "PDO_prev_yr_may_jun_jul")]
pdo <- scale(pdo)
pca_result <- prcomp(pdo, center = TRUE, scale. = TRUE)
summary(pca_result) # the PC1 axis of the analysis describes the variation in the PDO values very well (84% of variation)
plot(pca_result, type = "l")
pdo_pc <- pca_result$x[,1] # we will now use PC1 as our new PDO variable


# AAO values in may-july (previous to current year) are collinear, we'll group them by running a Principal Component (PC) analysis as with PDO.
aao_mjj <- clim[, c("AAO_may_jun_jul", "AAO_prev_yr_may_jun_jul")]
aao_mjj <- scale(aao_mjj)
pca_result <- prcomp(aao_mjj, center = TRUE, scale. = TRUE)
summary(pca_result) # the PC1 axis of the analysis describes the variation in the AAO values well (69% of variation)
plot(pca_result, type = "l")
aao_mjj_pc <- pca_result$x[,1] # we will now use PC1 as our new PDO variable


# ENSO values are also collinear: may-july and sep-nov of the same year need to be grouped
# current year
enso_cy <- clim[, c("ENSO_may_jun_jul", "ENSO_sep_oct_nov")]
enso_cy <- scale(enso_cy)
pca_result <- prcomp(enso_cy, center = TRUE, scale. = TRUE)
summary(pca_result) # the PC1 axis of the analysis describes the variation in the PDO values very well (84% of variation)
plot(pca_result, type = "l") # again, the PC1 axis of the analysis describes the variation in the ENSO values very well (92% of variation)
enso_cy_pc <- pca_result$x[,1] # we will now use PC1 as our new PDO variable
# previous year
enso_py <- clim[, c("ENSO_prev_yr_may_jun_jul", "ENSO_prev_yr_sep_oct_nov")]
enso_py <- scale(enso_py)
pca_result <- prcomp(enso_py, center = TRUE, scale. = TRUE)
summary(pca_result) # the PC1 axis of the analysis describes the variation in the PDO values very well (84% of variation)
plot(pca_result, type = "l") # again, the PC1 axis of the analysis describes the variation in the ENSO values very well (93% of variation)
enso_py_pc <- pca_result$x[,1] # we will now use PC1 as our new PDO variable

# rebuild the climate variable dataset with clean, combined variables
clean_clim <- data.frame(
  time = clim$time,
  
  # the combined variables
  enso_cy_pc = enso_cy_pc,
  enso_py_pc = enso_py_pc,
  pdo_pc = pdo_pc,
  aao_mjj_pc = aao_mjj_pc,
  
  # the unchanged variables
  AAO_sep_oct_nov = clim$AAO_sep_oct_nov,
  AAO_prev_yr_sep_oct_nov = clim$AAO_prev_yr_sep_oct_nov,
  temp_may_jun_jul = clim$temp_may_jun_jul,
  temp_sep_oct_nov = clim$temp_sep_oct_nov,
  temp_prev_may_jun_jul = clim$temp_prev_may_jun_jul,
  temp_prev_sep_oct_nov = clim$temp_prev_sep_oct_nov
)

m <- cor(clean_clim)
test <-  cor.mtest(clean_clim, conf.level = 0.95)
corrplot(m, p.mat = test$p, sig.level = 0.10, order = 'hclust')
# here none of the candidate variables are collinear (except two AAO variables but ignoring that for now)

# plot check
clean_clim_long <- clean_clim %>%
  pivot_longer(
    cols = -time,  # everything except 'time'
    names_to = "variable",
    values_to = "value"
  )
ggplot(clean_clim_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# we can now move on


## Process the data -----------------------------------------------------------

# check CJS assumptions by running a goodness-of-fit test
ch <- sbsh %>% dplyr::select(`2010`:`2024`) %>% 
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()
n <- rep(1, nrow(sbsh))

overall_CJS(ch, n) # the GOF test for the whole model is not significant, we can go ahead with the CJS model (Gimenez et al. 2017)


## Test reference models ------------------------------------------------------

# again, following Grosbois et al. 2008, we are first assessing the general shape of survival. is it constant, it is logarithmic etc.
# with the best shape, we can then add the environmental variables to test.


# process the data
sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2010)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clean_clim) # climate indices are assumed to influence survival Phi, not recapture probability

# add recapture group to the design data
sbsh_ddl$p <- merge_design.covariates(sbsh_ddl$p, recap_group)


# we will test 4 survival shapes, with 2 recapture shapes
# constant survival
Phi.dot <- list(formula = ~1)

# fully time-dependent survival
Phi.time <- list(formula = ~factor(time))

# linear time survival
Phi.linear.time <- list(formula = ~time)

# quadratic time survival
Phi.quad.time <- list(formula = ~time + I(time^2))


# constant recapture
p.dot <- list(formula = ~1)

# time-dependent recapture
p.time <- list(formula = ~factor(time))


# now test the 8 reference models
sbsh_ddl$Phi$time <- as.numeric(as.character(sbsh_ddl$Phi$time))
sbsh_ddl$p$time   <- as.numeric(as.character(sbsh_ddl$p$time))
# Define formulas
phi.formulas <- list(
  const     = list(formula = ~1),
  time      = list(formula = ~factor(time)),
  linear    = list(formula = ~time),
  quadratic = list(formula = ~time + I(time^2))
)

p.formulas <- list(
  const = list(formula = ~1),
  time  = list(formula = ~factor(time))
)

# Cross all combinations as per Table 5
model.combos <- list(
  list(phi = phi.formulas$time,      p = p.formulas$const),  # 1
  list(phi = phi.formulas$linear,    p = p.formulas$const),  # 2
  list(phi = phi.formulas$const,     p = p.formulas$const),  # 3
  list(phi = phi.formulas$quadratic, p = p.formulas$const),  # 4
  list(phi = phi.formulas$time,      p = p.formulas$time),   # 5
  list(phi = phi.formulas$linear,    p = p.formulas$time),   # 6
  list(phi = phi.formulas$const,     p = p.formulas$time),   # 7
  list(phi = phi.formulas$quadratic, p = p.formulas$time)    # 8
)

# Initialize a list to hold models
models <- list()

# Loop through model combinations
for (i in seq_along(model.combos)) {
  models[[paste0("model", i)]] <- mark(
    data = sbsh_process,
    ddl = sbsh_ddl,
    model.parameters = list(
      Phi = model.combos[[i]]$phi,
      p   = model.combos[[i]]$p
    ),
    filename = paste0("data/rmark_outputs/reference_models/", "reference_model_", i),
    model.name = paste0("model", i),
    delete = FALSE   # keep the model in memory so collect.models() can find it
    
  )
}

# compare the reference models with AICc
reference_model_set <- model.table(models, use.AIC = TRUE)
reference_model_set # this ranks the reference models by AICc, telling us the Phi~1 and p ~factor(time) are the best to use


## Test covariate models ------------------------------------------------------

# standardise climate covariates
clim_std <- clean_clim
clim_std[ , -1] <- scale(clean_clim[ , -1])

clim_long <- clim_std %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

ggplot(clim_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  labs(x = "year",
       y = "standardised value",
       color = "covariate")

# process the data
sbsh_process <- process.data(sbsh, model = "CJS", begin.time = 2010)
sbsh_ddl <- make.design.data(sbsh_process)

# add climate indices values to the design data
sbsh_ddl$Phi <- merge_design.covariates(sbsh_ddl$Phi, clim_std) # climate indices are assumed to influence survival Phi, not recapture probability

# add recapture group to the design data
sbsh_ddl$p <- merge_design.covariates(sbsh_ddl$p, recap_group)



# create sets of linear and quadratic survivals to test
glimpse(sbsh_ddl)

# current covariates
covariates <- c("enso_cy_pc", "enso_py_pc", 
                "pdo_pc", 
                "aao_mjj_pc", "AAO_sep_oct_nov", "AAO_prev_yr_sep_oct_nov", 
                "temp_may_jun_jul", "temp_sep_oct_nov", "temp_prev_may_jun_jul", "temp_prev_sep_oct_nov")

# make quadratic covariates
for (cov in covariates) {
  sq_name <- paste0(cov, "2")
  sbsh_ddl$Phi[[sq_name]] <- sbsh_ddl$Phi[[cov]]^2
}

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
                          filename = paste0("data/rmark_outputs/covariate_models/", "linear_covariate_model_", cov),
                          model.name = paste0("model_", cov)
                          )
  
  mod_quadratic  <- mark(data = sbsh_process, ddl = sbsh_ddl,
                         model.parameters = list(Phi = list(formula = f_quadratic),
                                                 p = list(formula = ~ factor(time))),
                         silent = TRUE,
                         filename = paste0("data/rmark_outputs/covariate_models/", "quadratic_covariate_model_", cov),
                         model.name = paste0("model_", cov)
                         )
  
  # store results
  results[[paste0(cov, "_linear")]] <- mod_linear
  results[[paste0(cov, "_quadratic")]] <- mod_quadratic
}
test_model_set <- model.table(results, use.AIC = TRUE) # check which model comes out on top, AIC-wise
test_model_set # these are our tests models, AIC values find Phi~temp_prev_mjj to be


## Hypothesis testing ---------------------------------------------------------

# now we need to test whether the top model explains a significant amount of variation in survival
# Grosbois suggests 3 tests based on likelihood-ratio.

# set up what we need first: 
# deviance values
Dev_cst <- reference_model_set[8,]$Deviance   # Fcst
Dev_co  <- test_model_set[1,]$Deviance        # Fco
Dev_t   <- reference_model_set[4,]$Deviance   # Ft
# Number of parameters
J <- test_model_set[1,]$npar                  # number of parameters in covariate model
n <- reference_model_set[4,]$npar             # number of survival estimates in Ft (i.e number of parameters because it's the fully-time-dependent model)



# 1. we'll test the following null hypothesis: " the climate covariate model fits the data as well as the time-dependent model. "
lrt <- Dev_co - Dev_t # Fco - Ft in Grosbois.
lrt
df <- n - J
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value 
# we can't reject the null hypothesis. this means we can move to the second test


# 2. we'll now test a second null hypothesis: " the climate covariate has no effect on survival. "
lrt <- Dev_cst - Dev_co # Fcst - Fco in Grosbois
lrt
df <- J - 1
p_value <- pchisq(lrt, df = df, lower.tail = FALSE)
p_value
# we reject this null hypothesis, the p-value is highly significant, so the climate covariate has an effect on survival. this means we can do the third test.


# 3. this is an alternate null hypothesis to 2. " the climate covariate has no effect on survival. ", in case your 1. was significant.
c_hat <- (Dev_co - Dev_t) / (n - J)
F_stat <- ((Dev_cst - Dev_co) / (J - 1)) / c_hat
df1 <- J - 1
df2 <- n - J
p_value <- pf(F_stat, df1, df2, lower.tail = FALSE)

cat("F-statistic =", F_stat, "\n")
cat("Degrees of freedom =", df1, "and", df2, "\n")
cat("p-value =", p_value, "\n")
# the output is corrected for multiple testing... not sure if it's worth me adding it given that my 2. already answers that null hypothesis.



# select based on information theory, AIC values...

# reference models AIC
AIC_ft   <- reference_model_set[4, "AIC"]  # model5 (Phi ~ factor(time) + p ~ factor(time))
AIC_fcst <- reference_model_set[8, "AIC"]  # model3 (Phi ~ 1)

info_theory <- test_model_set %>% 
  dplyr::select(Phi, AIC) %>% 
  mutate( # delta AIC
    delta_AIC_co_cst = AIC - AIC_fcst,
    delta_AIC_co_t = AIC - AIC_ft
         ) %>% 
  # we'll test linear terms first
  mutate(
    support_co_cst = ifelse(delta_AIC_co_cst < -2, "linear covariate support", "no support"),
    support_co_t = ifelse(delta_AIC_co_t <= 2, "linear covariate support", "no support"),
    #support_quad = ifelse(grepl("\\+.*2", info_theory$Phi) & ___) # i don't really get what how to test the quadratic terms but skipping for now.
  ) %>% 
  glimpse()
# i feel like there's something dodgy here because all models have significant support....


## Bayesian approach

# haven't gotten there yet.








# for later :: parameter averages
sbsh.mod.avg <- model.average(results, vcv = TRUE)
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