library(tidyverse)
library(janitor)
library(tictoc)
library(metadat)
library(metafor)
library(boot)

fit_selmodel <- function(dat, index = 1:nrow(dat), 
                         type = "stepfun", steps = 0.025) {
  
  # take subset of data
  boot_dat_cluster <- dat[index, ]
  
  # expand to one row per effect size
  boot_dat <- tidyr::unnest(boot_dat_cluster, data)
  
  # build run_selmodel
  run_sel_model <- function(dat, type, steps) {
    
    # initial random effects model
    RE_mod <- metafor::rma.uni(
      yi = yi, vi = vi, data = dat, 
      method = "ML"
    )
    
    # fit selection model
    res <- metafor::selmodel(
      RE_mod, 
      type = type, 
      steps = steps,
      skiphes = TRUE, # turn off SE calculation
      skiphet = TRUE # turn off heterogeneity test
    )
    
    # compile parameter estimates into a vector
    c(beta = res$beta[,1], 
      tau = sqrt(res$tau2),
      delta = if (type == "stepfun") res$delta[-1] else res$delta)
    
  }
  
  p <- 2L + length(steps)
  run_sel_model <- purrr::possibly(run_sel_model, otherwise = rep(NA_real_, p))
  
  # fit selection model, return vector
  run_sel_model(boot_dat, type = type, steps = steps)
  
}


# Example dataset

lehmann_dat <- 
  dat.lehmann2018 %>%
  clean_names() %>%
  mutate(study = str_split_fixed(short_title, pattern = "-", n = 2)[, 1]) %>%
  select(study, everything())


# Random effects model

RE_mod <- 
  rma.uni(yi, vi = vi, data = lehmann_dat, method = "ML") |>
  robust(cluster = study, clubSandwich = TRUE)
RE_mod

# Nest the data for each study

lehmann_nested <- nest_by(lehmann_dat, study, .key = "data")

fit_selmodel(lehmann_nested)



# Generate bootstraps

set.seed(20230222)

tic()

boots <- boot(
  data = lehmann_nested,
  statistic = fit_selmodel, steps = .025,
  R = 1999,
  parallel = "snow", ncpus = 8 # your mileage may vary
)

toc()

# For overall average ES
boot.ci(boots, type = "perc", index = 1)

# For heterogeneity
boot.ci(boots, type = "perc", index = 2)

# For selection weight
boot.ci(boots, type = "perc", index = 3)
