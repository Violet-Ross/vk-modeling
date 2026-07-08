library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = 8)

# Read in and prepare data
synth <- fitting_data
id <- synth[[1]]
x <- synth[[2]]
y <- synth[[3]]

J <- max(id)                  # individuals
N <- nrow(synth)

data_list <- list(
  N = N,
  J = J,
  id = id,
  x = x,
  y = y,
  sd = prior_sd
)

# Fit model
output <- capture.output(fit <- stan(
  file = "src/viral_bhm_vary_priors.stan",
  data = data_list,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  seed = 123, 
  refresh = 0
))

# save fit model
#saveRDS(fit, file = "throughput/viral_bhm.rds")