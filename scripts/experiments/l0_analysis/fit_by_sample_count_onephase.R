# fit_by_sample_count.R
# Fix l_0 = 4, prior_sd = 9. Simulate one dataset, split trajectories by the
# number of observed positive samples, then fit a single-stage BHM to each
# non-empty group.  Full Stan fit objects are saved as .rds files.

library(rstan)
library(readr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# ── Constants ─────────────────────────────────────────────────────────────────

LOD      <- 3
L_0      <- 5
PRIOR_SD <- 9
STAN_FILE <- "src/viral_bhm_vary_priors.stan"

# ── Functions (unchanged from run_condition.R) ────────────────────────────────

generate_data <- function(LOD, l_0) {
  num_people    <- 500
  num_days      <- 500
  prev          <- 0.3
  min_x         <- 0
  max_x         <- 50
  
  alpha         <- LOD
  mu_beta1      <- 2;   mu_beta2 <- -4;  mu_psi <- 5
  sigma_beta1   <- 0.1; sigma_beta2 <- 0.1; sigma_psi <- 0.1
  sigma_epsilon <- abs(rnorm(1, 0, sd = 0.5))
  
  pop_params <- data.frame(mu_beta1, mu_beta2, mu_psi,
                           sigma_beta1, sigma_beta2, sigma_psi, sigma_epsilon)
  
  infection_matrix <- matrix(0, nrow = num_people, ncol = num_days)
  is_infx          <- rbinom(num_people, size = 1, prev)
  x_vals           <- seq(from = min_x, to = max_x, length.out = max_x - min_x)
  
  for (row_ind in 1:num_people) {
    if (is_infx[row_ind] == 1) {
      beta1 <- rnorm(1, mu_beta1, sigma_beta1)
      beta2 <- rnorm(1, mu_beta2, sigma_beta2)
      psi   <- rnorm(1, mu_psi,   sigma_psi)
      
      y_vals <- sapply(x_vals, function(x) {
        obs_mu <- alpha + (beta1 * x) + (beta2 * (x - psi) * (x > psi))
        rnorm(1, obs_mu, sigma_epsilon)
      })
      
      pos_y_vals  <- y_vals[y_vals > 0]
      infx_length <- length(pos_y_vals)
      start_day   <- sample(1:(num_days - infx_length), 1)
      infection_matrix[row_ind, start_day:(start_day + infx_length - 1)] <- pos_y_vals
    }
  }
  
  mat        <- data.frame(infection_matrix)
  test_times <- seq(from = 1, to = ncol(mat), by = l_0)
  
  traj_values <- c(); index_values <- c(); time_values  <- c()
  person_values <- c(); date_values <- c()
  positive_rows <- 1
  
  for (row_ind in 1:nrow(mat)) {
    if (any(mat[row_ind, test_times] > 0)) {
      row          <- mat[row_ind, test_times]
      traj_indices <- which(row > 0)
      traj         <- unlist(row[traj_indices])
      
      traj_values   <- append(traj_values,   traj)
      index_values  <- append(index_values,  rep(positive_rows, length(traj)))
      time_values   <- append(time_values,   (seq(1, length(traj)) - 1) * l_0)
      person_values <- append(person_values, rep(row_ind, length(traj)))
      date_values   <- append(date_values,   test_times[traj_indices])
      
      positive_rows <- positive_rows + 1
    }
  }
  
  all_trajectories <- data.frame(index_values, time_values,
                                 traj_values,  person_values, date_values)
  list(all_trajectories = all_trajectories, pop_params = pop_params)
}


fit_model <- function(all_trajectories, prior_sd, stan_file, seed = 123) {
  id <- all_trajectories$index_values
  x  <- all_trajectories$time_values
  y  <- all_trajectories$traj_values
  
  # Re-index individuals to 1..J within this subset
  id_reindexed <- as.integer(factor(id))
  
  data_list <- list(
    N  = nrow(all_trajectories),
    J  = max(id_reindexed),
    id = id_reindexed,
    x  = x,
    y  = y,
    sd = prior_sd
  )
  
  stan(
    file    = stan_file,
    data    = data_list,
    chains  = 4,
    iter    = 4000,
    warmup  = 2000,
    seed    = seed,
    refresh = 0,
    cores   = 4
  )
}




# ── Simulate one dataset ──────────────────────────────────────────────────────

set.seed(42)
dat              <- generate_data(LOD, L_0)
all_trajectories <- dat$all_trajectories
pop_params       <- dat$pop_params

cat("True population parameters:\n")
print(pop_params)

# ── Split trajectories by number of positive observations ────────────────────

sample_counts <- tapply(all_trajectories$index_values,
                        all_trajectories$index_values,
                        length)

cat(sprintf("\nTotal trajectories: %d\n", length(sample_counts)))
cat("Distribution of sample counts:\n")
print(table(sample_counts))

dir.create("throughput/by_sample_count", showWarnings = FALSE, recursive = TRUE)

# ── Fit BHM for each group ────────────────────────────────────────────────────

for (n_samples in 1:4) {
  target_ids <- as.integer(names(sample_counts[sample_counts == n_samples]))
  
  if (length(target_ids) == 0) {
    cat(sprintf("\nGroup n=%d: empty — skipping.\n", n_samples))
    next
  }
  
  cat(sprintf("\nGroup n=%d: %d trajectories — fitting model...\n",
              n_samples, length(target_ids)))
  
  group_traj <- all_trajectories[all_trajectories$index_values %in% target_ids, ]
  
  fit <- fit_model(group_traj, PRIOR_SD, STAN_FILE, seed = n_samples)
  
  out_path <- sprintf("throughput/by_sample_count/fit_n%d_samples.rds", n_samples)
  saveRDS(fit, out_path)
  cat(sprintf("  Saved to %s\n", out_path))
}

cat("\nDone.\n")