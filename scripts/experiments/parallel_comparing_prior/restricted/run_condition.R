# Worker script — run once per (prior_sd, l_0) condition.
# l_0 restricted to {3, 4, 5}.
# Stores posterior means of population-level parameters alongside RMSE.

library(rstan)
library(readr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# ── Constants ─────────────────────────────────────────────────────────────────

LOD       <- 3
N_ITER    <- 5
STAN_FILE <- "src/viral_bhm_vary_priors.stan"

# True population-level parameter values
TRUE_ALPHA <- LOD
TRUE_BETA1 <- 2
TRUE_BETA2 <- -4
TRUE_PSI   <- 5

# ── Read condition for this task ──────────────────────────────────────────────

task_id    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid <- read_csv("throughput/param_grid_l0_restricted.csv", show_col_types = FALSE)
prior_sd   <- param_grid$prior_sd[task_id]
l_0        <- param_grid$l_0[task_id]

cat(sprintf("Task %d | prior_sd = %.2f | l_0 = %d\n", task_id, prior_sd, l_0))

# ── Functions ─────────────────────────────────────────────────────────────────

generate_data <- function(LOD, l_0) {
  num_people    <- 500
  num_days      <- 500
  prev          <- 0.3
  min_x         <- 0
  max_x         <- 50
  
  alpha       <- LOD
  mu_beta1    <- 2;   mu_beta2 <- -4;  mu_psi <- 5
  sigma_beta1 <- 0.1; sigma_beta2 <- 0.1; sigma_psi <- 0.1
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
  
  data_list <- list(
    N  = nrow(all_trajectories),
    J  = max(id),
    id = id,
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


compute_results <- function(fit, pop_params, LOD) {
  post <- rstan::extract(fit)
  
  # Posterior means of population-level parameters
  est_alpha <- mean(post$mu_alpha)
  est_beta1 <- mean(post$mu_beta1)
  est_beta2 <- mean(post$mu_beta2)
  est_psi   <- mean(post$mu_psi)
  
  data.frame(
    # Estimated parameters
    est_alpha  = est_alpha,
    est_beta1  = est_beta1,
    est_beta2  = est_beta2,
    est_psi    = est_psi,
    # True parameters
    true_alpha = LOD,
    true_beta1 = pop_params$mu_beta1,
    true_beta2 = pop_params$mu_beta2,
    true_psi   = pop_params$mu_psi,
    # RMSE
    alpha_rmse = sqrt(mean((post$mu_alpha - LOD)^2)),
    beta1_rmse = sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2)),
    beta2_rmse = sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2)),
    psi_rmse   = sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
  )
}

# ── Run N_ITER replicates for this condition ──────────────────────────────────

rows <- vector("list", N_ITER)

for (iter in seq_len(N_ITER)) {
  cat(sprintf("  iter %d / %d\n", iter, N_ITER))
  
  dat              <- generate_data(LOD, l_0)
  all_trajectories <- dat$all_trajectories
  pop_params       <- dat$pop_params
  
  fit1  <- fit_model(all_trajectories, prior_sd, STAN_FILE, seed = iter)
  post1 <- rstan::extract(fit1)
  
  # Time offset correction
  for (infx in seq_len(max(all_trajectories$index_values))) {
    offset       <- (LOD - mean(post1$alpha[, infx])) / mean(post1$beta1[, infx])
    infx_indices <- which(all_trajectories$index_values == infx)
    all_trajectories$time_values[infx_indices] <-
      all_trajectories$time_values[infx_indices] - offset
  }
  
  fit2    <- fit_model(all_trajectories, prior_sd, STAN_FILE, seed = iter + 1000)
  results <- compute_results(fit2, pop_params, LOD)
  
  rows[[iter]] <- data.frame(
    results,
    sd      = prior_sd,
    l0      = l_0,
    task_id = task_id,
    iter    = iter
  )
}

# ── Write this task's results ─────────────────────────────────────────────────

result   <- do.call(rbind, rows)
out_path <- sprintf("throughput/results_l0_restricted/result_task_%04d.csv", task_id)
write_csv(result, out_path)
cat(sprintf("Written to %s\n", out_path))