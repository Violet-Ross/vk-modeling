# fit_by_sample_count.R
# Fix prior_sd = 9. For each l_0 in {3, 4, 5}, simulate one dataset, split
# trajectories by the number of observed positive samples, then fit a two-stage
# BHM to each non-empty group.  Full Stan fit objects are saved as .rds files.

library(rstan)
library(readr)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

setwd("/Users/violetross/Desktop/FTAS/vk-modeling")

# ── Constants ─────────────────────────────────────────────────────────────────

LOD       <- 3
L_0_VALUES <- c(3, 4, 5)
PRIOR_SD  <- 9
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


two_stage_fit <- function(all_trajectories, prior_sd, stan_file,
                          LOD, seed = 1) {
  # Stage 1
  fit1  <- fit_model(all_trajectories, prior_sd, stan_file, seed = seed)
  post1 <- rstan::extract(fit1)
  
  # Time offset correction
  # post1$alpha and post1$beta1 are indexed by individual after re-indexing,
  # so we iterate over the re-indexed IDs that were actually passed to fit1
  id_reindexed <- as.integer(factor(all_trajectories$index_values))
  n_infx       <- max(id_reindexed)
  
  original_ids <- vapply(seq_len(n_infx), function(i) {
    all_trajectories$index_values[which(id_reindexed == i)[1]]
  }, double(1))
  
  offsets <- numeric(n_infx)
  for (infx in seq_len(n_infx)) {
    offsets[infx]  <- (LOD - mean(post1$alpha[, infx])) / mean(post1$beta1[, infx])
    infx_indices   <- which(id_reindexed == infx)
    all_trajectories$time_values[infx_indices] <-
      all_trajectories$time_values[infx_indices] - offsets[infx]
  }
  
  offsets_df <- data.frame(
    index_value = original_ids,
    offset      = offsets
  )
  
  # Stage 2
  fit2 <- fit_model(all_trajectories, prior_sd, stan_file, seed = seed + 1000)
  
  list(fit = fit2, offsets = offsets_df)
}

# ── Outer loop over l_0 values ────────────────────────────────────────────────

for (L_0 in L_0_VALUES) {
  
  cat(sprintf("\n══════════════════════════════════════════\n"))
  cat(sprintf("  l_0 = %d\n", L_0))
  cat(sprintf("══════════════════════════════════════════\n"))
  
  # ── Simulate one dataset ────────────────────────────────────────────────────
  
  set.seed(42)
  dat              <- generate_data(LOD, L_0)
  all_trajectories <- dat$all_trajectories
  pop_params       <- dat$pop_params
  
  cat("True population parameters:\n")
  print(pop_params)
  
  # ── Split trajectories by number of positive observations ──────────────────
  
  sample_counts <- tapply(all_trajectories$index_values,
                          all_trajectories$index_values,
                          length)
  
  cat(sprintf("\nTotal trajectories: %d\n", length(sample_counts)))
  cat("Distribution of sample counts:\n")
  print(table(sample_counts))
  
  out_dir <- sprintf("throughput/by_sample_count/l0_%d", L_0)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ── Fit BHM for each group ──────────────────────────────────────────────────
  
  # Summary table to record infection counts per group
  group_summary <- data.frame(
    n_samples    = integer(),
    n_infections = integer(),
    stringsAsFactors = FALSE
  )
  
  for (n_samples in 1:4) {
    target_ids <- as.integer(names(sample_counts[sample_counts == n_samples]))
    
    if (length(target_ids) == 0) {
      cat(sprintf("\nGroup n=%d: empty — skipping.\n", n_samples))
      next
    }
    
    group_traj <- all_trajectories[all_trajectories$index_values %in% target_ids, ]
    
    # Count the number of infections (unique individuals) in this group
    n_infections <- length(unique(group_traj$index_values))
    
    cat(sprintf("\nGroup n=%d: %d trajectories, %d infections — fitting model...\n",
                n_samples, length(target_ids), n_infections))
    
    # Record in summary
    group_summary <- rbind(group_summary, data.frame(
      n_samples    = n_samples,
      n_infections = n_infections,
      stringsAsFactors = FALSE
    ))
    
    result <- two_stage_fit(group_traj, PRIOR_SD, STAN_FILE, LOD, seed = n_samples)
    
    fit_path     <- sprintf("%s/fit_n%d_samples.rds", out_dir, n_samples)
    offsets_path <- sprintf("%s/offsets_n%d_samples.csv", out_dir, n_samples)
    
    saveRDS(result$fit, fit_path)
    write_csv(result$offsets, offsets_path)
    
    cat(sprintf("  Saved fit to %s\n", fit_path))
    cat(sprintf("  Saved offsets to %s\n", offsets_path))
  }
  
  # ── Save and display infection count summary ────────────────────────────────
  
  summary_path <- sprintf("%s/infection_counts_by_group.csv", out_dir)
  write_csv(group_summary, summary_path)
  
  cat(sprintf("\nInfection counts by group (l_0 = %d):\n", L_0))
  print(group_summary)
  cat(sprintf("Summary saved to %s\n", summary_path))
}

cat("\nDone.\n")
