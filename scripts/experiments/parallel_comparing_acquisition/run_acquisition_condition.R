# Worker script — one SLURM array task per (prior_sd, extra_pts) condition.

library(rstan)
library(tidyverse)
library(readr)

rstan_options(auto_write = TRUE)

# ── Constants ─────────────────────────────────────────────────────────────────

LOD       <- 3
N_ITER    <- 5
STAN_FILE <- "src/viral_bhm_vary_priors.stan"

# ── Read condition and fixed constants for this task ──────────────────────────

task_id    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
param_grid <- read_csv("throughput/acquisition_param_grid.csv", show_col_types = FALSE)
constants  <- readRDS("throughput/acquisition_constants.rds")

prior_sd  <- param_grid$prior_sd[task_id]
extra_pts <- param_grid$extra_pts[task_id]
l_0       <- constants$l_0

cat(sprintf("Task %d | prior_sd = %.2f | extra_pts = %d\n", task_id, prior_sd, extra_pts))

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
  
  list(all_trajectories = all_trajectories,
       pop_params       = pop_params,
       matrix           = mat)
}


fit_model <- function(data_df, prior_sd, stan_file, seed = 123) {
  data_list <- list(
    N  = nrow(data_df),
    J  = max(data_df[[1]]),
    id = data_df[[1]],
    x  = data_df[[2]],
    y  = data_df[[3]],
    sd = prior_sd
  )
  capture.output(
    fit <- stan(file    = stan_file,
                data    = data_list,
                chains  = 4,
                iter    = 4000,
                warmup  = 2000,
                seed    = seed,
                refresh = 0,
                cores   = 4)
  )
  fit
}


compute_rmse <- function(post, pop_params, LOD) {
  data.frame(
    alpha_rmse = sqrt(mean((post$mu_alpha - LOD)^2)),
    beta1_rmse = sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2)),
    beta2_rmse = sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2)),
    psi_rmse   = sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
  )
}


make_entropy_fns <- function(post, all_trajectories, mat) {
  
  compute_entropy <- function(x_star, infection_num) {
    alpha <- post$alpha[, infection_num]
    beta1 <- post$beta1[, infection_num]
    beta2 <- post$beta2[, infection_num]
    psi   <- post$psi[,   infection_num]
    mu    <- alpha + beta1 * x_star + beta2 * pmax(x_star - psi, 0)
    dens  <- density(mu)
    pk    <- dens$y / sum(dens$y)
    -sum(pk * log(pk + 1e-10))
  }
  
  test_points <- function(person_id, x_vals) {
    infection_start <- all_trajectories[all_trajectories$person_values == person_id, 5][1]
    unlist(mat[person_id, (infection_start + x_vals)])
  }
  
  entropy_maximizer <- function(infection_num) {
    infx               <- all_trajectories %>% filter(index_values == infection_num)
    first_interval_hit <- min(infx$time_values)
    final_interval_hit <- max(infx$time_values)
    closest_time_to_0  <- first_interval_hit - floor(first_interval_hit)
    
    beta1 <- post$beta1[, infection_num]
    beta2 <- post$beta2[, infection_num]
    psi   <- post$psi[,   infection_num]
    closest_time_to_end <- mean((beta2 * psi) / (beta1 + beta2))
    
    time_options <- seq(min(closest_time_to_0, closest_time_to_end),
                        max(closest_time_to_0, closest_time_to_end), 1)
    
    time_options <- time_options[!(time_options %in% infx$time_values)]
    
    if (length(time_options) == 0) return(c(NA, first_interval_hit))
    
    max_entropy <- -Inf; max_option <- NA
    for (option in time_options) {
      entropy <- compute_entropy(option, infection_num)
      if (entropy > max_entropy) { max_entropy <- entropy; max_option <- option }
    }
    c(max_option, first_interval_hit)
  }
  
  list(entropy_maximizer = entropy_maximizer, test_points = test_points)
}


run_condition <- function(prior_sd, extra_pts, LOD, l_0, stan_file, n_iter = 5) {
  rows <- vector("list", n_iter)
  
  for (iter in seq_len(n_iter)) {
    cat(sprintf("  iter %d / %d\n", iter, n_iter))
    
    dat              <- generate_data(LOD, l_0)
    all_trajectories <- dat$all_trajectories
    pop_params       <- dat$pop_params
    mat              <- dat$matrix
    
    # ── Fit 1: before offset ──────────────────────────────────────────────────
    fit1  <- fit_model(all_trajectories, prior_sd, stan_file, seed = iter)
    post1 <- rstan::extract(fit1)
    
    # ── Time offset correction ────────────────────────────────────────────────
    for (infx in seq_len(max(all_trajectories$index_values))) {
      offset       <- (LOD - mean(post1$alpha[, infx])) / mean(post1$beta1[, infx])
      infx_indices <- which(all_trajectories$index_values == infx)
      all_trajectories$time_values[infx_indices] <-
        all_trajectories$time_values[infx_indices] - offset
    }
    
    # ── Fit 2: after offset ───────────────────────────────────────────────────
    fitting_data <- all_trajectories
    fit2         <- fit_model(fitting_data, prior_sd, stan_file, seed = iter + 1000)
    post2        <- rstan::extract(fit2)
    
    # ── Acquisition loop ──────────────────────────────────────────────────────
    for (step in seq_len(extra_pts)) {
      cat(sprintf("    acquisition step %d / %d\n", step, extra_pts))
      
      fns <- make_entropy_fns(post2, all_trajectories, mat)
      
      for (infx_id in unique(all_trajectories$index_values)) {
        results            <- fns$entropy_maximizer(infx_id)
        samples_x          <- results[1]
        first_interval_hit <- results[2]
        
        if (is.na(samples_x)) next
        
        person_id <- all_trajectories %>%
          filter(index_values == infx_id) %>% pull(4) %>% first()
        
        samples_y <- fns$test_points(person_id, samples_x)
        
        if (length(samples_y) == 0) {
          warning(paste("test_points returned empty for infx_id =", infx_id,
                        "samples_x =", samples_x))
          next
        }
        
        for (j in seq_along(samples_x)) {
          fitting_data <- rbind(fitting_data,
                                data.frame(index_values  = infx_id,
                                           time_values   = samples_x[j] + first_interval_hit,
                                           traj_values   = samples_y[j],
                                           person_values = person_id,
                                           date_values   = NA))
        }
      }
      
      # Refit on augmented data
      all_trajectories <- fitting_data
      fit2  <- fit_model(fitting_data, prior_sd, stan_file, seed = iter + step * 1000)
      post2 <- rstan::extract(fit2)
    }
    
    # ── Record RMSE after final acquisition fit ───────────────────────────────
    rmse <- compute_rmse(post2, pop_params, LOD)
    rows[[iter]] <- data.frame(
      alpha_rmse = rmse$alpha_rmse,
      beta1_rmse = rmse$beta1_rmse,
      beta2_rmse = rmse$beta2_rmse,
      psi_rmse   = rmse$psi_rmse,
      sd         = prior_sd,
      extra_pts  = extra_pts,
      task_id    = task_id,
      iter       = iter
    )
  }
  
  do.call(rbind, rows)
}

# ── Run and save ──────────────────────────────────────────────────────────────

result   <- run_condition(prior_sd, extra_pts, LOD, l_0, STAN_FILE)
out_path <- sprintf("throughput/results_acquisition/result_task_%04d.csv", task_id)
write_csv(result, out_path)
cat(sprintf("Written to %s\n", out_path))