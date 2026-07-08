# plot_by_sample_count.R
# Reads saved .rds fits and infection count CSVs from fit_by_sample_count.R,
# then produces three plots per l_0 value:
#   - l0{N}_error_by_n.png    : RMSE by group, faceted by parameter
#   - l0{N}_count_by_n.png   : infection count per group
#   - l0{N}_posteriors.png   : posterior histograms for the n=3 group

library(rstan)
library(readr)
library(tidyverse)
library(patchwork)

setwd("/Users/violetross/Desktop/FTAS/vk-modeling")

# -- Constants -----------------------------------------------------------------

LOD        <- 3
L_0_VALUES <- c(3, 4, 5)

# True population parameters (fixed across all simulations)
pop_params <- data.frame(mu_beta1 = 2, mu_beta2 = -4, mu_psi = 5)

# -- Functions -----------------------------------------------------------------

compute_rmse <- function(fit, pop_params, LOD) {
  post <- rstan::extract(fit)
  data.frame(
    alpha_rmse = sqrt(mean((post$mu_alpha - LOD)^2)),
    beta1_rmse = sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2)),
    beta2_rmse = sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2)),
    psi_rmse   = sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
  )
}


make_plots <- function(fits, group_summary, pop_params, LOD, l_0, fig_dir) {
  
  # -- RMSE bar chart per parameter, faceted by n_samples --------------------
  rmse_rows <- lapply(names(fits), function(n) {
    rmse <- compute_rmse(fits[[n]], pop_params, LOD)
    rmse$num_samples <- as.integer(n)
    rmse
  })
  rmse_df <- do.call(rbind, rmse_rows)
  
  p_rmse <- rmse_df %>%
    pivot_longer(cols = c(beta1_rmse, beta2_rmse, psi_rmse),
                 names_to = "var", values_to = "rmse") %>%
    ggplot() +
    geom_col(aes(x = num_samples, y = rmse)) +
    facet_wrap(~var) +
    labs(title = sprintf("RMSE by group (l_0 = %d)", l_0),
         x = "Number of samples", y = "RMSE") +
    theme_minimal()
  
  ggsave(file.path(fig_dir, sprintf("l0%d_error_by_n.png", l_0)), p_rmse)
  
  # -- Infection count bar chart ---------------------------------------------
  p_count <- group_summary %>%
    ggplot() +
    geom_col(aes(x = n_samples, y = n_infections)) +
    labs(title = sprintf("Infections per group (l_0 = %d)", l_0),
         x = "Number of samples", y = "Number of infections") +
    theme_minimal()
  
  ggsave(file.path(fig_dir, sprintf("l0%d_count_by_n.png", l_0)), p_count)
  
  # -- Posterior histograms for n=3 group (if available) --------------------
  if ("3" %in% names(fits)) {
    post3 <- rstan::extract(fits[["3"]])
    p0 <- wrap_elements(~hist(post3$mu_beta1, breaks = 20,
                              main = "mu_beta1", xlab = ""))
    p1 <- wrap_elements(~hist(post3$mu_beta2, breaks = 20,
                              main = "mu_beta2", xlab = ""))
    p2 <- wrap_elements(~hist(post3$mu_psi,   breaks = 20,
                              main = "mu_psi",   xlab = ""))
    p_post <- p0 | p1 | p2
    ggsave(file.path(fig_dir, sprintf("l0%d_posteriors.png", l_0)),
           p_post, width = 12, height = 4)
  } else {
    cat(sprintf("  Note: n=3 group not available for l_0=%d; skipping posterior plot.\n", l_0))
  }
}

# -- Outer loop over l_0 values -----------------------------------------------

dir.create("fig", showWarnings = FALSE, recursive = TRUE)

for (L_0 in L_0_VALUES) {
  
  cat(sprintf("\n==========================================\n"))
  cat(sprintf("  Plotting l_0 = %d\n", L_0))
  cat(sprintf("==========================================\n"))
  
  in_dir <- sprintf("throughput/by_sample_count/l0_%d", L_0)
  
  # -- Load fits -------------------------------------------------------------
  fits <- list()
  for (n_samples in 2:4) {
    rds_path <- sprintf("%s/fit_n%d_samples.rds", in_dir, n_samples)
    if (file.exists(rds_path)) {
      fits[[as.character(n_samples)]] <- readRDS(rds_path)
      cat(sprintf("  Loaded %s\n", rds_path))
    } else {
      cat(sprintf("  %s not found -- skipping.\n", rds_path))
    }
  }
  
  if (length(fits) == 0) {
    cat(sprintf("  No fits found for l_0=%d; skipping.\n", L_0))
    next
  }
  
  # -- Load infection count summary ------------------------------------------
  summary_path  <- sprintf("%s/infection_counts_by_group.csv", in_dir)
  group_summary <- read_csv(summary_path, show_col_types = FALSE)
  
  # -- Produce plots ---------------------------------------------------------
  make_plots(fits, group_summary, pop_params, LOD, L_0, fig_dir = "fig")
  cat(sprintf("  Plots saved to fig/l0%d_*.png\n", L_0))
}

cat("\nDone.\n")



off2 <- read_csv("/Users/violetross/Desktop/FTAS/vk-modeling/throughput/by_sample_count/l0_4/offsets_n2_samples.csv")
off3 <- read_csv("/Users/violetross/Desktop/FTAS/vk-modeling/throughput/by_sample_count/l0_4/offsets_n3_samples.csv")
off4 <- read_csv("/Users/violetross/Desktop/FTAS/vk-modeling/throughput/by_sample_count/l0_4/offsets_n4_samples.csv")


