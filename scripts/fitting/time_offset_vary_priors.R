source("scripts/pooling/find_infx.R")

source("src/viral_bhm_vary_priors.R")

post <- extract(fit)

alpha_rmse <- sqrt(mean((post$mu_alpha - LOD)^2))
beta1_rmse <- sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2))
beta2_rmse <- sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2))
psi_rmse   <- sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
rmse_before_offset <- data.frame(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)

for(infx in 1:length(unique(all_trajectories$index_values))){
  offset <- (LOD - mean(post$alpha[, infx])) / mean(post$beta1[, infx])
  infx_indices <- which(all_trajectories$index_values == infx)
  all_trajectories$time_values[infx_indices] <- all_trajectories$time_values[infx_indices] + (-offset)
}

source("src/viral_bhm_vary_priors.R")

post2 <- extract(fit)

alpha_rmse <- sqrt(mean((post2$mu_alpha - LOD)^2))
beta1_rmse <- sqrt(mean((post2$mu_beta1 - pop_params$mu_beta1)^2))
beta2_rmse <- sqrt(mean((post2$mu_beta2 - pop_params$mu_beta2)^2))
psi_rmse   <- sqrt(mean((post2$mu_psi   - pop_params$mu_psi)^2))
rmse_after_offset <- data.frame(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)

