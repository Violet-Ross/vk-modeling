# locate infections
print("finding infections")
source("scripts/pooling/find_infx.R")

# fit BHM to infections with incorrect time stamps
print("fitting BHM before offset")
source("src/viral_bhm_vary_priors.R")

post <- extract(fit)

alpha_rmse <- sqrt(mean((post$mu_alpha - LOD)^2))
beta1_rmse <- sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2))
beta2_rmse <- sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2))
psi_rmse   <- sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
rmse_before_offset <- data.frame(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)

# correct the time stamps
offset_list <- c()
for(infx in 1:length(unique(all_trajectories$index_values))){
  offset <- (LOD - mean(post$alpha[, infx])) / mean(post$beta1[, infx])
  offset_list <- c(offset_list, offset)
  infx_indices <- which(all_trajectories$index_values == infx)
  all_trajectories$time_values[infx_indices] <- all_trajectories$time_values[infx_indices] + (-offset)
}

# fit BHM to infections with corrected time stamps
print("fitting BHM after offset")
fitting_data <- all_trajectories
source("src/viral_bhm_acquisition.R")
post <- extract(fit)

alpha_rmse <- sqrt(mean((post$mu_alpha - LOD)^2))
beta1_rmse <- sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2))
beta2_rmse <- sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2))
psi_rmse   <- sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
rmse_after_offset <- data.frame(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)

source("scripts/fitting/entropy_functions.R")

for(point_adding_step in 1:extra_pts){
  print(paste0("adding point number ", point_adding_step, " to each infx"))
  for(infx_id in unique(all_trajectories$index_values)){
    acquisition_func(infx_id)
  } 
  
  all_trajectories <- fitting_data
  
  source("src/viral_bhm_acquisition.R")
  post <- extract(fit)
}

alpha_rmse <- sqrt(mean((post$mu_alpha - LOD)^2))
beta1_rmse <- sqrt(mean((post$mu_beta1 - pop_params$mu_beta1)^2))
beta2_rmse <- sqrt(mean((post$mu_beta2 - pop_params$mu_beta2)^2))
psi_rmse   <- sqrt(mean((post$mu_psi   - pop_params$mu_psi)^2))
rmse_after_acquisition <- data.frame(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)


