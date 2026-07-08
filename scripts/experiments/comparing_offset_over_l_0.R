LOD <- 3
extract <- rstan::extract
l_0 <- 4
prior_sd <- 0.5
extra_pts <- 2

source("scripts/data_gen/new_matrix_vary_priors.R")

source("scripts/fitting/time_offset_acquisition.R")

improvement_list <- data.frame(alpha_rmse_improvement = numeric(),
                               beta1_rmse_improvement = numeric(),
                               beta2_rmse_improvement = numeric(),
                               psi_rmse_improvement = numeric())

improvement <- rmse_before_offset - rmse_after_offset

improvement_list <- rbind(improvement_list, improvement)

for(l_0 in 1:5){
  source("scripts/data_gen/new_matrix_vary_priors.R")
  source("scripts/fitting/time_offset_vary_priors.R")
  
  improvement <- rmse_before_offset - rmse_after_offset
  
  improvement_list <- rbind(improvement_list, improvement)
}

write_csv(improvement_list, "throughput/improvement_list.csv")


read_csv("throughput/improvement_list.csv")
