LOD <- 3
improvement_list <- data.frame(alpha_rmse_improvement = numeric(),
                               beta1_rmse_improvement = numeric(),
                               beta2_rmse_improvement = numeric(),
                               psi_rmse_improvement = numeric())

for(l_0 in 1:5){
  source("scripts/data_gen/new_matrix.R")
  source("scripts/fitting/time_offset.R")
  
  improvement <- rmse_before_offset - rmse_after_offset
  
  improvement_list <- rbind(improvement_list, improvement)
}
