LOD <- 3
extract <- rstan::extract
rmse_list <- data.frame(alpha_rmse = numeric(),
                        beta1_rmse = numeric(),
                        beta2_rmse = numeric(),
                        psi_rmse = numeric(),
                        sd = numeric(),
                        l0 = numeric())

total_start <- proc.time()

for(prior_sd in seq(from = 0.1, to = 10, by = 0.2)){
  for(l_0 in 1:15){
    print(paste("sd =", prior_sd, "| l_0 =", l_0))
    for(iter in 1:5){
      source("scripts/data_gen/new_matrix_vary_priors.R")
      source("scripts/fitting/time_offset_vary_priors.R")
      
      rmse_row <- c(unlist(rmse_after_offset[1,], use.names = F), prior_sd, l_0)
      rmse_list <- rbind(rmse_list, rmse_row)
    }
  }
  write_csv(rmse_list, "throughput/rmse_list_vary_priors.csv")
}

colnames(rmse_list) <- c("alpha_rmse",
                         "beta1_rmse",
                         "beta2_rmse",
                         "psi_rmse",
                         "sd",
                         "l0")

write_csv(rmse_list, "throughput/rmse_list_vary_priors.csv")

print(proc.time() - total_start)