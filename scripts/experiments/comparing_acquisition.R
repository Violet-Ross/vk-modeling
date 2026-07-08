LOD <- 3

#determine l_0
size <- 10000 

alpha <- LOD
mu_beta1 <- 2
mu_beta2 <- -4
mu_psi <- 5

sigma_beta1 <- 0.1
sigma_beta2 <- 0.1
sigma_psi <- 0.1

beta1 <- rnorm(size, mean = mu_beta1, sd = sigma_beta1)
beta2 <- rnorm(size, mean = mu_beta2, sd = sigma_beta2)
psi <- rnorm(size, mean = mu_psi, sd = sigma_psi)

l_i <- psi * (1 - (beta1 / beta2))

l_0 <- quantile(l_i, 0.05, names = F)

E.l_i <- mean(l_i)

summation <- 0
for(j in 0:(l_0 - 1)){
  summation <- summation + ceiling((l_i - j) / l_0) * (1 / l_0)
}
E.N <- mean(summation)

E.remaining_pts <- E.l_i - E.N

extract <- rstan::extract
rmse_list <- data.frame(alpha_rmse = numeric(),
                        beta1_rmse = numeric(),
                        beta2_rmse = numeric(),
                        psi_rmse = numeric(),
                        sd = numeric(),
                        l0 = numeric())

total_start <- proc.time()

for(prior_sd in seq(from = 0.1, to = 10, by = 0.2)){
  for(extra_pts in 1:floor(E.remaining_pts)){
    print(paste("sd =", prior_sd, "| extra_pts =", extra_pts))
    for(iter in 1:5){
      source("scripts/data_gen/new_matrix_vary_priors.R")
      source("scripts/fitting/time_offset_acquisition.R")
      
      rmse_row <- c(unlist(rmse_after_offset[1,], use.names = F), prior_sd, extra_pts)
      rmse_list <- rbind(rmse_list, rmse_row)
    }
  }
}

colnames(rmse_list) <- c("alpha_rmse",
                         "beta1_rmse",
                         "beta2_rmse",
                         "psi_rmse",
                         "sd",
                         "extra_pts")

write_csv(rmse_list, "throughput/rmse_list_acquisition.csv")

print(proc.time() - total_start)