num_people <- 500
num_days <- 500

prev <- 0.3
min_x <- 0
max_x <- 50

size <- 10000

# initialize matrix of 0s
infection_matrix <- matrix(0, nrow = num_people, ncol = num_days)

# determine which rows will be infected
is_infx <- rbinom(num_people, size = 1, prev)

# population level parameters (fixed)
alpha <- LOD
mu_beta1 <- 2
mu_beta2 <- -4
mu_psi <- 5

sigma_beta1 <- 0.1
sigma_beta2 <- 0.1
sigma_psi <- 0.1

sigma_epsilon <- abs(rnorm(1, 0, sd = 0.5)) # observation noise, drawn from prior

pop_params <- data.frame(mu_beta1, mu_beta2, mu_psi, 
                         sigma_beta1, sigma_beta2, sigma_psi, sigma_epsilon)

for(row_ind in 1:nrow(infection_matrix)){
  if(is_infx[row_ind] == 1){
    # individual level parameters
    beta1 <- rnorm(1, mu_beta1, sigma_beta1)
    beta2 <- rnorm(1, mu_beta2, sigma_beta2)
    psi <- rnorm(1, mu_psi, sigma_psi)
    
    x_vals = seq(from = min_x, to = max_x, length.out = max_x - min_x)
    
    # observation level
    y_vals <- c()
    for(x in x_vals){
      obs_mu = alpha + (beta1 * x) + (beta2 * (x - psi) * (x > psi))
      y <- rnorm(1, obs_mu, sigma_epsilon)
      y_vals <- c(y_vals, y)
    } 
    pos_y_vals <- y_vals[y_vals > 0]
    infx_length <- length(pos_y_vals)
    
    start_day <- sample(1:(num_days - infx_length), 1)
    
    infection_matrix[row_ind, start_day:(start_day + infx_length - 1)] <- pos_y_vals
  }
}


matrix <- data.frame(infection_matrix)
