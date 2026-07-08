compute_entropy <- function(x_star, infection_num) {
  alpha <- post$alpha[, infection_num]
  beta1 <- post$beta1[, infection_num]
  beta2 <- post$beta2[, infection_num]
  psi   <- post$psi[,   infection_num]
  
  mu <- alpha + beta1 * x_star + beta2 * pmax(x_star - psi, 0)
  
  dens <- density(mu)
  pk   <- dens$y / sum(dens$y)
  -sum(pk * log(pk + 1e-10))
}

test_points <- function(person_id, x_vals){
  person <- matrix[person_id,]
  infection_start <- all_trajectories[all_trajectories$person_values == person_id, 5][1]
  y_vals <- unlist(matrix[person_id, (infection_start + x_vals)])
  
  return(y_vals)
}

entropy_maximizer <- function(infection_num, N = 1){
  infx <- all_trajectories %>% filter(index_values == infection_num)
  
  first_interval_hit <- min(infx$time_values)
  final_interval_hit <- max(infx$time_values)
  
  closest_time_to_0 <- first_interval_hit - floor(first_interval_hit)

  beta1 <- post$beta1[, infection_num]
  beta2 <- post$beta2[, infection_num]
  psi   <- post$psi[,   infection_num]
  closest_time_to_end <- mean((beta2 * psi) / (beta1 + beta2))
  
  if(closest_time_to_0 >= closest_time_to_end){
    time_options <- seq(closest_time_to_end, closest_time_to_0, 1)
  }
  else{
    time_options <- seq(closest_time_to_0, closest_time_to_end, 1)
  }
  
  
  time_options <- time_options[!(time_options %in% infx$time_values)]
  
  max_entropy <- -Inf
  for(option in time_options){
    entropy <- compute_entropy(option, infection_num)
    if(entropy > max_entropy){
      max_entropy <- entropy
      max_option <- option
    }
  }
  
  actual_index <- max_option - first_interval_hit
  return(c(actual_index, first_interval_hit))
}

acquisition_func <- function(infx_id){
  results <- entropy_maximizer(infx_id)
  samples_x <- results[1]
  first_interval_hit <- results[2]
  person_id  <- all_trajectories %>%
    filter(index_values == infx_id) %>%
    pull(4) %>%
    first()
  
  samples_y <- test_points(person_id, samples_x)
  
  # Guard: skip if test_points returned nothing (out-of-range column index)
  if(length(samples_y) == 0){
    warning(paste("test_points returned empty for infx_id =", infx_id,
                  "samples_x =", samples_x))
    return(invisible(NULL))
  }
  
  for(j in 1:length(samples_x)){
    new_row <- data.frame(index_values  = infx_id,
                          time_values   = samples_x[j] + first_interval_hit,
                          traj_values   = samples_y[j],
                          person_values = person_id,
                          date_values   = NA)
    fitting_data <<- rbind(fitting_data, new_row)
  }
}