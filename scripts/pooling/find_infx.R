test_times <- seq(from = 1, to = ncol(matrix), by = l_0)

traj_values <- c()
index_values <- c()
time_values <- c()
person_values <- c()
date_values <- c()
positive_rows <- 1

for(row_ind in 1:nrow(matrix)){
  if(any(matrix[row_ind, test_times] > 0)){
    row <- matrix[row_ind,test_times]
    traj_indices <- which(row > 0)
    traj <- unlist(row[traj_indices])
    traj_values <- append(traj_values, traj)
    index_values <- append(index_values, rep(positive_rows, length(traj)))
    time_values <- append(time_values, (seq(1, length(traj)) - 1) * l_0)
    person_values <- append(person_values, rep(row_ind, length(traj)))
    #date_values <- append(date_values, traj_indices * l_0)
    date_values <- append(date_values, test_times[traj_indices]) # claude
    
    positive_rows <- positive_rows + 1
  }
}

all_trajectories <- data.frame(index_values, time_values, traj_values, person_values,
                               date_values)