# Run once before submitting the job array.
# Writes param_grid.csv which the array jobs index into.
# Restricted to l_0 in {3, 4, 5}.

library(readr)

param_grid <- expand.grid(
  prior_sd = seq(from = 0.1, to = 10, by = 0.2),
  l_0      = c(3, 4, 5)
)

dir.create("throughput/results_l0_restricted", showWarnings = FALSE, recursive = TRUE)

write_csv(param_grid, "throughput/param_grid_l0_restricted.csv")
cat(sprintf("Wrote %d conditions to throughput/param_grid_l0_restricted.csv\n", nrow(param_grid)))