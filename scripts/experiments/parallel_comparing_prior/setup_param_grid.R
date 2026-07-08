# Run this once before submitting the job array.
# Writes param_grid.csv which the array jobs index into.

library(readr)

param_grid <- expand.grid(
  prior_sd = seq(from = 0.1, to = 10, by = 0.2),
  l_0      = 1:15
)

dir.create("throughput/results", showWarnings = FALSE, recursive = TRUE)

write_csv(param_grid, "throughput/param_grid.csv")
cat(sprintf("Wrote %d conditions to throughput/param_grid.csv\n", nrow(param_grid)))
