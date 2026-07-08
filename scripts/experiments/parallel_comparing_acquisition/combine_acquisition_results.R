# Run after all array jobs finish.
# Combines per-task CSVs into a single results file.

library(readr)
library(dplyr)

param_grid   <- read_csv("throughput/acquisition_param_grid.csv", show_col_types = FALSE)
n_conditions <- nrow(param_grid)

result_files <- list.files("throughput/results_acquisition",
                           pattern = "result_task_.*\\.csv",
                           full.names = TRUE)

cat(sprintf("Found %d result files (expected %d)\n", length(result_files), n_conditions))

found_ids   <- as.integer(gsub(".*task_(\\d+)\\.csv", "\\1", result_files))
missing_ids <- setdiff(seq_len(n_conditions), found_ids)

if (length(missing_ids) > 0) {
  cat(sprintf("WARNING: %d tasks missing: %s\n",
              length(missing_ids), paste(missing_ids, collapse = ", ")))
} else {
  cat("All tasks accounted for.\n")
}

rmse_list <- bind_rows(lapply(result_files, read_csv, show_col_types = FALSE))
write_csv(rmse_list, "throughput/rmse_list_acquisition.csv")
cat(sprintf("Combined %d rows written to throughput/rmse_list_acquisition.csv\n",
            nrow(rmse_list)))