# Run this after all array jobs have finished.
# Combines the 750 per-task CSVs into a single results file.

library(readr)
library(dplyr)

result_files <- list.files("throughput/results", pattern = "result_task_.*\\.csv",
                            full.names = TRUE)

cat(sprintf("Found %d result files (expected 750)\n", length(result_files)))

# Flag any missing tasks
all_ids      <- 1:750
found_ids    <- as.integer(gsub(".*task_(\\d+)\\.csv", "\\1", result_files))
missing_ids  <- setdiff(all_ids, found_ids)

if (length(missing_ids) > 0) {
  cat(sprintf("WARNING: %d tasks missing: %s\n",
              length(missing_ids), paste(missing_ids, collapse = ", ")))
} else {
  cat("All 750 tasks accounted for.\n")
}

rmse_list <- bind_rows(lapply(result_files, read_csv, show_col_types = FALSE))
write_csv(rmse_list, "throughput/rmse_list_vary_priors.csv")
cat(sprintf("Combined %d rows written to throughput/rmse_list_vary_priors.csv\n",
            nrow(rmse_list)))
