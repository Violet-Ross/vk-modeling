# Run after all array jobs finish.
# Combines per-task CSVs into a single results file.

library(readr)
library(dplyr)

n_conditions <- 150  # 50 prior_sd x 3 l_0

result_files <- list.files("throughput/results_l0_restricted",
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

results <- bind_rows(lapply(result_files, read_csv, show_col_types = FALSE))
write_csv(results, "throughput/results_l0_restricted.csv")
cat(sprintf("Combined %d rows written to throughput/results_l0_restricted.csv\n",
            nrow(results)))