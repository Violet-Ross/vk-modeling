# Run once before submitting the job array.
# Computes E.remaining_pts from a large stable draw, then writes param_grid.csv.

library(readr)

set.seed(42)  # fix seed so grid is reproducible across runs

LOD <- 3

alpha    <- LOD
mu_beta1 <- 2;  mu_beta2 <- -4;  mu_psi <- 5
sigma_beta1 <- 0.1; sigma_beta2 <- 0.1; sigma_psi <- 0.1

size  <- 10000
beta1 <- rnorm(size, mean = mu_beta1, sd = sigma_beta1)
beta2 <- rnorm(size, mean = mu_beta2, sd = sigma_beta2)
psi   <- rnorm(size, mean = mu_psi,   sd = sigma_psi)

l_i <- psi * (1 - (beta1 / beta2))
l_0 <- quantile(l_i, 0.05, names = FALSE)
E.l_i <- mean(l_i)

summation <- 0
for (j in 0:(l_0 - 1)) {
  summation <- summation + ceiling((l_i - j) / l_0) * (1 / l_0)
}
E.N <- mean(summation)
E.remaining_pts <- E.l_i - E.N

max_extra_pts <- floor(E.remaining_pts)
cat(sprintf("l_0 = %.4f | E.remaining_pts = %.4f | max_extra_pts = %d\n",
            l_0, E.remaining_pts, max_extra_pts))

param_grid <- expand.grid(
  prior_sd   = seq(from = 0.1, to = 10, by = 0.2),
  extra_pts  = 1:max_extra_pts
)

dir.create("throughput/results_acquisition", showWarnings = FALSE, recursive = TRUE)

# Save both the grid and the fixed constants workers will need
write_csv(param_grid, "throughput/acquisition_param_grid.csv")
saveRDS(list(l_0 = l_0, E.remaining_pts = E.remaining_pts),
        "throughput/acquisition_constants.rds")

cat(sprintf("Wrote %d conditions to throughput/acquisition_param_grid.csv\n",
            nrow(param_grid)))