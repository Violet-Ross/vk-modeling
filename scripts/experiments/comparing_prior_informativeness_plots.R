library(patchwork)
library(latex2exp)
library(tidyverse)

## read in prior informativeness results

rmse_list_full <- read_csv("throughput/rmse_list_vary_priors.csv")

mu_beta1    <- 2;   mu_beta2 <- -4;  mu_psi <- 5
sigma_beta1 <- 0.1; sigma_beta2 <- 0.1; sigma_psi <- 0.1
size <- 100000

beta1 <- rnorm(size, mean = mu_beta1, sd = sigma_beta1)
beta2 <- rnorm(size, mean = mu_beta2, sd = sigma_beta2)
psi <- rnorm(size, mean = mu_psi, sd = sigma_psi)

l_i <- psi * (1 - (beta1 / beta2))

true_l_0 <- quantile(l_i, 0.05, names = F)

total_rmse_inform <- rmse_list_full %>%
  mutate(total_rmse = alpha_rmse + beta1_rmse + beta2_rmse + psi_rmse) 

## read in acquisition results
acquisition_results <- read_csv("/Users/violetross/Desktop/FTAS/vk-modeling/throughput/rmse_list_acquisition.csv")

total_rmse_acquire <- acquisition_results %>%
  mutate(total_rmse = alpha_rmse + beta1_rmse + beta2_rmse + psi_rmse)

## plots of prior informativeness results
p0 <- rmse_list_full %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = alpha_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(y = TeX("$\\l_0$"), title = TeX("$\\alpha$"), x = "std dev of prior")

p1 <- rmse_list_full %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = beta1_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_1$"), x = "std dev of prior")

p2 <- rmse_list_full %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = beta2_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_2$"), x = "std dev of prior")

p3 <- rmse_list_full %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = psi_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\psi$"), fill = "RMSE", x = "std dev of prior")

p0 | p1 | p2 | p3

ggsave("fig/rmse_heatmaps_full.png", 
       plot = (p0 | p1 | p2 | p3) & theme(plot.margin = margin(5, 15, 5, 15)), 
       width = 15, height = 4, dpi = 300)

p_sum <- total_rmse_inform %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = total_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(y = TeX("$\\l_0$"), x = "std dev of prior", fill = "total RMSE",
       title = "Summed RMSE")

p_sum
ggsave("fig/rmse_heatmap_total.png", plot = p_sum, width = 5, height = 4, dpi = 300)


## plots of acquisition results


p0 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = alpha_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(y = "extra points", title = TeX("$\\alpha$"), x = "std dev of prior")

p1 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = beta1_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_1$"), x = "std dev of prior")

p2 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = beta2_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_2$"), x = "std dev of prior")

p3 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = psi_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_inform$total_rmse,
                                                                            total_rmse_acquire))) +  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\psi$"), fill = "RMSE", x = "std dev of prior")

p0 | p1 | p2 | p3

ggsave("fig/rmse_acquisition.png", 
       plot = (p0 | p1 | p2 | p3) & theme(plot.margin = margin(5, 15, 5, 15)), 
       width = 15, height = 4, dpi = 300)

p_sum <-  total_rmse_acquire %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = total_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse, total_rmse_acquire$total_rmse))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(y = "extra points", x = "std dev of prior", fill = "total RMSE",
       title = "Summed RMSE")

p_sum

ggsave("fig/rmse_acquisition_total.png", plot = p_sum, width = 5, height = 4, dpi = 300)
