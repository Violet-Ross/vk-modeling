library(patchwork)
library(latex2exp)
library(tidyverse)

acquisition_results <- read_csv("/Users/violetross/Desktop/FTAS/vk-modeling/throughput/rmse_list_acquisition.csv")
mu_beta1    <- 2;   mu_beta2 <- -4;  mu_psi <- 5
sigma_beta1 <- 0.1; sigma_beta2 <- 0.1; sigma_psi <- 0.1


p0 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = alpha_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(y = "extra points", title = TeX("$\\alpha$"), x = "std dev of prior")

p1 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = beta1_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_1$"), x = "std dev of prior")

p2 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = beta2_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\beta_2$"), x = "std dev of prior")

p3 <- acquisition_results %>%
  ggplot() +
  geom_tile(aes(x = sd, y = extra_pts, fill = psi_rmse)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(rmse_list_full$alpha_rmse, rmse_list_full$beta1_rmse, rmse_list_full$beta2_rmse, rmse_list_full$psi_rmse))) +
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = TeX("$\\psi$"), fill = "RMSE", x = "std dev of prior")

p0 | p1 | p2 | p3

ggsave("fig/rmse_acquisition.png", 
       plot = (p0 | p1 | p2 | p3) & theme(plot.margin = margin(5, 15, 5, 15)), 
       width = 15, height = 4, dpi = 300)

total_rmse_acquire <- acquisition_results %>%
  mutate(total_rmse = alpha_rmse + beta1_rmse + beta2_rmse + psi_rmse)


p_sum <-  total_rmse %>%
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