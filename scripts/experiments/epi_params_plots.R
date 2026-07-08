library(patchwork)
library(latex2exp)
library(tidyverse)

## read in onephase results
onephase_results <- read_csv("throughput/rmse_list_vary_priors_onephase.csv")

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

## get means over iterations
means <- onephase_results %>%
  group_by(sd, l0) %>%
  summarise(mean_psi_rmse = mean(psi_rmse), 
            mean_peak_rmse = mean(peak_rmse),
            mean_duration_rmse = mean(duration_rmse))


## plots of onephase epi results (with fixed colorbar)
p1 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_psi_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(means$mean_psi_rmse, means$mean_duration_rmse, means$mean_peak_rmse))) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Duration of Proliferation Phase", x = "std dev of prior")

p2 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_peak_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(means$mean_psi_rmse, means$mean_duration_rmse, means$mean_peak_rmse))) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Peak Viral Load", x = "std dev of prior")

p3 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_duration_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(c(means$mean_psi_rmse, means$mean_duration_rmse, means$mean_peak_rmse))) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Duration of Clearance Phase", fill = "RMSE", x = "std dev of prior")

p1 | p2 | p3

ggsave("fig/rmse_epi_fixed_colorscale.png", 
       plot = (p1 | p2 | p3) & theme(plot.margin = margin(5, 15, 5, 15)), 
       width = 13, height = 4, dpi = 300)

## plots of onephase epi results (with varying colorbars)
p1 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_psi_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(means$mean_psi_rmse)) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom", axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Duration of Proliferation Phase", x = "std dev of prior", fill = "RMSE") +
  guides(fill = guide_colorbar(barwidth = 10))

p2 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_peak_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = range(means$mean_peak_rmse)) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom", axis.title.y = element_blank(), axis.line = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Peak Viral Load", x = "std dev of prior", fill = "") +
  guides(fill = guide_colorbar(barwidth = 10))

p3 <- means %>%
  ggplot() +
  geom_tile(aes(x = sd, y = l0, fill = mean_duration_rmse)) +
  geom_segment(aes(x = 0.01, xend = 9.95, y = true_l_0, yend = true_l_0), color = "black",
               linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, trans = "log", 
                       breaks = scales::log_breaks(n = 5),
                       limits = range(means$mean_duration_rmse)) +  
  theme_classic() +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, seq(2, 10, 2))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "bottom", axis.title.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Duration of Clearance Phase", fill = "", x = "std dev of prior") +
  guides(fill = guide_colorbar(barwidth = 10))

p1 | p2 | p3

ggsave("fig/rmse_epi_vary_colorscale.png", 
       plot = (p1 | p2 | p3) & theme(plot.margin = margin(5, 15, 5, 15)), 
       width = 13, height = 4, dpi = 300)
