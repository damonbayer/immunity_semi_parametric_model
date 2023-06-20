library(tidyverse)
source("src/immunity_semi_parametric_model.R")
library(latex2exp)
library(fs)

alpha_2_choices <- c(0.1, 0.3, 1, 3, 10)
max_dur_choice <- 100
alpha_1_star_choice <- 0.5

annotation_data <-
  tibble(
    x = c(0.0, 0.5),
    y = c(max_dur_choice * alpha_1_star_choice, max_dur_choice),
    label = TeX(c(
      r'($\alpha_0 \cdot \alpha_1$)',
      r'($\alpha_0$)'
    ))
  )

example_immune_durration_plot <-
  expand_grid(
    max_dur = max_dur_choice,
    alpha_2 = alpha_2_choices,
    alpha_1_star = alpha_1_star_choice
  ) %>%
  mutate(
    alpha_0 = log(max_dur),
    alpha_1 = 4^alpha_2 * log(alpha_1_star)
  ) %>%
  expand_grid(x = seq(0, 1, length.out = 1001)) %>%
  mutate(dur = exp(alpha_0 + alpha_1 * (x * (1 - x))^alpha_2)) %>%
  ggplot(aes(x, dur, color = factor(alpha_2), group = alpha_2)) +
  geom_hline(
    data = annotation_data,
    mapping = aes(yintercept = y),
    linetype = "dashed",
    color = "gray30"
  ) +
  geom_line() +
  geom_text(
    data = annotation_data,
    mapping = aes(x, y, label = label),
    parse = TRUE,
    inherit.aes = F, vjust = "inward", hjust = "inward"
  ) +
  theme_cowplot() +
  scale_x_continuous(TeX(r"($\delta$)"), labels = percent) +
  scale_y_continuous(TeX(r"($1 / \kappa(\delta)$)")) +
  scale_color_discrete(TeX(r'($\alpha_2$)')) +
  ggtitle(
    # label = TeX(r"(Example $1 / \kappa(\delta) = \exp \left{\alpha_0 + \alpha_1 \left[ \delta \left( 1 - \delta \right)^{\alpha_2} \right]\right}$)"),
    label = "Example Immunity Duration Curves",
    subtitle = TeX(sprintf(r'($\alpha_0 = %d$, $\alpha_1 = %.1f)', max_dur_choice, alpha_1_star_choice))
  )


iota_0 <- qlogis(0.05)
t_star <- 10

example_delta_t_plot <- 
  tibble(t_star = t_star,
         iota_0 = iota_0) %>% 
  expand_grid(t = seq(0, 30, by = 0.1),
              iota_1_star = 5 * 1:4) %>% 
  mutate(iota_1 = (qlogis(0.99) - qlogis(0.01)) / iota_1_star) %>% 
  mutate(delta_t = plogis(iota_0 + iota_1 * (t - t_star))) %>% 
  ggplot(aes(t, delta_t, color = factor(iota_1_star))) +
  geom_hline(yintercept = plogis(iota_0), linetype = "dashed", color = "gray50"
  ) +
  annotate(geom = "text", x = 0, y = plogis(iota_0), label = TeX(r"($logit(\iota_0)$)"), hjust = "left", vjust = "bottom") +
  geom_vline(xintercept = t_star, linetype = "dashed", color = "gray50"
  ) +
  annotate(geom = "text", x = t_star, y = 0.75, label = TeX(r"($t^*$)"), hjust = "right", vjust = "top") +
  geom_line() +
  scale_y_continuous(TeX(r"($\delta(t)$)"), labels = percent) +
  scale_color_discrete(TeX(r"($\iota_1^*$)"))

# CRPS
x <- 0.25

crps_data <- tibble(y = seq(-3, 3, length.out = 1000)) %>% 
  mutate(F_y = pnorm(y)) %>% 
  mutate(one = as.numeric(y >= x)) %>% 
  mutate(F_y_minus_one = F_y - one) %>% 
  mutate(F_y_minus_one_sqr = F_y_minus_one^2) %>% 
  mutate(x_min = 0 * (y < x) + F_y * (y >= x)) %>% 
  mutate(x_max = F_y * (y < x) + 1 * (y >= x))

crps_illustration_plot <- 
  crps_data %>%
  select(y, F_y, one) %>% 
  pivot_longer(-y) %>% 
  ggplot(aes(y, value, color = name)) +
  geom_line() +
  geom_ribbon(data = crps_data,
              mapping = aes(x = y, ymin = x_min, ymax = x_max), alpha = 0.1, inherit.aes = F) +
  cowplot::theme_cowplot() +
  scale_y_continuous(NULL) +
  theme(legend.position = "none") +
  ggtitle("CRPS Illustration")


# Save figures ------------------------------------------------------------

save_plot_target_asp(filename = path(defense_figure_dir, "example_delta_t_plot", ext = "pdf"), plot = example_delta_t_plot, base_asp = slide_target_asp / 2)


save_plot_target_asp(filename = path(defense_figure_dir, "example_immune_durration_plot", ext = "pdf"), plot = example_immune_durration_plot, base_asp = slide_target_asp / 2)

save_plot(filename = path(defense_figure_dir, "crps_illustration_plot", ext = "pdf"), plot = crps_illustration_plot, base_asp = slide_target_asp/2)
