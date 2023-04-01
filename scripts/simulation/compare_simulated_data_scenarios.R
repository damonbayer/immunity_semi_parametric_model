source("src/immunity_semi_parametric_model.R")

dat_1 <- read_csv_as_draws("data/simulation/compare_inference_known_variant_proportion/simulated_data.csv") %>%
  tidy_predictive() %>% 
  mutate(simulation = "old")

dat_2 <- read_csv_as_draws("data/simulation/compare_inference_known_variant_proportion_2/simulated_data.csv") %>% tidy_predictive() %>% 
  mutate(simulation = "new")


bind_rows(dat_1, dat_2) %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = simulation, group = .width)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5)

gq_1 <- read_csv_as_draws("data/simulation/compare_inference_known_variant_proportion/true_generated_quantities.csv") %>%
  tidy_generated_quantities() %>%
  mutate(simulation = "old") %>%
  filter(name == "prop_variant_2")

gq_2 <- read_csv_as_draws("data/simulation/compare_inference_known_variant_proportion_2/true_generated_quantities.csv") %>%
  tidy_generated_quantities() %>% 
  mutate(simulation = "new") %>%
  filter(name == "prop_variant_2")


bind_rows(gq_1, gq_2) %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, color = simulation)) +
  geom_line() +
  scale_y_continuous("Proportion Variant 2", labels = percent)
