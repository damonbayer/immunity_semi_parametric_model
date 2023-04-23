library(tidyverse)
source("src/immunity_semi_parametric_model.R")

dat <- read_csv_as_draws("data/simulation/compare_inference_unknown_variant_proportion_ihr/simulated_data.csv") %>% tidy_predictive()

gq <- read_csv_as_draws("data/simulation/compare_inference_unknown_variant_proportion_ihr/true_generated_quantities.csv") %>% tidy_generated_quantities()

dat %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon() +
  my_theme

gq %>% filter(name == "S")

gq %>%
  filter(!is.na(time)) %>%
  filter(str_detect(name, "₁|₂|₀", negate = F)) %>%
  mutate(compartment = str_extract(name, "^[^₀₁₂]*")) %>%
  mutate(suffix = str_extract(name, "[₀₁₂]+")) %>%
  ggplot(aes(time, value, color = suffix)) +
  facet_wrap(~compartment, scales = "free_y") +
  geom_line()

gq %>%
  filter(!is.na(time)) %>%
  filter(str_detect(name, "₁|₂|₀", negate = T)) %>%
  filter(name != "prop_variant_2") %>%
  ggplot(aes(time, value)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line()


gq %>%
  filter(name == "prop_variant_2") %>%
  ggplot(aes(time, value)) +
  geom_line() +
  geom_point()


gq %>%
  filter(name == "prop_variant_2") %>% 
  group_by(gt1 = value > 0.01,
           gt99 = value > 0.99) %>% 
  summarize(min(time))
