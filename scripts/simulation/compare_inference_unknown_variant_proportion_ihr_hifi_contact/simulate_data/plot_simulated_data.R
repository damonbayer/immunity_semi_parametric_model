library(tidyverse)
source("src/immunity_semi_parametric_model.R")

variant_2_import_time <- 40
first_obs_time <- 20

dat <- read_csv_as_draws("data/simulation/compare_inference_unknown_variant_proportion_ihr_hifi_contact/simulated_data.csv") %>%
  tidy_predictive() %>%
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - 1 + (time - 1) / 7, time)) %>%
  mutate(time = time - (first_obs_time - 1))

dat %>% group_by(name) %>% summarize(min(time), max(time))

gq <- read_csv_as_draws("data/simulation/compare_inference_unknown_variant_proportion_ihr_hifi_contact/true_generated_quantities.csv") %>%
  tidy_generated_quantities() %>%
  distinct(name, time, value) %>%
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time + (time - 1) / 7, time + 20)) %>%
  mutate(time = time - (first_obs_time - 1))

gq %>%
  filter(!is.na(time)) %>%
  filter(str_detect(name, "prop_variant_2")) %>%
  group_by(name) %>%
  summarize(min(time), max(time))

gq %>%
  filter(!is.na(time)) %>%
  filter(str_detect(name, "prop_variant_2")) %>%
  ggplot(aes(time, value, color = name)) +
  geom_point(alpha = 0.5)


dat %>%
  ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon() +
  my_theme


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
  filter(!str_detect(name, "prop_variant_2")) %>%
  ggplot(aes(time, value)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line()

# revise this
# also find closest to 5%
gq %>%
  filter(str_detect(name, "seq_prop_variant_2")) %>%
  group_by(gt1 = value > 0.01,
           gt99 = value > 0.99) %>%
  summarize(min(time))


gq %>%
  filter(str_detect(name, "seq_prop_variant_2")) %>%
  arrange(abs(value - 0.05))
