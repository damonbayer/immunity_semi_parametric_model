source("src/immunity_semi_parametric_model.R")

library(fs)


all_dat <- 
  tibble(file_path =dir_ls(path("data", "simulation"), recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "simulated_data.csv") %>% 
  mutate(simulation = map_chr(path_split(file_path), ~pluck(.x, 3))) %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_predictive())) %>% 
  select(-file_path) %>% 
  unnest(dat)

all_generated_quantities <- 
  tibble(file_path =dir_ls(path("data", "simulation"), recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "true_generated_quantities.csv") %>% 
  mutate(simulation = map_chr(path_split(file_path), ~pluck(.x, 3))) %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_generated_quantities())) %>% 
  select(-file_path) %>% 
  unnest(dat)

all_dat %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = simulation, group = .width)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5)

all_generated_quantities %>% 
  filter(!is.na(time)) %>% 
  distinct(simulation, time, name, value) %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, color = simulation)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line(alpha = 0.5)

target_simulation_name <- "compare_inference_unknown_variant_proportion_ihr"

# First day over 1%
first_time_new_var <- 
  all_generated_quantities %>% 
  filter(simulation == target_simulation_name,
         name == "prop_variant_2") %>% 
  distinct(time, name, value) %>% 
  filter(value >= 0.01) %>% 
  pull(time) %>% 
  pluck(1)

time_5_percent <-
  all_generated_quantities %>% 
  filter(simulation == target_simulation_name,
         name == "prop_variant_2") %>% 
  distinct(time, name, value) %>%
  arrange(abs(value - 0.05)) %>% 
  pull(time) %>% 
  pluck(1)

# First day over 99%
first_time_sat <- 
  all_generated_quantities %>% 
  filter(simulation == target_simulation_name,
         name == "prop_variant_2") %>% 
  distinct(time, name, value) %>% 
  filter(value >= 0.99) %>% 
  pull(time) %>% 
  pluck(1)


first_time_second_hosp_wave <- 
  all_generated_quantities %>% 
  filter(simulation == target_simulation_name,
         name == "H") %>% 
  distinct(time, name, value) %>% 
  filter(sign(value - lag(value)) == -1) %>% 
  filter(time - lag(time) > 1) %>% 
  pull(time) %>% 
  pluck(1)

time_second_hosp_peak <-
  all_generated_quantities %>% 
  filter(simulation == target_simulation_name,
         name == "H",
         time >= first_time_second_hosp_wave) %>% 
  distinct(time, name, value) %>% 
  filter(value == max(value)) %>% 
  pull(time) %>% 
  pluck(1)


# Final Date Range
c(first_time_new_var, max(first_time_sat, time_second_hosp_peak))
