library(tidyverse)
source("src/immunity_semi_parametric_model.R")

context <- "sensitivity_analysis_cross_immunity_transmissibility"
source(path("scripts", "simulation", context, "shared_constants.txt"))

file_path_tbl <- 
  tibble(file_path = dir_ls(path("data", "simulation", context),
                            regexp = "[.]csv$")) %>% 
  mutate(file_name = file_path %>% 
           path_file() %>% 
           path_ext_remove()) %>% 
  mutate(file_type = str_remove_all(file_name, "_takeover_speed=.+"),
         takeover_speed = str_extract(file_name, "(?<=_takeover_speed=).+"))

all_dat <- 
  file_path_tbl %>% 
  filter(file_type == "simulated_data") %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_predictive())) %>% 
  select(-starts_with("file")) %>% 
  unnest(dat) %>% 
  mutate(variant_2_import_time = case_when(
    takeover_speed == "slow" ~ variant_2_import_time_slow,
    takeover_speed == "medium" ~ variant_2_import_time_medium,
    takeover_speed == "fast" ~ variant_2_import_time_fast,
  )) %>% 
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - 1 + (time - 1) / 7, time)) %>%
  mutate(time = time - (first_obs_time - 1))

all_generated_quantities <- 
  file_path_tbl %>% 
  filter(file_type == "true_generated_quantities") %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_generated_quantities())) %>% 
  select(-starts_with("file")) %>% 
  unnest(dat) %>% 
  mutate(variant_2_import_time = case_when(
    takeover_speed == "slow" ~ variant_2_import_time_slow,
    takeover_speed == "medium" ~ variant_2_import_time_medium,
    takeover_speed == "fast" ~ variant_2_import_time_fast,
  )) %>% 
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - first_obs_time + (time - 1) / 7, time))

popsize_tbl <- 
  file_path_tbl %>% 
  filter(file_type == "true_generated_quantities") %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_format_draws_time())) %>% 
  select(-starts_with("file")) %>% 
  unnest(dat) %>% 
  filter(!is.na(time)) %>% 
  distinct(takeover_speed, name, value, time) %>% 
  mutate(time = if_else(str_ends(name, "_mean"), time, time - 1)) %>% 
  filter(str_ends(name, "_compartment")) %>% 
  filter(!str_detect(name, "₀|₁|₂")) %>% 
  filter(name != "C_compartment") %>% 
  group_by(takeover_speed, time) %>% 
  summarize(popsize = sum(value), .groups = "drop")

max_t_tbl <- 
  all_generated_quantities %>% 
  filter(name == "seq_prop_variant_2") %>% 
  distinct(takeover_speed, time, value) %>% 
  mutate(gt_1 = value > 0.01,
         gt_5 = value > 0.05,
         gt_99 = value > 0.99) %>% 
  filter(gt_1 | gt_5 | gt_99) %>% 
  group_by(takeover_speed, gt_1, gt_5, gt_99) %>% 
  summarise(time = min(time), .groups = "drop") %>% 
  pivot_longer(starts_with("gt_"), names_prefix = "gt_", names_to = "percent_cutoff", names_transform = ~as.numeric(.x) / 100) %>% 
  filter(value == T) %>% 
  group_by(takeover_speed, percent_cutoff) %>% 
  summarize(time = min(time), .groups = "drop") %>% 
  mutate(name = fct_recode(as.character(percent_cutoff), first_max_t = "0.01", last_max_t = "0.99", offset_time = "0.05")) %>% 
  select(-percent_cutoff) %>% 
  pivot_wider(names_from = name, values_from = time) %>% 
  mutate(weeks_to_takeover = last_max_t - first_max_t) %>% 
  mutate(across(ends_with("max_t"), floor)) %>% 
  mutate(total_models = last_max_t - first_max_t + 1)


# Plots -------------------------------------------------------------------

all_dat %>% 
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = takeover_speed, color = takeover_speed, group = .width)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  geom_vline(data = max_t_tbl %>% 
               select(takeover_speed, ends_with("max_t")) %>% 
               pivot_longer(-takeover_speed) %>% 
               select(takeover_speed, time = value) %>% 
               expand_grid(all_dat %>% distinct(name)),
             mapping = aes(xintercept = time, color = takeover_speed),
             alpha = 0.5,
             linetype = "dashed")

all_generated_quantities %>% 
  filter(!is.na(time)) %>% 
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = takeover_speed, color = takeover_speed, group = .width)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5)

all_generated_quantities %>% 
  filter(is.na(time)) %>% 
  distinct(takeover_speed, name, value) %>% 
  ggplot(aes(value, takeover_speed, color = takeover_speed)) +
  facet_wrap(~name, scales = "free_x") +
  geom_point()

popsize_tbl %>% 
  ggplot(aes(time, popsize, color = takeover_speed)) +
  geom_line() +
  scale_y_continuous(labels = comma)

all_generated_quantities %>% 
  filter(name == "seq_prop_variant_2") %>% 
  distinct(takeover_speed, time, value) %>% 
  ggplot(aes(time, value, color = takeover_speed)) + 
  geom_line() +
  geom_vline(data = max_t_tbl %>% 
               select(takeover_speed, ends_with("max_t")) %>% 
               pivot_longer(-takeover_speed) %>% 
               select(takeover_speed, time = value) %>% 
               expand_grid(all_dat %>% distinct(name)),
             mapping = aes(xintercept = time, color = takeover_speed),
             alpha = 0.5,
             linetype = "dashed")
