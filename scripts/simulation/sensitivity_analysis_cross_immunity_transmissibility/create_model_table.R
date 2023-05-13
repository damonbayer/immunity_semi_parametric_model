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


computed_shared_constants <- 
  max_t_tbl %>% 
  select(takeover_speed,
         offset_time,
         weeks_to_takeover) %>% 
  mutate(logistic_growth_time_offset = round(offset_time),
         time_look_for_second_wave = ceiling(offset_time),
         time_to_saturation = round(weeks_to_takeover)) %>% 
  select(-offset_time, - weeks_to_takeover) %>% 
  pivot_longer(-takeover_speed) %>% 
  mutate(text = glue("{name}_{takeover_speed} = {formatC(value, format = 'f', digits = 1)}")) %>% 
  arrange(name) %>% 
  pull(text)

model_table <- 
  bind_rows(
  tibble(immunity_model = "constant",
         `R₀_model` = "gmrf",
         prior_takeover_speed = NA),
  tibble(immunity_model = "gmrf",
         `R₀_model` = "constant",
         prior_takeover_speed = NA),
expand_grid(immunity_model = "seq-informed",
            `R₀_model` = "constant",
            prior_takeover_speed = c("slow", "medium", "fast", "wide"))) %>% 
  mutate(CDR_model = "constant") %>% 
expand_grid(max_t_tbl %>% 
              mutate(max_t = map2(first_max_t, last_max_t, ~.x:.y)) %>% 
              select(data_takeover_speed = takeover_speed, max_t) %>% 
              unnest(max_t)) %>% 
  mutate(sim_id = 1) %>% 
  mutate(fit_id = 1:n()) %>% 
  select(fit_id, sim_id, max_t, data_takeover_speed, everything())

write_lines(computed_shared_constants, path("scripts", "simulation", context, "computed_shared_constants", ext = "txt"))
write_csv(model_table, path("scripts", "simulation", context, "model_table", ext = "csv"))
