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

simulated_data <- 
  file_path_tbl %>% 
  filter(file_type == "simulated_data") %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_format_draws_time())) %>% 
  select(takeover_speed, dat) %>% 
  unnest(dat) %>% 
  filter(.draw == 1) %>% 
  select(takeover_speed, name, time, value) %>% 
  mutate(variant_2_import_time = case_when(
    takeover_speed == "slow" ~ variant_2_import_time_slow,
    takeover_speed == "medium" ~ variant_2_import_time_medium,
    takeover_speed == "fast" ~ variant_2_import_time_fast,
  )) %>% 
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - 1 + (time - 1) / 7, time)) %>%
  mutate(time = time - (first_obs_time - 1)) %>% 
  mutate(name = fct_relevel(name, data_order)) %>% 
  bind_rows(.,
            filter(., str_detect(name, "seq")) %>% 
              group_by(takeover_speed, time, variant_2_import_time) %>% 
              summarize(value = sum(value), .groups = "drop") %>% 
              mutate(name = "data_new_seq_total")) %>% 
  filter(name != "data_new_seq_variant_1")

model_table <- read_csv(path("scripts", "simulation", context, "model_table.csv"))

model_details <- 
  model_table %>% 
  group_by(data_takeover_speed) %>% 
  summarize(min_max_t = min(max_t),
            max_max_t = max(max_t)) %>% 
  mutate(first_forecast_target_t = min_max_t + 1,
         last_forecast_target_t = max_max_t + 4) %>% 
  rename(takeover_speed = data_takeover_speed)

simulated_binned_data_plots_tbl <- 
  simulated_data %>% 
  left_join(model_details) %>% 
  filter(time <= last_forecast_target_t) %>% 
  group_by(takeover_speed, first_forecast_target_t, last_forecast_target_t, variant_2_import_time) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(figure = pmap(list(takeover_speed, first_forecast_target_t, last_forecast_target_t, data), ~{
    ggplot(..4, aes(time, value)) +
      facet_wrap(~name, scales = "free_y", labeller = my_sim_labeller_fn) +
      annotate("rect", xmin = ..2, xmax = ..3, ymin = -Inf, ymax = Inf, alpha = 0.5) +
      geom_point() +
      geom_line() +
      scale_x_continuous("Time") +
      scale_y_continuous("Count", label = comma) +
      ggtitle("Simulated Data",
              subtitle = glue("{str_to_title(..1)} Novel Variant Takeover Speed"))
  })) %>% 
  mutate(figure_dims = map(figure, gg_facet_dims)) %>%
  unnest_wider(figure_dims) %>%
  rename(n_row = ROW, n_col = COL) %>% 
  mutate(file_path = path(manuscript_figure_dir, glue("simulated_binned_data_{takeover_speed}_plot"), ext = "pdf"))

simulated_binned_data_plots_tbl %>% 
  select(file_path, figure, n_col, n_row) %>% 
  as.list() %>% 
  pwalk(~save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 1, base_height = 3))

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

proportion_novel_variant_simulated_data_plot <- 
  all_generated_quantities %>% 
  filter(name %in% c("prop_variant_2", "seq_prop_variant_2")) %>% 
  select(takeover_speed, time, value) %>% 
  filter(time <= max(model_details$last_forecast_target_t)) %>% 
  ggplot(aes(time, value, color = takeover_speed, group = takeover_speed)) +
  geom_vline(data = all_generated_quantities %>%
               distinct(takeover_speed, variant_2_import_time) %>% 
               mutate(variant_2_import_time = variant_2_import_time - first_obs_time + 1),
             mapping = aes(xintercept = variant_2_import_time, color = takeover_speed),
             linetype = "dashed",
             alpha = 0.5,
             show.legend = F) +
  geom_line() +
  scale_y_continuous(TeX(r"(Novel Variant Proportion $\left( \delta(t) \right)$)"),
                     labels = percent) +
  scale_x_continuous("Time") +
  scale_color_discrete("Novel Variant Takeover Speed", labels = str_to_title) +
  theme(legend.position = "bottom") +
  ggtitle("Novel Variant Proportion for Simulated Datasets", subtitle = "Dashed line indicates first novel variant import time")

save_plot(filename = path(manuscript_figure_dir, "proportion_novel_variant_simulated_data_plot", ext = "pdf"),
          plot = proportion_novel_variant_simulated_data_plot,
          base_height = 4,
          base_asp = 1.75)
