library(tidyverse)
library(glue)
library(fs)
source("src/immunity_semi_parametric_model.R")
context <- path("simulation", "compare_inference_known_variant_proportion")
all_models_dir <- path("results", context)
sim_id <- 1

dat_tidy <-
  path("data", context, "simulated_data", ext = "csv") %>%
  read_csv_as_draws() %>%
  filter(.iteration == sim_id) %>%
  tidy_format_draws_time() %>%
  select(name, time, value) %>%
  mutate(time = time - min(time) + 1) %>% 
  mutate(name = str_remove(name, "data+_"))

tidy_posterior_predictive <-
  dir_ls(all_models_dir, recurse = 2) %>% 
  enframe(name = NULL, value = "file_path") %>% 
  mutate(file_name = file_path %>% 
           path_file() %>% 
           path_ext_remove()) %>% 
  filter(str_detect(file_name, "^tidy_posterior_predictive_max_t=\\d+$")) %>% 
  mutate(model_name = map_chr(path_split(file_path), ~pluck(.x, 4)), .before = 1) %>% 
  mutate(max_t = file_name %>% 
           str_extract("(?<=^tidy_posterior_predictive_max_t=)\\d+") %>% 
           as.integer()) %>% 
  mutate(tidy_posterior_predictive = map(file_path, read_csv)) %>% 
  unnest(tidy_posterior_predictive) %>% 
  select(-c(file_path, file_name)) %>% 
  mutate(name = str_remove(name, "data+_"))

tidy_posterior_predictive_score <- 
  dir_ls(all_models_dir, recurse = 2) %>% 
  enframe(name = NULL, value = "file_path") %>% 
  filter(file_path %>% 
           path_file() %>% 
           path_ext_remove() %>% 
           str_detect("^posterior_predictive_score_max_t=\\d+$")) %>% 
  mutate(posterior_predictive_score = map(file_path, read_csv)) %>% 
  select(-file_path) %>% 
  unnest(posterior_predictive_score) %>% 
  select(-time, -max_t) %>% 
  pivot_longer(cols = -c(model, target_type, weeks_ahead)) %>% 
  mutate(target_type = str_remove(target_type, "data+_")) %>% 
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Horizon"))

all_target_types <- unique(tidy_posterior_predictive_score$target_type)

plot_forecast_comparison <- function(target_type) {
  tmp_tidy_posterior_predictive <- 
    tidy_posterior_predictive %>% 
    filter(name == target_type,
           weeks_ahead %in% c(0, 1, 4))
  
  tmp_dat_tidy <-
    dat_tidy %>% 
    filter(name == target_type)
  
  ggplot(mapping = aes(time, value)) +
    facet_grid(weeks_ahead ~ model_name, scale = "free_y",
               labeller = labeller(
                 weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
                 model_name = as_labeller(~str_replace_all(.x, "_", " ")))) +
    geom_lineribbon(data = tmp_tidy_posterior_predictive,
                    mapping = aes(ymin =.lower, ymax = .upper), step = "mid") +
    geom_point(data = tmp_tidy_posterior_predictive %>%
                 group_by(weeks_ahead) %>%
                 summarize(min_time = min(time),
                           max_time = max(time)) %>%
                 expand_grid(tmp_dat_tidy) %>%
                 filter(time <= max_time)) +
    scale_y_continuous(name = "Value", labels = comma) +
    scale_x_continuous(name = "Time") +
    my_theme +
    ggtitle(glue("Forecast Comparison for {str_to_title(str_replace_all(target_type, '_', ' '))}"))
}

plot_forecast_metrics_comparison <- function(target_target_type) {
  tidy_posterior_predictive_score %>%
    filter(model != "SEIRS",
           target_type == target_target_type) %>% 
    ggplot(aes(model, value, color = model)) +
    facet_grid(name ~ forecast_horizon, scales = "free_y") +
    geom_boxplot() +
    scale_y_continuous(name = "Value", labels = comma) +
    scale_x_discrete(name = "Model") +
    theme(legend.position = "bottom",
          axis.text.x = element_blank()) +
    scale_color_discrete(name = "Model") +
    ggtitle(glue("Forecast Metrics Comparison for {str_to_title(str_replace_all(target_target_type, '_', ' '))}"))
}

forecast_comparison_plots <-
  tibble(target_type = all_target_types,
         file_path = path("figures", glue("forecast_comparison_{all_target_types}_plot"), ext = "pdf"),
         figure = map(all_target_types, plot_forecast_comparison))

forecast_metrics_comparison_plots <-
  tibble(target_type = all_target_types,
         file_path = path("figures", glue("forecast_metrics_comparison_{all_target_types}_plot"), ext = "pdf"),
         figure = map(all_target_types, plot_forecast_metrics_comparison))


walk2(forecast_comparison_plots$file_path, forecast_comparison_plots$figure,
      ~save_plot(filename = .x,
                plot = .y,
                ncol = 4,
                nrow = 3))

walk2(forecast_metrics_comparison_plots$file_path, forecast_metrics_comparison_plots$figure,
      ~save_plot(filename = .x,
                 plot = .y,
                 ncol = 4,
                 nrow = 7,
                 base_height = 2))


Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))

str_c(forecast_comparison_plots$file_path, collapse = " ") %>% 
  str_c(path("figures", "all_forecast_comparison_plots", ext = "pdf"), sep = " ") %>% 
  system2("pdfunite", args = .)

str_c(forecast_metrics_comparison_plots$file_path, collapse = " ") %>% 
  str_c(path("figures", "all_forecast_metrics_comparison_plots", ext = "pdf"), sep = " ") %>% 
  system2("pdfunite", args = .)