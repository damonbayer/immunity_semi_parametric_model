library(tidyverse)
library(glue)
library(fs)
library(furrr)
source("src/immunity_semi_parametric_model.R")

options("parallelly.fork.enable" = TRUE)
plan(multicore)

experiment_name <-  "compare_inference_unknown_variant_proportion"
context <- path("simulation", experiment_name)
all_models_dir <- path("results", context)
sim_id <- 1

# Loading Data ------------------------------------------------------------
dat_tidy <-
  path("data", context, "simulated_data", ext = "csv") %>%
  read_csv_as_draws() %>%
  filter(.iteration == sim_id) %>%
  tidy_format_draws_time() %>%
  select(name, time, value) %>%
  mutate(time = time - min(time) + 1) %>%
  bind_rows(.,
            filter(., str_starts(name, "data_new_cases")) %>% 
              pivot_wider(names_from = name, values_from = value) %>% 
              mutate(name = "data_new_cases",
                     value = data_new_cases_variant_1 + data_new_cases_variant_2) %>% 
              select(time, name, value)
  ) %>% 
  mutate(name = str_remove(name, "data+_"))

results_tbl <- 
  dir_ls(all_models_dir, recurse = 2, type = "file") %>% 
  enframe(name = NULL, value = "file_path") %>% 
  filter(path_ext(file_path) == "csv") %>% 
  mutate(split_path = path_split(file_path)) %>% 
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>% 
  mutate(model_name = map_chr(split_path, ~pluck(.x, 4))) %>% 
  mutate(result_type = map_chr(split_path, ~pluck(.x, 5))) %>% 
  filter(str_starts(result_type, "tidy_") | str_detect(result_type, "score")) %>%
  mutate(distribution = str_extract(result_type, "prior|posterior")) %>% 
  mutate(result_type = str_remove(result_type, ".*(posterior|prior)_")) %>% 
  mutate(max_t = file_name %>% 
           str_extract("(?<=max_t=)\\d+") %>% 
           as.integer()) %>% 
  select(model_name, result_type, distribution, max_t, file_path)

tidy_predictive_tbl <-
  results_tbl %>% 
  filter(result_type == "predictive") %>% 
  mutate(tidy_predictive = future_map(file_path, read_csv)) %>% 
  unnest(tidy_predictive) %>% 
  select(-result_type, -file_path) %>% 
  mutate(name = str_remove(name, "data+_"))

tidy_generated_quantities_tbl <- 
  results_tbl %>% 
  filter(result_type == "generated_quantities") %>% 
  mutate(tidy_generated_quantities = future_map(file_path, read_csv)) %>% 
  unnest(tidy_generated_quantities) %>% 
  select(-result_type, -file_path) %>% 
  mutate(name = str_remove(name, "data+_"))

tidy_posterior_predictive_score_tbl <- 
  results_tbl %>% 
  filter(result_type == "predictive_score") %>% 
  select(-max_t) %>% 
  mutate(tidy_predictive_score = future_map(file_path, read_csv)) %>% 
  unnest(tidy_predictive_score) %>% 
  select(-c(model_name, result_type, file_path, time, max_t)) %>% 
  pivot_longer(cols = -c(distribution, model, target_type, weeks_ahead)) %>% 
  mutate(target_type = str_remove(target_type, "data+_")) %>% 
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Horizon"))

all_target_types <- unique(tidy_predictive_tbl$name)
all_model_names <- unique(tidy_generated_quantities_tbl$model_name)


# Plot functions ----------------------------------------------------------

plot_forecast_comparison <- function(target_type) {
  tmp_tidy_posterior_predictive <- 
    tidy_predictive_tbl %>% 
    filter(distribution == "posterior",
           name == target_type,
           weeks_ahead %in% c(0, 4, 8))
  
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
  tidy_posterior_predictive_score_tbl %>% 
    filter(target_type == target_target_type,
           weeks_ahead %in% c(1, 4, 8)) %>% 
    ggplot(aes(model, value, color = model)) +
    facet_grid(name ~ forecast_horizon, scales = "free_y") +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_y_continuous(name = "Value", labels = comma) +
    scale_x_discrete(name = "Model") +
    theme(legend.position = "bottom",
          axis.text.x = element_blank()) +
    scale_color_discrete(name = "Model") +
    ggtitle(glue("Forecast Metrics Comparison for {str_to_title(str_replace_all(target_target_type, '_', ' '))}"))
}

plot_scalar_generated_quantities_by_forecast_time <- function(target_model_name) {
  tidy_generated_quantities_tbl %>% 
    filter(is.na(time), .width == 0.8) %>% 
    filter(model_name == target_model_name) %>% 
    ggplot(aes(max_t, value, ymin = .lower, ymax = .upper, color = distribution)) +
    facet_wrap(~name, scales = "free_y") +
    geom_interval(alpha = 0.5) +
    scale_y_continuous("Value") +
    scale_x_continuous("Forecast Time") +
    ggtitle(glue("Posterior Parameter Distirbutions by Forecast time for {target_model_name}"))
}

plot_vector_generated_quantities_by_forecast_time <- function(target_model_name) {
  tidy_generated_quantities_tbl %>% 
    filter(!is.na(time),
           distribution == "posterior") %>% 
    filter(model_name == target_model_name) %>% 
    distinct(name, time, value, max_t) %>% 
    ggplot(aes(time, value, color = max_t, group = max_t)) +
    facet_wrap(~name, scales = "free_y") +
    scale_color_viridis_c("Forecast Time") +
    geom_line(alpha = 1) +
    scale_x_continuous("Time") +
    scale_y_continuous("Value") +
    ggtitle(glue("Posterior Median by Forecast time for {target_model_name}")) +
    theme(legend.position = "bottom")
}

plot_single_posterior_predictive <- function(target_model_name, target_max_t) {
  tmp_tidy_posterior_predictive <- 
    tidy_predictive_tbl %>% 
    filter(distribution == "posterior",
           model_name == target_model_name,
           max_t == target_max_t)
  
  tmp_dat_tidy <-
    dat_tidy %>% 
    filter(time <= max(tmp_tidy_posterior_predictive$time)) %>% 
    filter(name %in% unique(tmp_tidy_posterior_predictive$name)) %>% 
    mutate(data_type = if_else(time <= target_max_t, "observed", "future"))
  
  ggplot(mapping = aes(time, value)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(data = tmp_tidy_posterior_predictive,
                    mapping = aes(ymin = .lower, ymax = .upper), step = "mid") +
    geom_point(data = tmp_dat_tidy, mapping = aes(shape = data_type)) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    my_theme +
    scale_shape_discrete(name = "Data Type", labels = str_to_title) +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("{target_model_name} Posterior Predictive at t = {target_max_t}"))
}

plot_single_generated_quantities <- function(target_model_name, target_max_t) {
  tidy_generated_quantities_tbl %>% 
    filter(!is.na(time),
           model_name == target_model_name,
           max_t == target_max_t) %>% 
    ggplot(mapping = aes(time, value, ymin = .lower, ymax = .upper, color = distribution, fill = distribution, group = .width)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(alpha = 0.25) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    # my_theme +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("{target_model_name} Generated Quantities at t = {target_max_t}"))
}


# Create Figures ----------------------------------------------------------
augment_figure_tbl <- function(figure_tbl) {
  figure_tbl %>% 
    mutate(figure_dims = future_map(figure, gg_facet_dims)) %>% 
    unnest_wider(figure_dims) %>% 
    rename(n_row = ROW, n_col = COL) %>% 
    select(file_path, figure, n_col, n_row, everything())
}


forecast_comparison_plots <-
  tibble(target_type = all_target_types) %>% 
  mutate(file_path = path("figures", experiment_name, glue("forecast_comparison_{target_type}_plot"), ext = "pdf"),
         figure = future_map(target_type, plot_forecast_comparison)) %>% 
  augment_figure_tbl()

forecast_metrics_comparison_plots <- 
  tibble(target_type = all_target_types) %>% 
  mutate(file_path = path("figures", experiment_name, glue("forecast_metrics_comparison_{target_type}_plot"), ext = "pdf"),
         figure = future_map(target_type, plot_forecast_metrics_comparison)) %>% 
  augment_figure_tbl()

scalar_generated_quantities_by_forecast_time_plots <- 
  tibble(target_model_name = all_model_names) %>% 
  mutate(file_path = path("figures", experiment_name, glue("scalar_generated_quantities_by_forecast_time_{target_model_name}_plot"), ext = "pdf"),
         figure = future_map(target_model_name, plot_scalar_generated_quantities_by_forecast_time)) %>% 
  augment_figure_tbl()

vector_generated_quantities_by_forecast_time_plots <- 
  tibble(target_model_name = all_model_names) %>% 
  mutate(file_path = path("figures", experiment_name, glue("vector_generated_quantities_by_forecast_time_{target_model_name}_plot"), ext = "pdf"),
         figure = future_map(target_model_name, plot_vector_generated_quantities_by_forecast_time)) %>% 
  mutate(figure_dims = future_map(figure, gg_facet_dims)) %>% 
  augment_figure_tbl()

single_generated_quantities_plots <- 
  tidy_predictive_tbl %>% 
  distinct(model_name, max_t) %>% 
  rename_with(~str_c("target_", .x)) %>% 
  mutate(file_path = path("figures", experiment_name, glue("single_generated_quantities_{target_model_name}_{target_max_t}_plot"), ext = "pdf"),
         figure = map2(target_model_name, target_max_t, ~plot_single_generated_quantities(target_model_name = .x, target_max_t = .y))) %>% 
  augment_figure_tbl()

single_posterior_predictive_plots <- 
  tidy_predictive_tbl %>% 
  distinct(model_name, max_t) %>% 
  rename_with(~str_c("target_", .x)) %>% 
  mutate(file_path = path("figures", experiment_name, glue("single_posterior_predictive_{target_model_name}_{target_max_t}_plot"), ext = "pdf"),
         figure = map2(target_model_name, target_max_t, ~plot_single_posterior_predictive(target_model_name = .x, target_max_t = .y))) %>% 
  augment_figure_tbl()


# Save Figures ------------------------------------------------------------
dir_create(path("figures", experiment_name))
dir_create(path("figures", experiment_name, all_model_names))

all_plot_tbl_names <- ls()[str_ends(ls(), "_plots")]
single_plot_tbl_names <- all_plot_tbl_names[str_starts(all_plot_tbl_names, "single", negate = F)]
non_single_plot_tbl_names <- all_plot_tbl_names[str_starts(all_plot_tbl_names, "single", negate = T)]

# Save all figures
future_walk(all_plot_tbl_names, ~pwalk(as.list(get(.x)), ~save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, device = cairo_pdf)))

# Merge all figures
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))

non_single_plot_tbl_names %>% 
  walk(~{
    uncollected_plot_paths <- get(.x)$file_path
    collected_plot_path <- path("figures", experiment_name, str_c("all_", .x), ext = "pdf")
    
    str_c(uncollected_plot_paths, collapse = " ") %>% 
      str_c(collected_plot_path, sep = " ") %>% 
      system2("pdfunite", args = .)
    
    file_delete(uncollected_plot_paths)
  })

merge_single_plots <- function(single_plot_tbl_name) {
  single_plot_tbl <- get(single_plot_tbl_name)
  single_plot_tbl %>% 
    group_by(target_model_name) %>% 
    group_walk(~{
      uncollected_plot_paths <- .x$file_path
      collected_plot_path <- path("figures", experiment_name, .y, str_c("all_", single_plot_tbl_name, "_", .y), ext = "pdf")
      
      str_c(uncollected_plot_paths, collapse = " ") %>%
        str_c(collected_plot_path, sep = " ") %>%
        system2("pdfunite", args = .)
    })
  
  file_delete(single_plot_tbl$file_path)
}

get(single_plot_tbl_names[1])$file_path

walk(single_plot_tbl_names, merge_single_plots)
