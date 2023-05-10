library(tidyverse)
library(glue)
library(fs)
library(furrr)
source("src/immunity_semi_parametric_model.R")
county_id <- 30
options("parallelly.fork.enable" = TRUE)
plan(multicore)

experiment_name <- "ba1_immunity_contact_ihr"
context <- path("real_data", experiment_name)
data_dir <- path("data/real_data/ba1")
all_models_dir <- path("results", context)

source(path(data_dir, "shared_constants.txt"))
date_time_0 <- ymd(date_time_0)
time_to_date <- function(time) date_time_0 + time  * 7
date_to_time <- function(date) as.numeric((date - date_time_0) / 7)

# Loading Data ------------------------------------------------------------
dat_tidy <- 
  bind_rows(
    path(data_dir, "cdph_data", ext = "csv") %>% 
      read_csv() %>% 
      filter(id == county_id) %>% 
      select(date = end_date,
             data_new_cases = cases,
             data_new_deaths = deaths,
             data_hospitalizations = hospitalized,
             data_icu = icu) %>% 
      pivot_longer(-date),
    path(data_dir, "seq_data", ext = "csv") %>% 
      read_csv() %>% 
      filter(id == county_id) %>% 
      select(date,
             data_new_seq_variant_2 = lineage_count,
             data_new_seq_variant_1 = other_count) %>% 
      pivot_longer(-date)) %>% 
  mutate(name = str_remove(name, "data+_"))


true_peak_dat <-
  dat_tidy %>%
  mutate(time = date_to_time(date)) %>% 
  select(-date) %>% 
  filter(time > time_look_for_second_wave,
         !str_detect(name, "variant")) %>%
  group_by(name) %>%
  filter(value== max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-name, names_to = "peak_type")

results_tbl <-
  dir_ls(all_models_dir, recurse = 2, type = "file") %>%
  enframe(name = NULL, value = "file_path") %>%
  filter(path_ext(file_path) == "csv") %>%
  mutate(split_path = path_split(file_path)) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  mutate(model_name = map_chr(split_path, ~pluck(.x, 4))) %>%
  mutate(result_type = map_chr(split_path, ~pluck(.x, 5))) %>%
  filter(str_starts(result_type, "tidy_") | str_detect(result_type, "score") | str_detect(result_type, "posterior_samples")) %>%
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
  mutate(name = str_remove(name, "data+_")) %>%
  mutate(name = fct_reorder(name, str_detect(name, "mean") * 1 + str_detect(name, "compartment") * 2))

tidy_posterior_predictive_score_tbl <-
  results_tbl %>%
  filter(result_type == "predictive_score") %>%
  select(-max_t) %>%
  mutate(tidy_predictive_score = future_map(file_path, read_csv)) %>%
  unnest(tidy_predictive_score) %>%
  select(-c(model_name, result_type, file_path, date, max_t)) %>%
  pivot_longer(cols = -c(distribution, model, target_type, weeks_ahead)) %>%
  mutate(target_type = str_remove(target_type, "data+_")) %>%
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Horizon"))

tidy_posterior_peak <-
  results_tbl %>%
  filter(result_type == "peak") %>%
  mutate(peak = future_map(file_path, read_csv)) %>%
  unnest(peak) %>%
  select(-result_type, -file_path) %>%
  mutate(name = str_remove(name, "data+_"))

all_target_types <- unique(tidy_predictive_tbl$name)
all_model_names <- unique(tidy_generated_quantities_tbl$model_name)
all_peak_types <- unique(tidy_posterior_peak$peak_type)

# Plot functions ----------------------------------------------------------

plot_forecast_comparison <- function(target_type) {
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    filter(distribution == "posterior",
           name == target_type,
           weeks_ahead %in% c(0, 1, 2, 4))

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(name == target_type)

  ggplot(mapping = aes(date, value)) +
    facet_grid(weeks_ahead ~ model_name, scale = "free_y",
               labeller = labeller(
                 weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
                 model_name = as_labeller(~str_replace_all(.x, "_", " ")))) +
    geom_lineribbon(data = tmp_tidy_posterior_predictive,
                    mapping = aes(ymin =.lower, ymax = .upper), step = "mid") +
    geom_point(data = tmp_tidy_posterior_predictive %>%
                 group_by(weeks_ahead) %>%
                 summarize(min_date = min(date),
                           max_date = max(date)) %>%
                 expand_grid(tmp_dat_tidy) %>%
                 filter(date <= max_date)) +
    scale_y_continuous(name = "Value", labels = comma) +
    scale_x_date(name = "Date") +
    my_theme +
    ggtitle(glue("Forecast Comparison for {str_to_title(str_replace_all(target_type, '_', ' '))}"))
}

plot_forecast_metrics_comparison <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(target_type == target_target_type,
           weeks_ahead %in% c(1, 2, 4)) %>%
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
  available_generated_quantities_names <-
    tidy_generated_quantities_tbl %>%
    filter(model_name == target_model_name) %>%
    filter(name != "prop_variant_2") %>%
    pull(name) %>%
    unique()

  tidy_generated_quantities_tbl %>%
    mutate(max_date = time_to_date(max_t)) %>% 
    filter(is.na(date), .width == 0.8) %>%
    filter(model_name == target_model_name) %>%
    ggplot(aes(max_date, value, ymin = .lower, ymax = .upper, color = distribution)) +
    facet_wrap(~name, scales = "free_y") +
    geom_interval(alpha = 0.5) +
    scale_y_continuous("Value") +
    scale_x_date("Forecast Date") +
    ggtitle(glue("Posterior Parameter Distirbutions by Forecast time for {target_model_name}"))
}

plot_vector_generated_quantities_by_forecast_time <- function(target_model_name) {
  tidy_generated_quantities_tbl %>%
    filter(!is.na(date),
           distribution == "posterior") %>%
    filter(model_name == target_model_name) %>%
    distinct(name, date, value, max_t) %>%
    ggplot(aes(date, value, color = max_t, group = max_t)) +
    facet_wrap(~name, scales = "free_y") +
    scale_color_viridis_c("Forecast Date") +
    geom_line(alpha = 1) +
    scale_x_date("Date") +
    scale_y_continuous("Value") +
    ggtitle(glue("Posterior Median by Forecast time for {target_model_name}")) +
    theme(legend.position = "bottom")
}

plot_single_posterior_predictive <- function(target_model_name, target_max_t) {
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    filter(distribution == "posterior",
           model_name == target_model_name,
           max_t == target_max_t) %>%
    bind_rows(tidy_generated_quantities_tbl %>%
                filter(distribution == "posterior",
                       model_name == target_model_name,
                       max_t == target_max_t,
                       name == "prop_variant_2"))

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(date <= max(tmp_tidy_posterior_predictive$date)) %>%
    filter(name %in% unique(tmp_tidy_posterior_predictive$name)) %>%
    mutate(data_type = if_else(date <= target_max_t, "observed", "future"))

  ggplot(mapping = aes(date, value)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(data = tmp_tidy_posterior_predictive,
                    mapping = aes(ymin = .lower, ymax = .upper), step = "mid") +
    geom_point(data = tmp_dat_tidy, mapping = aes(shape = data_type)) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    my_theme +
    scale_shape_discrete(name = "Data Type", labels = str_to_title) +
    scale_y_continuous(labels = comma) +
    scale_x_date("Date") +
    ggtitle(glue("{target_model_name} Posterior Predictive at t = {target_max_t}"))
}

plot_single_generated_quantities <- function(target_model_name, target_max_t) {
  tidy_generated_quantities_tbl %>%
    filter(!is.na(date),
           model_name == target_model_name,
           max_t == target_max_t,
           .width == 0.8) %>%
    ggplot(mapping = aes(date, value, ymin = .lower, ymax = .upper, color = distribution, fill = distribution)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(alpha = 0.25) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    # my_theme +
    scale_x_date("Date") +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("{target_model_name} Generated Quantities at t = {target_max_t}"),
            subtitle = "80% Credible Intervals")
}

plot_compare_posterior_generated_quantities <- function(target_max_t) {
  tidy_generated_quantities_tbl %>%
    filter(!is.na(date),
           distribution == "posterior",
           max_t == target_max_t,
           .width == 0.8) %>%
    ggplot(mapping = aes(date, value, ymin = .lower, ymax = .upper, color = model_name, fill = model_name, group = .width)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(alpha = 0.25) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    # my_theme +
    scale_x_date("Date") +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("Posterior Generated Quantities at t = {target_max_t}"),
            subtitle = "80% Credible intervals")
}

plot_peak_assessment <- function(target_peak_type){
  tmp_tidy_posterior_peak <-
    tidy_posterior_peak %>%
    filter(peak_type == target_peak_type)

  tmp_true_peak_dat <-
    true_peak_dat %>%
    filter(peak_type == target_peak_type) %>%
    expand_grid(tmp_tidy_posterior_peak %>%
                  distinct(max_t))

  ggplot(mapping = aes(max_t, value)) +
    facet_grid(name ~ model_name,
               scales = "free_y",
               labeller = labeller(
                 model_name = as_labeller(~str_replace_all(.x, "_", " ")))) +
    geom_interval(data = tmp_tidy_posterior_peak,
                  mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_tidy_posterior_peak, shape = "plus") +
    geom_point(data = tmp_true_peak_dat) +
    my_theme +
    scale_y_continuous("Value") +
    scale_x_continuous("Forecast Date") +
    ggtitle(glue("Posterior Peak {str_to_title(target_peak_type)}"))
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

compare_posterior_generated_quantities_plots <-
  tidy_predictive_tbl %>%
  distinct(max_t) %>%
  rename_with(~str_c("target_", .x)) %>%
  mutate(file_path = path("figures", experiment_name, glue("compare_posterior_generated_quantities_{target_max_t}_plot"), ext = "pdf"),
         figure = map(target_max_t, ~plot_compare_posterior_generated_quantities(target_max_t = .x))) %>%
  augment_figure_tbl()

single_posterior_predictive_plots <-
  tidy_predictive_tbl %>%
  distinct(model_name, max_t) %>%
  rename_with(~str_c("target_", .x)) %>%
  mutate(file_path = path("figures", experiment_name, glue("single_posterior_predictive_{target_model_name}_{target_max_t}_plot"), ext = "pdf"),
         figure = map2(target_model_name, target_max_t, ~plot_single_posterior_predictive(target_model_name = .x, target_max_t = .y))) %>%
  augment_figure_tbl()

peak_assessment_plots <-
  tibble(target_peak_type = all_peak_types) %>%
  mutate(file_path = path("figures", experiment_name, glue("peak_assessment_{target_peak_type}_plot"), ext = "pdf"),
         figure = future_map(target_peak_type, plot_peak_assessment)) %>%
  augment_figure_tbl()

# Save Figures ------------------------------------------------------------
dir_create(path("figures", experiment_name))
dir_create(path("figures", experiment_name, all_model_names))

all_plot_tbl_names <- ls()[str_ends(ls(), "_plots")]
single_plot_tbl_names <- all_plot_tbl_names[str_detect(all_plot_tbl_names, "single", negate = F)]
non_single_plot_tbl_names <- all_plot_tbl_names[str_detect(all_plot_tbl_names, "single", negate = T)]

# Save all figures
future_walk(all_plot_tbl_names, ~pwalk(as.list(get(.x)), ~save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, device = cairo_pdf)))

# Merge all figures
if (Sys.info()['sysname'] == "Darwin") {
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))
}

non_single_plot_tbl_names %>%
  walk(~{
    uncollected_plot_paths <- get(.x)$file_path
    collected_plot_path <- path("figures", experiment_name, str_c("all_", .x), ext = "pdf")

    str_c(uncollected_plot_paths, collapse = " ") %>%
      str_c(collected_plot_path, sep = " ") %>%
      system2("pdfunite", args = .)
  })

merge_single_plots_fn <- function(single_plot_tbl_name) {
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
}

walk(single_plot_tbl_names, merge_single_plots_fn)

# Delete unmerged figures -------------------------------------------------
walk(all_plot_tbl_names, ~file_delete(get(.x)$file_path))
