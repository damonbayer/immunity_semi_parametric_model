target_data_id <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))
library(tidyverse)
library(glue)
library(fs)
library(furrr)
source("src/immunity_semi_parametric_model.R")
options("parallelly.fork.enable" = TRUE)
plan(multicore)

experiment_name <- "sensitivity_analysis_cross_immunity_transmissibility"
context <- path("simulation", experiment_name)
results_dir <- path("results", context)
data_dir <- path("data", context)
simulation_dir <- path("scripts", context)

source(path(simulation_dir, "shared_constants.txt"))
source(path(simulation_dir, "computed_shared_constants.txt"))


# Loading Data ------------------------------------------------------------
model_table_subset <- read_csv(path(simulation_dir, "model_table.csv")) %>% 
  filter(as.integer(factor(data_takeover_speed)) == target_data_id) %>% 
  mutate(across(c(immunity_model, `R₀_model`, CDR_model),
                ~case_when(.x == "seq-informed" ~ "seq",
                           .x == "constant" ~ "const"))) %>% 
  mutate(model_description = glue("R₀: {`R₀_model`}\nimm: {immunity_model}\nprior: {prior_takeover_speed}") %>% as.character(),
         model_description_one_line = glue("R0_{`R₀_model`}_imm_{immunity_model}_prior_{prior_takeover_speed}") %>% as.character())

target_data_takeover_speed <- model_table_subset[["data_takeover_speed"]][[1]]
target_sim_id <- model_table_subset[["sim_id"]][[1]]

variant_2_import_time <- case_when(
  target_data_takeover_speed == "slow" ~ variant_2_import_time_slow,
  target_data_takeover_speed == "medium" ~ variant_2_import_time_medium,
  target_data_takeover_speed == "fast" ~ variant_2_import_time_fast,
)

time_look_for_second_wave <- case_when(
  target_data_takeover_speed == "slow" ~ time_look_for_second_wave_slow,
  target_data_takeover_speed == "medium" ~ time_look_for_second_wave_medium,
  target_data_takeover_speed == "fast" ~ time_look_for_second_wave_fast,
)


dat_tidy <- 
  path(data_dir, str_c("simulated_data_takeover_speed=", target_data_takeover_speed), ext = "csv") %>%
  read_csv_as_draws() %>% 
  filter(.iteration == target_sim_id) %>% 
  tidy_format_draws_time() %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time - first_obs_time + 1)) %>%
  select(time, name, value) %>% 
  mutate(name = str_remove(name, "data+_"))

true_peak_dat <-
  dat_tidy %>%
  filter(time > time_look_for_second_wave,
         !str_detect(name, "variant")) %>%
  group_by(name) %>%
  filter(value == max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-name, names_to = "peak_type")

results_tbl <-
  dir_ls(results_dir, recurse = 2, type = "file") %>% 
  enframe(name = NULL, value = "file_path") %>%
  filter(path_ext(file_path) == "csv") %>% 
  mutate(split_path = path_split(file_path)) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  mutate(result_type = map_chr(split_path, ~pluck(.x, 4))) %>% 
  filter(str_starts(result_type, "tidy_") | str_detect(result_type, "score") | str_detect(result_type, "posterior_samples") | str_detect(result_type, "duration")) %>%
  mutate(distribution = str_extract(result_type, "prior|posterior")) %>%
  mutate(result_type = str_remove(result_type, ".*(posterior|prior)_")) %>%
  mutate(fit_id = file_name %>%
           str_extract("(?<=fit_id=)\\d+") %>%
           as.integer()) %>%
  select(fit_id, result_type, distribution, file_path) %>% 
  right_join(model_table_subset %>% select(fit_id)) %>% 
  arrange(fit_id)

tidy_duration_tbl <- 
  results_tbl %>%
  filter(result_type == "duration") %>% 
  mutate(tidy_duration = future_map(file_path, read_csv)) %>%
  unnest(tidy_duration) %>%
  select(-result_type, -file_path, -distribution) %>% 
  mutate(across(c(wall, compute), ~seconds(.) |> seconds_to_period() %>% as.duration()))

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
  mutate(tidy_predictive_score = future_map(file_path, read_csv)) %>%
  unnest(tidy_predictive_score) %>% 
  select(-c(result_type, file_path, fit_id, time)) %>% 
  pivot_longer(cols = -c(distribution, model, target_type, weeks_ahead)) %>% 
  mutate(target_type = str_remove(target_type, "data+_")) %>%
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Horizon")) %>% 
  rename(fit_id = model)

tidy_posterior_peak <-
  results_tbl %>%
  filter(result_type == "peak") %>%
  mutate(peak = future_map(file_path, read_csv)) %>%
  unnest(peak) %>%
  select(-result_type, -file_path) %>%
  mutate(name = str_remove(name, "data+_"))

all_target_types <- unique(tidy_predictive_tbl$name)
all_model_description_one_lines <- unique(model_table_subset$model_description_one_line)
all_peak_types <- unique(tidy_posterior_peak$peak_type)

model_description_key <- model_table_subset %>% 
  distinct(model_description_one_line, model_description) %>% 
  deframe()

# model_labels <- 
#   tibble(name = all_model_names) %>% 
#   mutate(value = name %>% 
#            str_remove_all("SEIRS_") %>% 
#            str_split("_") %>% 
#            map_chr(~{
#              str_c(
#                .x[seq(1, length(.x), by = 2)],
#                .x[seq(2, length(.x), by = 2)],
#                sep = ": ",
#                collapse = "\n")}) %>% 
#            str_replace_all("immunity", "imm")) %>% 
#   deframe()

forecast_metrics_to_keep <-
  tibble(metric = unique(tidy_posterior_predictive_score_tbl$name) %>% str_replace("log_score", "logs")) %>%
  mutate(ends_binom = str_ends(metric, "_nbinom")) %>% 
  mutate(clean_metric = str_remove(metric, "_nbinom")) %>%
  add_count(clean_metric, name = "replicates") %>% 
  mutate(replicates = as.logical(replicates - 1)) %>% 
  filter(!replicates | ends_binom) %>% 
  pull(metric) %>% 
  sort()


# Plot functions ----------------------------------------------------------

plot_forecast_comparison <- function(target_type) {
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    filter(distribution == "posterior",
           name == target_type,
           weeks_ahead %in% c(0, 1, 2, 4)) %>% 
    left_join(model_table_subset, by = "fit_id")

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(name == target_type)

  ggplot(mapping = aes(time, value)) +
    facet_grid(weeks_ahead ~ model_description, scale = "free_y",
               labeller = labeller(
                 weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")))) +
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
    ggtitle(glue("Forecast Comparison for {str_to_title(str_replace_all(target_type, '_', ' '))}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))
}

plot_forecast_metrics_comparison <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(target_type == target_target_type,
           weeks_ahead %in% c(1, 2, 4),
           name %in% forecast_metrics_to_keep) %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    ggplot(aes(model_description, value, color = model_description)) +
    facet_grid(name ~ forecast_horizon, scales = "free_y") +
    geom_boxplot(show.legend = F) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(name = "Value", labels = comma) +
    scale_x_discrete(name = "Model") +
    theme(legend.position = "bottom") +
    ggtitle(glue("Forecast Metrics Comparison for {str_to_title(str_replace_all(target_target_type, '_', ' '))}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))
}

plot_scalar_generated_quantities_by_forecast_time <- function(target_model_description_one_line) {
  target_model_description <- unname(model_description_key[target_model_description_one_line])
  
  available_generated_quantities_names <-
    tidy_generated_quantities_tbl %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    filter(model_description == target_model_description) %>%
    filter(name != "prop_variant_2") %>%
    pull(name) %>%
    unique()

  tidy_generated_quantities_tbl %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    mutate(max_time = max_t) %>%
    filter(is.na(time), .width == 0.8) %>%
    filter(model_description == target_model_description) %>%
    ggplot(aes(max_time, value, ymin = .lower, ymax = .upper, color = distribution)) +
    facet_wrap(~name, scales = "free_y") +
    geom_interval(alpha = 0.5) +
    scale_y_continuous("Value") +
    scale_x_continuous("Forecast Time") +
    ggtitle(glue("Posterior Parameter Distirbutions by Forecast time for {str_replace_all(target_model_description, '\n', ', ')}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))
}

plot_vector_generated_quantities_by_forecast_time <- function(target_model_description_one_line) {
  target_model_description <- unname(model_description_key[target_model_description_one_line])
  
  tidy_generated_quantities_tbl %>%
    filter(!is.na(time),
           distribution == "posterior") %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    filter(model_description == target_model_description) %>%
    distinct(name, time, value, max_t) %>%
    ggplot(aes(time, value, color = max_t, group = max_t)) +
    facet_wrap(~name, scales = "free_y") +
    scale_color_viridis_c("Forecast Time") +
    geom_line(alpha = 1) +
    scale_x_continuous("Time") +
    scale_y_continuous("Value") +
    ggtitle(glue("Posterior Median by Forecast time for {str_replace_all(target_model_description, '\n', ', ')}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data")) +
    theme(legend.position = "bottom")
}

plot_single_posterior_predictive <- function(target_model_description_one_line, target_max_t) {
  target_model_description <- unname(model_description_key[target_model_description_one_line])
  
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    filter(distribution == "posterior",
           model_description == target_model_description,
           max_t == target_max_t) %>%
    bind_rows(tidy_generated_quantities_tbl %>%
                left_join(model_table_subset, by = "fit_id") %>% 
                filter(distribution == "posterior",
                       model_description == target_model_description,
                       max_t == target_max_t,
                       name == "prop_variant_2"))

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
    scale_x_continuous("Time") +
    ggtitle(glue("{str_replace_all(target_model_description, '\n', ', ')} Posterior Predictive at t = {target_max_t}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))
}

plot_single_generated_quantities <- function(target_model_description_one_line, target_max_t) {
  target_model_description <- unname(model_description_key[target_model_description_one_line])
  
  tidy_generated_quantities_tbl %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    filter(!is.na(time),
           model_description == target_model_description,
           max_t == target_max_t,
           .width == 0.8) %>%
    ggplot(mapping = aes(time, value, ymin = .lower, ymax = .upper, color = distribution, fill = distribution)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(alpha = 0.25) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    scale_x_continuous("Time") +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("{str_replace_all(target_model_description, '\n', ', ')} Generated Quantities at t = {target_max_t}"),
            subtitle = "{str_to_title(target_data_takeover_speed)} Takeover Data, 80% Credible Intervals")
}

plot_compare_posterior_generated_quantities <- function(target_max_t) {
  tidy_generated_quantities_tbl %>%
    left_join(model_table_subset, by = "fit_id") %>% 
    filter(!is.na(time),
           distribution == "posterior",
           max_t == target_max_t,
           .width == 0.8) %>%
    ggplot(mapping = aes(time, value, ymin = .lower, ymax = .upper, color = model_description, fill = model_description, group = .width)) +
    facet_wrap(.~name, scales = "free_y") +
    geom_lineribbon(alpha = 0.25) +
    geom_vline(xintercept = target_max_t, linetype = "dashed") +
    # my_theme +
    scale_x_continuous("Time") +
    scale_y_continuous(labels = comma) +
    ggtitle(glue("Posterior Generated Quantities at t = {target_max_t}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data, 80% Credible Intervals, 80% Credible intervals"))
}

plot_peak_assessment <- function(target_peak_type){
  tmp_tidy_posterior_peak <-
    tidy_posterior_peak %>%
    filter(peak_type == target_peak_type) %>% 
    left_join(model_table_subset, by = "fit_id")

  tmp_true_peak_dat <-
    true_peak_dat %>%
    filter(peak_type == target_peak_type) %>%
    expand_grid(tmp_tidy_posterior_peak %>%
                  distinct(max_t))

  ggplot(mapping = aes(max_t, value)) +
    facet_grid(name ~ model_description,
               scales = "free_y") +
    geom_interval(data = tmp_tidy_posterior_peak,
                  mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = tmp_tidy_posterior_peak, shape = "plus") +
    geom_point(data = tmp_true_peak_dat) +
    my_theme +
    scale_y_continuous("Value") +
    scale_x_continuous("Forecast Date") +
    ggtitle(glue("Posterior Peak {str_to_title(target_peak_type)}"),
            subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))
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
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("forecast_comparison_{target_type}_plot"), ext = "pdf"),
         figure = future_map(target_type, plot_forecast_comparison)) %>%
  augment_figure_tbl()

forecast_metrics_comparison_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("forecast_metrics_comparison_{target_type}_plot"), ext = "pdf"),
         figure = future_map(target_type, plot_forecast_metrics_comparison)) %>%
  augment_figure_tbl()

scalar_generated_quantities_by_forecast_time_plots <-
  tibble(target_model_description_one_line = all_model_description_one_lines) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("scalar_generated_quantities_by_forecast_time_{target_model_description_one_line}_plot"), ext = "pdf"),
         figure = future_map(target_model_description_one_line, plot_scalar_generated_quantities_by_forecast_time)) %>%
  augment_figure_tbl()

vector_generated_quantities_by_forecast_time_plots <-
  tibble(target_model_description_one_line = all_model_description_one_lines) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("vector_generated_quantities_by_forecast_time_{target_model_description_one_line}_plot"), ext = "pdf"),
         figure = future_map(target_model_description_one_line, plot_vector_generated_quantities_by_forecast_time)) %>%
  mutate(figure_dims = future_map(figure, gg_facet_dims)) %>%
  augment_figure_tbl()

single_generated_quantities_plots <-
  tidy_predictive_tbl %>%
  left_join(model_table_subset, by = "fit_id") %>% 
  distinct(model_description_one_line, max_t) %>%
  rename_with(~str_c("target_", .x)) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("single_generated_quantities_{target_model_description_one_line}_{target_max_t}_plot"), ext = "pdf"),
         figure = map2(target_model_description_one_line, target_max_t, ~plot_single_generated_quantities(target_model_description_one_line = .x, target_max_t = .y))) %>%
  augment_figure_tbl()

compare_posterior_generated_quantities_plots <-
  tidy_predictive_tbl %>%
  left_join(model_table_subset, by = "fit_id") %>% 
  distinct(max_t) %>%
  rename_with(~str_c("target_", .x)) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("compare_posterior_generated_quantities_{target_max_t}_plot"), ext = "pdf"),
         figure = map(target_max_t, ~plot_compare_posterior_generated_quantities(target_max_t = .x))) %>%
  augment_figure_tbl()

single_posterior_predictive_plots <-
  tidy_predictive_tbl %>%
  left_join(model_table_subset, by = "fit_id") %>% 
  distinct(model_description_one_line, max_t) %>%
  rename_with(~str_c("target_", .x)) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("single_posterior_predictive_{target_model_description_one_line}_{target_max_t}_plot"), ext = "pdf"),
         figure = map2(target_model_description_one_line, target_max_t, ~plot_single_posterior_predictive(target_model_description_one_line = .x, target_max_t = .y))) %>%
  augment_figure_tbl()

peak_assessment_plots <-
  tibble(target_peak_type = all_peak_types) %>%
  mutate(file_path = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), glue("peak_assessment_{target_peak_type}_plot"), ext = "pdf"),
         figure = future_map(target_peak_type, plot_peak_assessment)) %>%
  augment_figure_tbl()

compute_time_plot <- 
  tidy_duration_tbl %>% 
  mutate(compute_hours = as.numeric(compute) / 60 / 60) %>% 
  left_join(model_table_subset, by = "fit_id") %>% 
  ggplot(aes(max_t, compute_hours, color = model_description)) +
  geom_line() +
  geom_point() +
  scale_x_continuous("Forecast Time") +
  scale_y_continuous("Compute Hours") +
  scale_color_discrete("Model") +
  ggtitle("Computation Time",
          subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data"))

# Save Figures ------------------------------------------------------------
dir_create(path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}")))
dir_create(path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), all_model_description_one_lines))


# Save the compute time plot
save_plot(filename = path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), "compute_time", ext = "pdf"),
          plot = compute_time_plot, device = cairo_pdf)

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
    collected_plot_path <- path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), str_c("all_", .x), ext = "pdf")

    str_c(uncollected_plot_paths, collapse = " ") %>%
      str_c(collected_plot_path, sep = " ") %>%
      system2("pdfunite", args = .)
  })

merge_single_plots_fn <- function(single_plot_tbl_name) {
  single_plot_tbl <- get(single_plot_tbl_name)
  single_plot_tbl %>%
    group_by(target_model_description_one_line) %>%
    group_walk(~{
      uncollected_plot_paths <- .x$file_path
      collected_plot_path <- path("figures", experiment_name, glue("takeover_speed_{target_data_takeover_speed}"), .y, str_c("all_", single_plot_tbl_name, "_", .y), ext = "pdf")

      str_c(uncollected_plot_paths, collapse = " ") %>%
        str_c(collected_plot_path, sep = " ") %>%
        system2("pdfunite", args = .)
    })
}

walk(single_plot_tbl_names, merge_single_plots_fn)

# Delete unmerged figures -------------------------------------------------
walk(all_plot_tbl_names, ~file_delete(get(.x)$file_path))
