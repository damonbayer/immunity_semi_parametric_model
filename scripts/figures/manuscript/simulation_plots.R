# Keep individual forecast plots
# combine all crps and
library(tidyverse)
library(glue)
library(fs)
library(furrr)
library(ggh4x)
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
model_table <- read_csv(path(simulation_dir, "model_table.csv")) %>%
  mutate(across(
    c(immunity_model, `R₀_model`, CDR_model),
    ~ case_when(
      .x == "gmrf" ~ "GMRF",
      .x == "seq-informed" ~ "Genetic",
      .x == "constant" ~ "Constant",
      TRUE ~ .x
    )
  )) %>%
  filter(is.na(prior_takeover_speed) | prior_takeover_speed == "wide") %>%
  mutate(
    model_description_one_line = glue("R0_{`R₀_model`}_imm_{immunity_model}") %>% as.character(),
    # model_description = glue("$R_0(t)$: {`R₀_model`}, $1/\\kappa(t)$: {immunity_model}") %>% TeX(output = "character")
    model_description = glue("atop(R[0](t):~{`R₀_model`}, 1 / kappa(t):~{immunity_model})")
  ) %>% 
  mutate(
    variant_2_import_time = case_when(
      data_takeover_speed == "slow" ~ variant_2_import_time_slow,
      data_takeover_speed == "medium" ~ variant_2_import_time_medium,
      data_takeover_speed == "fast" ~ variant_2_import_time_fast),
    time_look_for_second_wave = case_when(
      data_takeover_speed == "slow" ~ time_look_for_second_wave_slow,
      data_takeover_speed == "medium" ~ time_look_for_second_wave_medium,
      data_takeover_speed == "fast" ~ time_look_for_second_wave_fast,
    )) %>% 
  mutate(data_takeover_speed = fct_rev(data_takeover_speed))

target_sim_id <- model_table[["sim_id"]][[1]]

dat_tidy <- 
  tibble(file_path = dir_ls(data_dir)) %>% 
  filter(str_detect(path_file(file_path), "simulated_data_takeover_speed=\\w+\\.csv")) %>% 
  mutate(data_takeover_speed = str_extract(path_file(file_path), "(?<=simulated_data_takeover_speed=)\\w+(?=\\.csv)")) %>% 
  mutate(data = map(file_path,
                    ~read_csv_as_draws(.x) %>%
                      filter(.iteration == target_sim_id) %>%
                      tidy_format_draws_time())) %>% 
  select(-file_path) %>% 
  unnest(data) %>% 
  mutate(
    variant_2_import_time = case_when(
      data_takeover_speed == "slow" ~ variant_2_import_time_slow,
      data_takeover_speed == "medium" ~ variant_2_import_time_medium,
      data_takeover_speed == "fast" ~ variant_2_import_time_fast),
    time_look_for_second_wave = case_when(
      data_takeover_speed == "slow" ~ time_look_for_second_wave_slow,
      data_takeover_speed == "medium" ~ time_look_for_second_wave_medium,
      data_takeover_speed == "fast" ~ time_look_for_second_wave_fast,
    )) %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time - first_obs_time + 1
  )) %>% 
  select(data_takeover_speed, time, name, value, variant_2_import_time, time_look_for_second_wave)

true_peak_dat <-
  dat_tidy %>%
  filter(
    time > time_look_for_second_wave,
    !str_detect(name, "variant")
  ) %>%
  select(-c(variant_2_import_time, time_look_for_second_wave)) %>% 
  group_by(name, data_takeover_speed) %>%
  filter(value == max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-c(name, data_takeover_speed), names_to = "peak_type")

results_tbl <-
  dir_ls(results_dir, recurse = 2, type = "file") %>%
  enframe(name = NULL, value = "file_path") %>%
  filter(path_ext(file_path) == "csv") %>%
  mutate(split_path = path_split(file_path)) %>%
  mutate(file_name = file_path %>% path_file() %>% path_ext_remove()) %>%
  mutate(result_type = map_chr(split_path, ~ pluck(.x, 4))) %>%
  filter(str_starts(result_type, "tidy_") | str_detect(result_type, "score") | str_detect(result_type, "posterior_samples") | str_detect(result_type, "duration")) %>%
  mutate(distribution = str_extract(result_type, "prior|posterior")) %>%
  mutate(result_type = str_remove(result_type, ".*(posterior|prior)_")) %>%
  mutate(fit_id = file_name %>%
    str_extract("(?<=fit_id=)\\d+") %>%
    as.integer()) %>%
  select(fit_id, result_type, distribution, file_path) %>%
  right_join(model_table %>% select(fit_id)) %>%
  arrange(fit_id)

tidy_duration_tbl <-
  results_tbl %>%
  filter(result_type == "duration") %>%
  mutate(tidy_duration = future_map(file_path, read_csv)) %>%
  unnest(tidy_duration) %>%
  select(-result_type, -file_path, -distribution) %>%
  mutate(across(c(wall, compute), ~ seconds(.) |>
    seconds_to_period() %>%
    as.duration()))

tidy_predictive_tbl <-
  results_tbl %>%
  filter(result_type == "predictive") %>%
  mutate(tidy_predictive = future_map(file_path, read_csv)) %>%
  unnest(tidy_predictive) %>%
  select(-result_type, -file_path)

tidy_generated_quantities_tbl <-
  results_tbl %>%
  filter(result_type == "generated_quantities") %>%
  mutate(tidy_generated_quantities = future_map(file_path, read_csv)) %>%
  unnest(tidy_generated_quantities) %>%
  select(-result_type, -file_path) %>%
  mutate(name = fct_reorder(name, str_detect(name, "mean") * 1 + str_detect(name, "compartment") * 2))

tidy_posterior_predictive_score_tbl <-
  results_tbl %>%
  filter(result_type == "predictive_score") %>%
  mutate(tidy_predictive_score = future_map(file_path, read_csv)) %>%
  unnest(tidy_predictive_score) %>%
  select(-c(result_type, file_path, fit_id, time)) %>%
  pivot_longer(cols = -c(distribution, model, target_type, weeks_ahead)) %>%
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Forecast Horizon")) %>%
  rename(fit_id = model)

tidy_posterior_peak <-
  results_tbl %>%
  filter(result_type == "peak") %>%
  mutate(peak = future_map(file_path, read_csv)) %>%
  unnest(peak) %>%
  select(-result_type, -file_path)

tidy_posterior_peak_score <-
  results_tbl %>%
  filter(result_type == "peak_score") %>%
  mutate(peak_score = future_map(file_path, read_csv)) %>%
  unnest(peak_score) %>%
  select(-result_type, -file_path) %>%
  select(-fit_id) %>%
  pivot_longer(cols = -c(distribution, model, target_type))

all_target_types <- unique(tidy_predictive_tbl$name)
all_target_types <- all_target_types[str_detect(all_target_types, "variant", negate = T)]
all_model_description_one_lines <- unique(model_table$model_description_one_line)
all_peak_types <- unique(tidy_posterior_peak$peak_type)
all_data_takeover_speeds <- unique(model_table$data_takeover_speed)

model_description_key <- model_table %>%
  distinct(model_description_one_line, model_description) %>%
  deframe()

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
plot_forecast_comparison <- function(target_type, target_data_takeover_speed) {
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    filter(
      distribution == "posterior",
      name == target_type,
      weeks_ahead %in% c(1, 2, 4)
    ) %>%
    left_join(model_table, by = "fit_id") %>%
    filter(is.na(prior_takeover_speed) | prior_takeover_speed == "wide",
           data_takeover_speed == target_data_takeover_speed)

  tmp_dat_tidy <-
    dat_tidy %>%
    filter(name == target_type,
           data_takeover_speed == target_data_takeover_speed)

  ggplot(mapping = aes(time, value)) +
    facet_grid(weeks_ahead ~ model_description,
      scale = "free_y",
      labeller = labeller(
        weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
        model_description = label_parsed
      )
    ) +
    geom_lineribbon(
      data = tmp_tidy_posterior_predictive,
      mapping = aes(ymin = .lower, ymax = .upper), step = "mid"
    ) +
    geom_point(data = tmp_tidy_posterior_predictive %>%
      group_by(weeks_ahead) %>%
      summarize(
        min_time = min(time),
        max_time = max(time)
      ) %>%
      expand_grid(tmp_dat_tidy) %>%
      filter(time <= max_time)) +
    scale_y_continuous(name = my_sim_labeller[target_type], labels = comma) +
    scale_x_continuous(name = "Time") +
    my_theme +
    ggtitle(glue("Forecast Comparison for {my_sim_labeller[target_type]}"),
      subtitle = glue("{str_to_title(target_data_takeover_speed)} Takeover Data")
    )
}

plot_crps_comparison <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(
      target_type == target_target_type,
      weeks_ahead %in% c(1, 2, 4),
      name == "crps_nbinom"
    ) %>%
    left_join(model_table) %>%
    select(weeks_ahead, value, max_t, forecast_horizon, model_description, data_takeover_speed) %>%
    
    ggplot(aes(max_t, value, color = model_description)) +
    facet_grid2(forecast_horizon~data_takeover_speed,
                       independent = "y",
               scales = "free",
               labeller = labeller(data_takeover_speed = ~glue("{str_to_title(.x)} Takeover Data"))
    ) +
    stat_smooth(
      geom = "line",
      formula = y ~ 1,
      method = "lm",
      se = F,
      linetype = "dashed",
      alpha = 1.0
    ) +
    geom_line(alpha = 0.8) +
    geom_point(alpha = 0.8) +
    scale_x_continuous("Forecast Time") +
    scale_y_continuous("CRPS", labels = comma) +
    scale_color_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for {my_sim_labeller[target_target_type]}")) +
    theme(
      legend.position = "bottom",
      legend.text.align = 0
    ) #+
    # guides(col = guide_legend(ncol = 1))
}

plot_crps_comparison_boxplot <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(
      target_type == target_target_type,
      weeks_ahead %in% c(1, 2, 4),
      name == "crps_nbinom"
    ) %>%
    left_join(model_table) %>%
    select(weeks_ahead, value, max_t, forecast_horizon, model_description, data_takeover_speed) %>%
    ggplot(aes(model_description, value, color = model_description)) +
    facet_grid2(forecast_horizon~data_takeover_speed,
                independent = "y",
                scales = "free",
                labeller = labeller(data_takeover_speed = ~glue("{str_to_title(.x)} Takeover Data"))
    ) +
    geom_boxplot(show.legend = F) +
    geom_beeswarm(alpha = 0.5, show.legend = F) +
    scale_y_continuous("CRPS", labels = comma) +
    scale_x_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for {my_sim_labeller[target_target_type]}")) +
    theme(
      legend.position = "bottom",
      legend.text.align = 0
    )
}
  

plot_peak_assessment <- function(target_peak_type) {
  true_peak_time <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == "time"
    ) %>% 
    mutate(data_takeover_speed = fct_rev(data_takeover_speed)) %>% 
    select(data_takeover_speed, max_t = value)
  
  true_peak_value <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == target_peak_type
    ) %>% 
    mutate(data_takeover_speed = fct_rev(data_takeover_speed)) %>% 
    select(data_takeover_speed, value)
  
  tmp_tidy_posterior_peak <-
    tidy_posterior_peak %>%
    filter(peak_type == target_peak_type) %>%
    left_join(model_table, by = "fit_id")
  
  tmp_peak_annotation <- 
    tmp_tidy_posterior_peak %>% 
    group_by(model_description, data_takeover_speed) %>% 
    summarize(max_value = max(.upper),
              min_t = min(max_t),
              .groups = "drop") %>%
    left_join(true_peak_time) %>% 
    left_join(true_peak_value) %>% 
    mutate(label_hline = glue("True Peak {str_to_title(target_peak_type)}\n"),
           label_vline = "True Peak Time ") %>% 
    group_by(data_takeover_speed) %>% 
    mutate(max_value = max(max_value)) %>% 
    filter(data_takeover_speed == "fast") %>% 
    filter(model_description == model_description[2])
  
  tmp_true_peak_dat <-
    true_peak_dat %>%
    filter(peak_type == target_peak_type) %>%
    expand_grid(tmp_tidy_posterior_peak %>%
                  distinct(max_t)) %>%
    mutate(data_takeover_speed = fct_rev(data_takeover_speed)) %>% 
    filter(name %in% tmp_tidy_posterior_peak$name)
  
  ggplot(tmp_tidy_posterior_peak,
         aes(max_t, value)) +
    facet_grid2(data_takeover_speed~model_description,
                labeller = labeller(data_takeover_speed = ~glue("{str_to_title(.x)} Takeover Data"), model_description = label_parsed),
                scales = "free",
                independent = "x"
    ) +
    geom_vline(data = true_peak_time, mapping = aes(xintercept = max_t), linetype = "dashed") +
    geom_text(data = tmp_peak_annotation, mapping = aes(x = max_t, y = max_value, label = label_vline, hjust = "inward", vjust = "inward")) +
    geom_hline(data = true_peak_value, mapping = aes(yintercept = value), linetype = "dashed") +
    geom_text(data = tmp_peak_annotation, mapping = aes(x = min_t, y = value, label = label_hline, hjust = "inward", vjust = "inward")) +
    geom_pointinterval(mapping = aes(ymin = .lower, ymax = .upper)) +
    scale_y_continuous(glue("Peak {str_to_title(target_peak_type)}"), labels = comma) +
    scale_x_continuous("Forecast Time") +
    ggtitle(glue("Posterior Peak Hospital Occupancy {str_to_title(target_peak_type)}"))
}

plot_peak_crps <- function(target_peak_type) {
  true_peak_time <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == "time"
    ) %>% 
    mutate(data_takeover_speed = fct_rev(data_takeover_speed)) %>%
    select(data_takeover_speed, max_t = value)
  
  tidy_posterior_peak_score %>%
    filter(target_type == target_peak_type) %>%
    filter(name == "crps") %>%
    rename(fit_id = model) %>%
    left_join(model_table) %>%
    ggplot(aes(max_t, value, color = model_description)) +
    facet_wrap(~data_takeover_speed,
               scales = "free", ncol = 1,
               labeller = labeller(data_takeover_speed = ~glue("{str_to_title(.x)} Takeover Data"), model_description = label_parsed)) +
    stat_smooth(
      geom = "line",
      formula = y ~ 1,
      method = "lm",
      se = F,
      linetype = "dashed",
      alpha = 1
    ) +
    geom_text(data = true_peak_time, mapping = aes(x = max_t, y = Inf, label = "\nTrue Peak Time "), inherit.aes = F, vjust = "inward", hjust = "inward") +
    geom_vline(data = true_peak_time, mapping = aes(xintercept = max_t), linetype = "dashed") +
    geom_line(alpha = 0.8) +
    geom_point(alpha = 0.8) +
    geom_point() +
    scale_x_continuous("Forecast Time") +
    scale_y_continuous("CRPS", labels = comma) +
    scale_color_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for Peak {my_sim_labeller[target_target_type]} {str_to_title(target_peak_type)}")) +
    theme(legend.position = "bottom")
    # theme(legend.text.align = 0)
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
  expand_grid(target_type = all_target_types,
         target_data_takeover_speed = all_data_takeover_speeds) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("simulated_forecast_comparison_{target_type}_{target_data_takeover_speed}_plot"), ext = "pdf"),
    figure = future_map2(target_type, target_data_takeover_speed, plot_forecast_comparison)
  ) %>%
  augment_figure_tbl()

crps_comparison_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("simulated_crps_comparison_{target_type}_plot"), ext = "pdf"),
    figure = future_map(target_type, plot_crps_comparison)
  ) %>%
  augment_figure_tbl()

crps_comparison_boxplot_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("simulated_crps_comparison_boxplot_{target_type}_plot"), ext = "pdf"),
    figure = future_map(target_type, plot_crps_comparison_boxplot)
  ) %>%
  augment_figure_tbl()

peak_assessment_plots <-
  tibble(peak_type = all_peak_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("simulated_peak_assessment_{peak_type}_plot"), ext = "pdf"),
    figure = future_map(peak_type, plot_peak_assessment)
  ) %>%
  augment_figure_tbl()

peak_crps_plots <-
  tibble(peak_type = all_peak_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("simulated_peak_crps_{peak_type}_{target_data_takeover_speed}_plot"), ext = "pdf"),
    figure = future_map(peak_type, plot_peak_crps)
  ) %>%
  augment_figure_tbl()

# Save figures ------------------------------------------------------------
forecast_comparison_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 1.75))

crps_comparison_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4))

crps_comparison_boxplot_plots %>% 
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4))

peak_assessment_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4))

peak_crps_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2.75))
