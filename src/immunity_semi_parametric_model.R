library(tidyverse)
library(tidybayes)
library(posterior)
library(scoringutils)
library(cowplot)
library(scales)

theme_set(theme_minimal_grid())


my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  scale_color_brewer(name = "Credible Interval Width",
                     labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE),
         color = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

my_theme_no_fill <- list(
  guides(fill = guide_legend(reverse = TRUE),
         color = guide_legend(reverse = TRUE),
         linetype = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

brewer_line_color <- "#08519c"

read_csv_as_draws <- function(file) {
  file %>% 
    read_csv() %>% 
    rename_with(~str_c(".", .x), c(iteration, chain)) %>% 
    as_draws()
}

tidy_format_draws_time <- function(draws) {
  draws %>% 
    pivot_longer(-starts_with("."),
               names_to = "name_raw") %>% 
    mutate(name = name_raw %>% 
             str_extract(".+(?=\\[\\d+\\])"),
           time = name_raw %>% 
             str_extract("(?<=\\[)\\d+(?=\\])") %>% 
             as.integer()) %>% 
    mutate(name = if_else(is.na(time), name_raw, name))
}

meta_columns <-
  c(
    "lp",
    "n_steps",
    "is_accept",
    "acceptance_rate",
    "log_density",
    "hamiltonian_energy",
    "hamiltonian_energy_error",
    "max_hamiltonian_energy_error",
    "tree_depth",
    "numerical_error",
    "step_size",
    "nom_step_size"
  )

tidy_predictive <- function(predictive, ci_widths = c(0.5, 0.8, 0.95)) {
  predictive %>% 
    tidy_format_draws_time() %>% 
    select(name, time, value) %>% 
    group_by(name, time) %>% 
    median_qi(.width = ci_widths)
}

tidy_generated_quantities <- function(generated_quantities, ci_widths = c(0.5, 0.8, 0.95)) {
  generated_quantities %>% 
    tidy_format_draws_time() %>% 
    select(name, time, value) %>% 
    group_by(name, time) %>% 
    median_qi(.width = ci_widths) %>% 
    mutate(time = if_else(str_ends(name, "_mean"), time, time - 1))
}

tidy_prop_variant_2 <- function(generated_quantities, ci_widths = c(0.5, 0.8, 0.95)) {
  generated_quantities %>% 
    select(starts_with("dur_waning_α"), starts_with("dur_waning_shape")) %>% 
    expand_grid(prop_variant_2 = seq(0, 1, length.out = 100)) %>% 
    mutate(dur_waning = exp(`dur_waning_α₀` + `dur_waning_α₁` * prop_variant_2^dur_waning_shape * (1 - prop_variant_2)^dur_waning_shape)) %>% 
    select(prop_variant_2, dur_waning) %>% 
    group_by(prop_variant_2) %>% 
    median_qi(.width = ci_widths)
}

score_with_override <- function (data, metrics = NULL, override = NULL, ...) {
  check_data <- check_forecasts(data)
  data <- check_data$cleaned_data
  prediction_type <- check_data$prediction_type
  forecast_unit <- check_data$forecast_unit
  target_type <- check_data$target_type
  metrics <- scoringutils:::check_metrics(metrics)
  
  if (!is.null(override)) {
    prediction_type <- override
    target_type <- override
  }
  
  if (target_type == "binary") {
    scores <- scoringutils:::score_binary(data = data, forecast_unit = forecast_unit,
                                          metrics = metrics)
  }
  if (prediction_type == "quantile") {
    scores <- scoringutils:::score_quantile(data = data, forecast_unit = forecast_unit,
                                            metrics = metrics, ...)
  }
  if (prediction_type %in% c("integer", "continuous") && (target_type !=
                                                          "binary")) {
    scores <- scoringutils:::score_sample(data = data, forecast_unit = forecast_unit,
                                          metrics = metrics, prediction_type = prediction_type)
  }
  return(scores)
}

gg_facet_dims <- function(p) apply(X = ggplot_build(p)$layout$layout[,c("ROW", "COL")], MARGIN = 2, max)
