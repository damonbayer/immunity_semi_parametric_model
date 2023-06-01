library(tidyverse)
library(glue)
library(fs)
library(furrr)
library(ggh4x)
library(ggbeeswarm)
source("src/immunity_semi_parametric_model.R")
options("parallelly.fork.enable" = TRUE)
plan(multicore)

experiment_name <- "ba1_immunity_contact_cdr"
context <- path("real_data", experiment_name)
results_dir <- path("results", context)
data_dir <- path("data", "real_data", "ba1")
simulation_dir <- path("scripts", context)

source(path(data_dir, "shared_constants.txt"))

date_time_0 <- ymd(date_time_0)
time_to_date <- function(time) date_time_0 + time  * 7

# Loading Data ------------------------------------------------------------
model_table <- 
  read_csv(path(simulation_dir, "model_table.csv")) %>%
  filter(immunity_model != "seq-informed") %>% 
  mutate(across(
    c(immunity_model, `R₀_model`, CDR_model),
    ~ case_when(
      .x == "gmrf" ~ "GMRF",
      .x == "seq-informed" ~ "Genetic",
      .x == "seq-informed-bin" ~ "Genetic",
      .x == "constant" ~ "Constant",
      TRUE ~ .x
    )
  )) %>%
  mutate(
    model_description_one_line = glue("R0_{`R₀_model`}_imm_{immunity_model}") %>% as.character(),
    # model_description = glue("$R_0(t)$: {`R₀_model`}, $1/\\kappa(t)$: {immunity_model}") %>% TeX(output = "character")
    model_description = glue("atop(R[0](t):~{`R₀_model`}, 1 / kappa(t):~{immunity_model})")
  ) %>% 
  mutate(variant_2_import_time = variant_2_import_time,
         time_look_for_second_wave = time_look_for_second_wave) %>% 
  # mutate(county = if_else(county == "California", county, glue("{county} County")) %>% fct_rev(county)) %>% 
  filter(!(immunity_model == "Genetic" & `R₀_model` == "Genetic"))

all_county_ids <- unique(model_table$county_id)

dat_tidy <- 
  bind_rows(
    path(data_dir, "cdph_data", ext = "csv") %>% 
      read_csv() %>% 
      filter(id %in% all_county_ids) %>% 
      select(county,
             time,
             date = end_date,
             data_new_cases = cases,
             data_new_deaths = deaths,
             data_hospitalizations = hospitalized,
             data_icu = icu) %>% 
      pivot_longer(-c(time, date, county)),
    path(data_dir, "seq_data", ext = "csv") %>% 
      read_csv() %>% 
      filter(id %in% all_county_ids) %>% 
      select(county,
             time,
             date,
             data_new_seq_variant_2 = lineage_count,
             data_new_seq_variant_1 = other_count) %>% 
      pivot_longer(-c(county, time, date))) %>% 
  mutate(variant_2_import_time = variant_2_import_time,
         time_look_for_second_wave = time_look_for_second_wave) %>% 
  mutate(county = fct_rev(county))

true_peak_dat <-
  dat_tidy %>%
  select(-date) %>% 
  filter(time > time_look_for_second_wave,
         !str_detect(name, "variant")) %>%
  select(-c(variant_2_import_time, time_look_for_second_wave)) %>% 
  group_by(name, county) %>%
  filter(value == max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-c(name, county), names_to = "peak_type") %>% 
  filter(name == "data_hospitalizations")

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
  select(-c(result_type, file_path, fit_id, date)) %>%
  pivot_longer(cols = -c(distribution, model, target_type, weeks_ahead)) %>%
  mutate(forecast_horizon = str_c(weeks_ahead, " Week Horizon")) %>%
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
all_counties <- unique(model_table$county)

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
plot_forecast_comparison <- function(target_type, target_county) {
  tmp_tidy_posterior_predictive <-
    tidy_predictive_tbl %>%
    filter(
      distribution == "posterior",
      name == target_type,
      weeks_ahead %in% c(1, 2, 4)
    ) %>%
    left_join(model_table, by = "fit_id") %>%
    filter(county == target_county)
  
  tmp_dat_tidy <-
    dat_tidy %>%
    filter(name == target_type,
           county == target_county)
  
  ggplot(mapping = aes(date, value)) +
    facet_grid(weeks_ahead ~ model_description,
               scale = "free_y",
               labeller = labeller(
                 weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead")),
                 model_description = label_parsed
               )
    ) +
    geom_lineribbon(
      data = tmp_tidy_posterior_predictive,
      mapping = aes(ymin = .lower, ymax = .upper), step = "mid",
      color = brewer_line_color
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
    scale_x_date(name = "Date", date_labels = "%b '%y", date_breaks = "1 month", guide = guide_axis(angle = 90)) +
    my_theme +
    ggtitle(glue("Forecast Comparison for {my_sim_labeller[target_type]}"),
            subtitle = glue("{county_labeller(target_county)} Data")
    )
}

plot_crps_comparison <- function(target_target_type) {
  tmp_dat <- 
    tidy_posterior_predictive_score_tbl %>%
    filter(
      target_type == target_target_type,
      weeks_ahead %in% c(1, 2, 4),
      name == "crps_nbinom"
    ) %>%
    left_join(model_table) %>%
    mutate(max_date = time_to_date(max_t)) %>% 
    select(weeks_ahead, value, max_date, forecast_horizon, model_description, county)
  
  ggplot(tmp_dat, aes(max_date, value, color = model_description)) +
    facet_grid2(forecast_horizon~county,
                independent = "y",
                scales = "free",
                labeller = labeller(county = ~glue("{county_labeller(.x)} Data"))
    ) +
    # stat_smooth(
    #   geom = "line",
    #   formula = y ~ 1,
    #   method = "lm",
    #   se = F,
    #   linetype = "dashed",
    #   alpha = 1.0
    # ) +
    geom_line() +
    geom_point() +
    scale_x_date(name = "Forecast Date", date_labels = "%b %d", breaks = unique(tmp_dat$max_date)) +
    scale_y_continuous("CRPS", labels = comma) +
    scale_color_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for {my_sim_labeller[target_target_type]}")) +
    theme(
      legend.position = "bottom",
      legend.text.align = 0
    )
}

plot_crps_comparison_boxplot <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(
      target_type == target_target_type,
      weeks_ahead %in% c(1, 2, 4),
      name == "crps_nbinom"
    ) %>%
    left_join(model_table) %>%
    select(weeks_ahead, value, max_t, forecast_horizon, model_description, county) %>%
    ggplot(aes(model_description, value, color = model_description)) +
    facet_grid2(forecast_horizon~county,
                independent = "y",
                scales = "free",
                labeller = labeller(county = ~glue("{county_labeller(.x)} Data"))
    ) +
    geom_boxplot(show.legend = F) +
    geom_beeswarm(alpha = 0.5, show.legend = F) +
    scale_y_continuous("CRPS", labels = comma) +
    scale_x_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for {my_sim_labeller[target_target_type]}"))
}

plot_crps_comparison_dotplot <- function(target_target_type) {
  tidy_posterior_predictive_score_tbl %>%
    filter(
      target_type == target_target_type,
      # weeks_ahead %in% c(1, 2, 4),
      name == "crps_nbinom"
    ) %>%
    left_join(model_table) %>%
    group_by(weeks_ahead, model_description, county) %>% 
    summarize(mean_crps = mean(value)) %>% 
    ggplot(aes(weeks_ahead, mean_crps, color = model_description)) +
    facet_wrap(~county,
               scales = "free",
               ncol = 1,
               labeller = labeller(county = ~glue("{county_labeller(.x)} Data"))
    ) +
    geom_line(linetype = "dashed", alpha = 0.5) +
    geom_point(size = 3) +
    scale_x_continuous("Forecast Horizon (Weeks)") +
    scale_y_continuous("Mean CRPS", labels = comma) +
    scale_color_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for {my_sim_labeller[target_target_type]}")) +
    theme(legend.position = "bottom")
}
  
plot_peak_assessment <- function(target_peak_type) {
  true_peak_time <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == "time"
    ) %>% 
    mutate(county = fct_rev(county)) %>% 
    select(county, max_t = value) %>% 
    mutate(max_date = time_to_date(max_t))
  
  true_peak_value <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == target_peak_type
    ) %>% 
    mutate(county = fct_rev(county)) %>% 
    select(county, value)
  
  tmp_tidy_posterior_peak <-
    tidy_posterior_peak %>%
    filter(peak_type == target_peak_type) %>%
    left_join(model_table, by = "fit_id") %>% 
    mutate(max_date = time_to_date(max_t))
  
  tmp_peak_annotation <- 
    tmp_tidy_posterior_peak %>% 
    group_by(model_description, county) %>% 
    summarize(max_value = max(.upper),
              min_t = min(max_t),
              .groups = "drop") %>%
    left_join(true_peak_time) %>% 
    left_join(true_peak_value) %>% 
    mutate(min_date = time_to_date(min_t)) %>% 
    mutate(label_hline = glue("True Peak {if_else(target_peak_type == 'time', 'Date', str_to_title(target_peak_type))}\n"),
           label_vline = "True Peak Date ") %>% 
    group_by(county) %>% 
    mutate(max_value = max(max_value)) %>% 
    filter(county == "California") %>% 
    filter(model_description == model_description[2])
  
  tmp_true_peak_dat <-
    true_peak_dat %>%
    filter(peak_type == target_peak_type) %>%
    expand_grid(tmp_tidy_posterior_peak %>%
                  distinct(max_t)) %>%
    mutate(county = fct_rev(county)) %>% 
    filter(name %in% tmp_tidy_posterior_peak$name) %>% 
    mutate(max_date = time_to_date(max_t))
  
  ggplot(tmp_tidy_posterior_peak %>% 
           filter(.width %in% c(0.8, 0.95)),
         aes(max_date, value)) +
    facet_grid2(county~model_description,
                labeller = labeller(county = ~glue("{county_labeller(.x)} Data"), model_description = label_parsed),
                scales = "free",
                independent = "x"
    ) +
    geom_vline(data = true_peak_time, mapping = aes(xintercept = max_date), linetype = "dashed") +
    geom_text(data = tmp_peak_annotation, mapping = aes(x = max_date, y = max_value, label = label_vline, hjust = "inward", vjust = "inward")) +
    geom_hline(data = true_peak_value, mapping = aes(yintercept = value), linetype = "dashed") +
    geom_text(data = tmp_peak_annotation, mapping = aes(x = min_date, y = value, label = label_hline, hjust = "inward", vjust = "inward")) +
    geom_pointinterval(mapping = aes(ymin = .lower, ymax = .upper)) +
    scale_y_continuous(glue("Peak {str_to_title(target_peak_type)}"), labels = comma) +
    scale_x_date(name = "Forecast Date", date_labels = "%b %d", breaks = unique(tmp_tidy_posterior_peak$max_date)) +
    ggtitle(glue("Posterior Peak Hospital Occupancy {str_to_title(target_peak_type)}"))
}

plot_peak_crps <- function(x = NULL) {
  true_peak_time <-
    true_peak_dat %>%
    filter(
      name == "data_hospitalizations",
      peak_type == "time"
    ) %>% 
    mutate(county = fct_rev(county)) %>%
    select(county, max_t = value) %>% 
    mutate(max_date = time_to_date(max_t))
  
  tmp_dat <- 
    tidy_posterior_peak_score %>%
    filter(name == "crps") %>%
    rename(fit_id = model) %>%
    left_join(model_table) %>%
    mutate(max_date = time_to_date(max_t))
  
  ggplot(tmp_dat, aes(max_date, value, color = model_description)) +
    facet_grid2(target_type~county,
                scales = "free", independent = "y",
               labeller = labeller(county = ~glue("{county_labeller(.x)} Data"),
                                   target_type = ~glue("Peak {str_to_title(.x)}"))) +
    # stat_smooth(
    #   geom = "line",
    #   formula = y ~ 1,
    #   method = "lm",
    #   se = F,
    #   linetype = "dashed",
    #   alpha = 1
    # ) +
    geom_text(data = true_peak_time %>% 
                filter(county == "California") %>% 
                mutate(target_type = "time"),
              mapping = aes(x = max_date, y = 1.5, label = "\nTrue Peak Date "), inherit.aes = F, vjust = "inward", hjust = "inward") +
    geom_vline(data = true_peak_time, mapping = aes(xintercept = max_date), linetype = "dashed") +
    geom_line() +
    geom_point() +
    scale_x_date(name = "Forecast Date", date_labels = "%b %d", breaks = unique(tmp_dat$max_date)) +
    scale_y_continuous("CRPS", labels = comma) +
    scale_color_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for Peak Hospitalization")) +
    theme(legend.position = "bottom")
}

plot_peak_crps_boxplot <- function(x = NULL) {
  tidy_posterior_peak_score %>%
    filter(name == "crps") %>%
    rename(fit_id = model) %>%
    left_join(model_table) %>%
    ggplot(aes(model_description, value, color = model_description)) +
    facet_grid2(target_type~county,
                scales = "free", independent = "y",
                labeller = labeller(county = ~glue("{county_labeller(.x)} Data"),
                                    target_type = ~glue("Peak {str_to_title(.x)}"))) +
    geom_boxplot(show.legend = F) +
    geom_beeswarm(alpha = 0.5, show.legend = F) +
    scale_y_continuous("CRPS", labels = comma) +
    scale_x_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for Peak Hospitalization")) +
    theme(
      legend.position = "bottom",
      legend.text.align = 0
    )
}

plot_peak_crps_dotplot <- function(x = NULL) {
  tidy_posterior_peak_score %>%
    filter(name == "crps") %>%
    rename(fit_id = model) %>%
    left_join(model_table) %>%
    group_by(model_description, county, target_type) %>% 
    summarize(mean_crps = mean(value)) %>% 
    ggplot(aes(model_description, mean_crps)) +
    facet_grid2(target_type~county,
                scales = "free", independent = "y",
                labeller = labeller(county = ~glue("{county_labeller(.x)} Data"),
                                    target_type = ~glue("Peak {str_to_title(.x)}"))) +
    geom_point(size = 4) +
    scale_y_continuous("Mean CRPS", labels = comma) +
    scale_x_discrete("Model", labels = label_parse()) +
    ggtitle(glue("Continuous Ranked Probability Score for Peak Hospitalization"))
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
         target_county = all_counties) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_forecast_comparison_{target_type}_{target_county}_plot"), ext = "pdf"),
    figure = future_map2(target_type, target_county, plot_forecast_comparison)
  ) %>%
  augment_figure_tbl()

crps_comparison_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_crps_comparison_{target_type}_plot"), ext = "pdf"),
    figure = future_map(target_type, plot_crps_comparison)
  ) %>%
  augment_figure_tbl()

crps_comparison_boxplot_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_crps_comparison_boxplot_{target_type}_plot"), ext = "pdf"),
    figure = future_map(target_type, plot_crps_comparison_boxplot)
  ) %>%
  augment_figure_tbl()

crps_comparison_dotplot_plots <-
  tibble(target_type = all_target_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_crps_comparison_dotplot_{target_type}_plot"), ext = "pdf"),
    figure = future_map(target_type, plot_crps_comparison_dotplot)
  ) %>%
  augment_figure_tbl()

peak_assessment_plots <-
  tibble(peak_type = all_peak_types) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_peak_assessment_{peak_type}_plot"), ext = "pdf"),
    figure = future_map(peak_type, plot_peak_assessment)
  ) %>%
  augment_figure_tbl()

peak_crps_plots <-
  tibble(x = 1) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_peak_crps_plot"), ext = "pdf"),
    figure = future_map(x, plot_peak_crps)
  ) %>%
  augment_figure_tbl()

peak_crps_boxplot_plots <-
  tibble(x = 1) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_peak_crps_boxplot_plot"), ext = "pdf"),
    figure = future_map(x, plot_peak_crps_boxplot)
  ) %>%
  augment_figure_tbl()

peak_crps_dotplot_plots <-
  tibble(x = 1) %>%
  mutate(
    file_path = path(manuscript_figure_dir, glue("real_data_peak_crps_dotplot_plot"), ext = "pdf"),
    figure = future_map(x, plot_peak_crps_dotplot)
  ) %>%
  augment_figure_tbl()

# Save figures ------------------------------------------------------------
forecast_comparison_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 1.75, base_height = 2.25))

crps_comparison_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2, base_height = 2.5))

crps_comparison_boxplot_plots %>% 
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2, base_heigh = 2.5))

crps_comparison_dotplot_plots %>% 
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2.25))

peak_assessment_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2, base_height = 2.5))

peak_crps_plots %>%
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2, base_height = 2.5))

peak_crps_boxplot_plots %>% 
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 2, base_height = 2.5))

peak_crps_dotplot_plots %>% 
  as.list() %>%
  pwalk(~ save_plot(filename = ..1, plot = ..2, ncol = ..3, nrow = ..4, base_asp = 1.5))
