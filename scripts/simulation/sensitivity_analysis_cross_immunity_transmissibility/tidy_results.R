target_fit_id <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 241, as.integer(commandArgs(trailingOnly=T[1])))
source("src/immunity_semi_parametric_model.R")
library(fs)

experiment_name <- "sensitivity_analysis_cross_immunity_transmissibility"
context <- path("simulation", experiment_name)
results_dir <- path("results", context)
data_dir <- path("data", context)
simulation_dir <- path("scripts", context)

source(path(simulation_dir, "shared_constants.txt"))
source(path(simulation_dir, "computed_shared_constants.txt"))

fit_info <- read_csv(path(simulation_dir, "model_table.csv")) %>% filter(fit_id == target_fit_id)


chain_to_exclude <- read_csv(path(simulation_dir, "chain_to_exclude.csv")) %>% filter(fit_id == target_fit_id) %>% pull(chain_to_exclude)

if (length(chain_to_exclude) == 0) {
  chain_to_exclude <- -1  
}

variant_2_import_time <- case_when(
  fit_info[["data_takeover_speed"]] == "slow" ~ variant_2_import_time_slow,
  fit_info[["data_takeover_speed"]] == "medium" ~ variant_2_import_time_medium,
  fit_info[["data_takeover_speed"]] == "fast" ~ variant_2_import_time_fast,
)

time_look_for_second_wave <- case_when(
  fit_info[["data_takeover_speed"]] == "slow" ~ time_look_for_second_wave_slow,
  fit_info[["data_takeover_speed"]] == "medium" ~ time_look_for_second_wave_medium,
  fit_info[["data_takeover_speed"]] == "fast" ~ time_look_for_second_wave_fast,
)

dat_tidy <- 
  path(data_dir, str_c("simulated_data_takeover_speed=", fit_info[["data_takeover_speed"]]), ext = "csv") %>%
  read_csv_as_draws() %>% 
  filter(.iteration == fit_info[["sim_id"]]) %>% 
  tidy_format_draws_time() %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time - first_obs_time + 1)) %>%
  select(time, name, true_value = value)


prior_generated_quantities_path <-
  path(results_dir,
       "prior_generated_quantities_csv",
       str_c("prior_generated_quantities_fit_id=", fit_info[["fit_id"]]),
       ext = "csv")

posterior_generated_quantities_path <-
  path(results_dir,
       "posterior_generated_quantities_csv",
       str_c("posterior_generated_quantities_fit_id=", fit_info[["fit_id"]]),
       ext = "csv")

prior_predictive_path <-
  path(results_dir,
       "prior_predictive_csv",
       str_c("prior_predictive_fit_id=", fit_info[["fit_id"]]),
       ext = "csv")

posterior_predictive_path <-
  path(results_dir,
       "posterior_predictive_csv",
       str_c("posterior_predictive_fit_id=", fit_info[["fit_id"]]),
       ext = "csv")


seq_hospitalizations_path <- path(results_dir,
     "posterior_predictive_seq_hospitalizations_csv",
     str_c("posterior_predictive_seq_hospitalizations_fit_id=", fit_info[["fit_id"]]),
     ext = "csv")

seq_hospitalization_mean_path <- 
  path(results_dir,
     "posterior_generated_quantities_seq_hospitalization_mean_csv",
     str_c("posterior_generated_quantities_seq_hospitalization_mean_fit_id=", fit_info[["fit_id"]]),
     ext = "csv")


prior_generated_quantities <- read_csv_as_draws(prior_generated_quantities_path) %>% 
  filter(.chain != chain_to_exclude)
posterior_generated_quantities <- read_csv_as_draws(posterior_generated_quantities_path) %>% 
  filter(.chain != chain_to_exclude)
prior_predictive <- read_csv_as_draws(prior_predictive_path) %>% 
  filter(.chain != chain_to_exclude)
posterior_predictive <- read_csv_as_draws(posterior_predictive_path) %>% 
  filter(.chain != chain_to_exclude)
seq_hospitalizations <- read_csv_as_draws(seq_hospitalizations_path) %>% 
  filter(.chain != chain_to_exclude)
seq_hospitalization_mean <- read_csv_as_draws(seq_hospitalization_mean_path) %>% 
  filter(.chain != chain_to_exclude)

tidy_prior_generated_quantities <- 
  prior_generated_quantities %>% 
  tidy_generated_quantities() %>% 
  mutate(time = if_else(name == "ρ_cases_t", time + 1, time)) %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  select(name, time, everything())

tidy_posterior_generated_quantities <- 
  posterior_generated_quantities %>% 
  tidy_generated_quantities() %>% 
  mutate(time = if_else(name == "ρ_cases_t", time + 1, time)) %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  select(name, time, everything())

tidy_prior_predictive <- 
  prior_predictive %>% 
  tidy_predictive() %>%
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>%
  mutate(weeks_ahead = time - fit_info[["max_t"]]) %>% 
  select(name, time, everything())

tidy_posterior_predictive <- 
  posterior_predictive %>% 
  tidy_predictive() %>%
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>%
  mutate(weeks_ahead = time - fit_info[["max_t"]]) %>% 
  select(name, time, everything())

# Peak --------------------------------------------------------------------
est_peak_time <-
  seq_hospitalization_mean %>%
  tidy_format_draws_time() %>% 
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  filter(time > time_look_for_second_wave) %>% 
  select(.draw, value, name, time) %>%
  group_by(.draw, name) %>%
  filter(value - lag(value) > 0) %>%
  filter(value == max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  select(.draw, name, peak_time = time) %>%
  mutate(name = name %>%
           str_remove("_mean") %>%
           str_c("data_", .))

posterior_peak <-
  seq_hospitalizations %>%
  tidy_format_draws_time() %>%
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  select(.draw, name, value, time) %>% 
  right_join(est_peak_time %>% rename(time = peak_time)) %>%
  select(.draw, name, value, time) %>%
  pivot_longer(-c(.draw, name), names_to = "peak_type") %>% 
  mutate(name = str_remove(name, "seq_"))

tidy_posterior_peak <- 
  posterior_peak %>% 
  select(-starts_with(".")) %>% 
  group_by(name, peak_type) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

# Score -------------------------------------------------------------------
phi_draws <-
  posterior_generated_quantities %>%
  select(starts_with("."), starts_with("ϕ_")) %>%
  tidy_format_draws_time() %>%
  select(-time) %>%
  mutate(param_type = "phi") %>%
  mutate(target_type = name %>%
           str_remove_all("ϕ_") %>%
           if_else(. == "new_seq", "new_seq_variant_1", .)) %>%
  select(-starts_with("name")) %>%
  pivot_wider(names_from = param_type, values_from = value) %>%
  bind_rows(., filter(., target_type == "new_seq_variant_1") %>% mutate(target_type = "new_seq_variant_2")) %>%
  arrange(.draw, target_type) %>%
  select(sample = .draw, target_type, phi) %>%
  mutate(phi = if_else(phi > 1e4, 1e4, phi))

mean_draws <-
  posterior_generated_quantities %>%
  select(starts_with("."), matches("mean")) %>%
  tidy_format_draws_time() %>%
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  filter(name != "new_seq_mean") %>%
  mutate(param_type = "mean") %>%
  mutate(target_type =  str_remove_all(name, "_mean")) %>%
  select(-starts_with("name")) %>%
  pivot_wider(names_from = param_type, values_from = value) %>%
  select(sample = .draw, time, target_type, mean)


mean_phi_draws <- left_join(mean_draws, phi_draws) %>%
  mutate(target_type = str_c("data_", target_type))

posterior_predictive_score_setup <-
  posterior_predictive %>%
  tidy_format_draws_time() %>%
  mutate(time = if_else(str_detect(name, "seq"),
                        variant_2_import_time - first_obs_time + (time - 1) / 7,
                        time)) %>% 
  select(sample = .draw, prediction = value, name, time) %>%
  mutate(weeks_ahead = time - fit_info[["max_t"]]) %>%
  filter(weeks_ahead > 0) %>%
  left_join(dat_tidy) %>%
  mutate(model = fit_info[["fit_id"]]) %>%
  select(time, model, target_type = name, weeks_ahead, prediction, sample, true_value)

posterior_predictive_score_nbinom <-
  posterior_predictive_score_setup %>%
  left_join(mean_phi_draws) %>%
  drop_na() %>% 
  mutate(nbinom = pmap(
    list(true_value, phi, mean),
    ~c(
      crps_nbinom = scoringRules::crps_nbinom(y = ..1, size = ..2, mu = ..3),
      logs_nbinom = scoringRules::logs_nbinom(y = ..1, size = ..2, mu = ..3),
      dss_nbinom = scoringRules::dss_nbinom(y = ..1, size = ..2, mu = ..3)))) %>%
  unnest_wider(nbinom) %>%
  group_by(time, model, target_type, weeks_ahead) %>%
  summarize(across(matches("nbinom"), mean), .groups = "drop")

posterior_predictive_score_kd <-
  posterior_predictive_score_setup %>%
  score_with_override(override = "continuous") %>%
  as_tibble()

posterior_predictive_score <- full_join(posterior_predictive_score_kd,
                                        posterior_predictive_score_nbinom)

true_peak_dat <-
  dat_tidy %>%
  filter(time > time_look_for_second_wave,
         !str_detect(name, "variant")) %>%
  group_by(name) %>%
  rename(value = true_value) %>% 
  filter(value == max(value)) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-name, names_to = "peak_type") %>% 
  filter(name == "data_hospitalizations")

posterior_peak_score <-
  posterior_peak %>%
  rename(prediction = value) %>% 
  left_join(true_peak_dat) %>% 
  mutate(model = fit_info[["fit_id"]]) %>%
  select(model, target_type = peak_type, prediction, sample = .draw, true_value = value) %>% 
  score_with_override(override = "continuous") %>%
  as_tibble()

dir_create(path(results_dir, "tidy_prior_generated_quantities"))
dir_create(path(results_dir, "tidy_posterior_generated_quantities"))
dir_create(path(results_dir, "tidy_prior_predictive"))
dir_create(path(results_dir, "tidy_posterior_predictive"))
dir_create(path(results_dir, "posterior_predictive_score"))
dir_create(path(results_dir, "tidy_posterior_peak"))
dir_create(path(results_dir, "posterior_peak_score"))

write_csv(tidy_prior_generated_quantities, file = path(results_dir, "tidy_prior_generated_quantities", str_c("tidy_prior_generated_quantities_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(tidy_posterior_generated_quantities, file = path(results_dir, "tidy_posterior_generated_quantities", str_c("tidy_posterior_generated_quantities_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(tidy_prior_predictive, file = path(results_dir, "tidy_prior_predictive", str_c("tidy_prior_predictive_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(tidy_posterior_predictive, file = path(results_dir, "tidy_posterior_predictive", str_c("tidy_posterior_predictive_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(posterior_predictive_score, file = path(results_dir, "posterior_predictive_score", str_c("posterior_predictive_score_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(tidy_posterior_peak, file = path(results_dir, "tidy_posterior_peak", str_c("tidy_posterior_peak_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
write_csv(posterior_peak_score, file = path(results_dir, "posterior_peak_score", str_c("posterior_peak_score_fit_id=", fit_info[["fit_id"]]), ext = "csv"))
