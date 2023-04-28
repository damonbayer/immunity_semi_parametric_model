time_look_for_second_wave <- 22
results_dir <- path("results", context, model_name)
sim_id <- 1

variant_2_import_time <- 40
first_obs_time <- 20

dat_tidy <- 
  path("data", context, "simulated_data", ext = "csv") %>%
  read_csv_as_draws() %>%
  filter(.iteration == sim_id) %>%
  tidy_format_draws_time() %>%
  select(name, time, value) %>%
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - 1 + (time - 1) / 7, time)) %>% 
  mutate(time = time - (first_obs_time - 1)) %>% 
  rename(true_value = value) %>%
  bind_rows(.,
            filter(., str_starts(name, "data_new_cases")) %>%
              pivot_wider(names_from = name, values_from = true_value) %>%
              mutate(name = "data_new_cases",
                     true_value = data_new_cases_variant_1 + data_new_cases_variant_2) %>%
              select(time, name, true_value)
  )

prior_generated_quantities_path <-
  path(results_dir,
       "prior_generated_quantities_csv",
       str_c("prior_generated_quantities_max_t=", target_max_t),
       ext = "csv")

posterior_generated_quantities_path <-
  path(results_dir,
       "posterior_generated_quantities_csv",
       str_c("posterior_generated_quantities_max_t=", target_max_t),
       ext = "csv")

prior_predictive_path <-
  path(results_dir,
       "prior_predictive_csv",
       str_c("prior_predictive_max_t=", target_max_t),
       ext = "csv")

posterior_predictive_path <-
  path(results_dir,
       "posterior_predictive_csv",
       str_c("posterior_predictive_max_t=", target_max_t),
       ext = "csv")


prior_generated_quantities <- read_csv_as_draws(prior_generated_quantities_path)
posterior_generated_quantities <- read_csv_as_draws(posterior_generated_quantities_path)
prior_predictive <- read_csv_as_draws(prior_predictive_path)
posterior_predictive <- read_csv_as_draws(posterior_predictive_path)

tidy_prior_generated_quantities <- tidy_generated_quantities(prior_generated_quantities) %>% mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - first_obs_time + (time - 1) / 7, time))
tidy_posterior_generated_quantities <- tidy_generated_quantities(posterior_generated_quantities) %>% mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - first_obs_time + (time - 1) / 7, time))
tidy_prior_predictive <- tidy_predictive(prior_predictive) %>%
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - first_obs_time + (time - 1) / 7, time)) %>% 
  mutate(weeks_ahead = time - target_max_t)
tidy_posterior_predictive <- tidy_predictive(posterior_predictive) %>%
  mutate(time = if_else(str_detect(name, "seq"), variant_2_import_time - first_obs_time + (time - 1) / 7, time)) %>% 
  mutate(weeks_ahead = time - target_max_t)


# Peak --------------------------------------------------------------------
est_peak_time <- 
  posterior_generated_quantities %>% 
  select(starts_with("."), matches("mean")) %>% 
  select(-matches("seq")) %>% 
  tidy_format_draws_time() %>% 
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

tidy_posterior_peak <- 
  posterior_predictive %>%
  tidy_format_draws_time() %>% 
  select(.draw, name, value, time) %>% 
  right_join(est_peak_time %>% rename(time = peak_time)) %>% 
  select(name, value, time) %>% 
  pivot_longer(-name, names_to = "peak_type") %>% 
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
  select(sample = .draw, target_type, phi) %>% mutate(phi = if_else(phi > 1e4, 1e4, phi))

mean_draws <-
  posterior_generated_quantities %>%
  select(starts_with("."), matches("mean")) %>%
  tidy_format_draws_time() %>%
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
  select(sample = .draw, prediction = value, name, time) %>%
  mutate(weeks_ahead = time - target_max_t) %>%
  filter(weeks_ahead > 0) %>%
  left_join(dat_tidy) %>%
  mutate(max_t = target_max_t) %>%
  mutate(model = model_name) %>%
  select(time, model, target_type = name, max_t, weeks_ahead, prediction, sample, true_value)

posterior_predictive_score_nbinom <-
  posterior_predictive_score_setup %>%
  left_join(mean_phi_draws) %>%
  mutate(nbinom = pmap(
    list(true_value, phi, mean),
    ~c(
      crps_nbinom = scoringRules::crps_nbinom(y = ..1, size = ..2, mu = ..3),
      logs_nbinom = scoringRules::logs_nbinom(y = ..1, size = ..2, mu = ..3),
      dss_nbinom = scoringRules::dss_nbinom(y = ..1, size = ..2, mu = ..3)))) %>%
  unnest_wider(nbinom) %>%
  group_by(time, model, target_type, max_t, weeks_ahead) %>%
  summarize(across(matches("nbinom"), mean), .groups = "drop")

posterior_predictive_score_kd <-
  posterior_predictive_score_setup %>%
  score_with_override(override = "continuous") %>%
  as_tibble()

posterior_predictive_score <- full_join(posterior_predictive_score_kd,
                                        posterior_predictive_score_nbinom)


dir_create(path(results_dir, "tidy_prior_generated_quantities"))
dir_create(path(results_dir, "tidy_posterior_generated_quantities"))
dir_create(path(results_dir, "tidy_prior_predictive"))
dir_create(path(results_dir, "tidy_posterior_predictive"))
dir_create(path(results_dir, "posterior_predictive_score"))
dir_create(path(results_dir, "tidy_posterior_peak"))

write_csv(tidy_prior_generated_quantities, file = path(results_dir, "tidy_prior_generated_quantities", str_c("tidy_prior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_generated_quantities, file = path(results_dir, "tidy_posterior_generated_quantities", str_c("tidy_posterior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_prior_predictive, file = path(results_dir, "tidy_prior_predictive", str_c("tidy_prior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_predictive, file = path(results_dir, "tidy_posterior_predictive", str_c("tidy_posterior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(posterior_predictive_score, file = path(results_dir, "posterior_predictive_score", str_c("posterior_predictive_score_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_peak, file = path(results_dir, "tidy_posterior_peak", str_c("tidy_posterior_peak_max_t=", target_max_t), ext = "csv"))