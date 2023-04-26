target_max_t <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 28, as.integer(commandArgs(trailingOnly=T[1])))
source("src/immunity_semi_parametric_model.R")
library(fs)
context <- path("simulation", "compare_inference_unknown_variant_proportion_ihr_contact")
model_name <- "SEIRS_ihr_gmrf_R0_gmrf_immunity_gmrf"
results_dir <- path("results", context, model_name)
sim_id <- 1

dat_tidy <-
  path("data", context, "simulated_data", ext = "csv") %>%
  read_csv_as_draws() %>%
  filter(.iteration == sim_id) %>%
  tidy_format_draws_time() %>%
  select(name, time, value) %>%
  mutate(time = time - min(time) + 1) %>%
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

tidy_prior_generated_quantities <- tidy_generated_quantities(prior_generated_quantities)
tidy_posterior_generated_quantities <- tidy_generated_quantities(posterior_generated_quantities)
tidy_prior_predictive <- tidy_predictive(prior_predictive) %>%
  mutate(weeks_ahead = time - target_max_t)
tidy_posterior_predictive <- tidy_predictive(posterior_predictive) %>%
  mutate(weeks_ahead = time - target_max_t)

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

write_csv(tidy_prior_generated_quantities, file = path(results_dir, "tidy_prior_generated_quantities", str_c("tidy_prior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_generated_quantities, file = path(results_dir, "tidy_posterior_generated_quantities", str_c("tidy_posterior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_prior_predictive, file = path(results_dir, "tidy_prior_predictive", str_c("tidy_prior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_predictive, file = path(results_dir, "tidy_posterior_predictive", str_c("tidy_posterior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(posterior_predictive_score, file = path(results_dir, "posterior_predictive_score", str_c("posterior_predictive_score_max_t=", target_max_t), ext = "csv"))
