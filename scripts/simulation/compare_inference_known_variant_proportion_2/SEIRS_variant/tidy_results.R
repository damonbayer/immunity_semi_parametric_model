target_max_t <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 28, as.integer(commandArgs(trailingOnly=T[1])))
source("src/immunity_semi_parametric_model.R")
library(fs)
context <- path("simulation", "compare_inference_known_variant_proportion_2")
model_name <- "SEIRS_variant"
results_dir <- path("results", context, model_name) 
sim_id <- 1

dat_tidy <- 
  path("data", context, "simulated_data", ext = "csv") %>% 
  read_csv_as_draws() %>% 
  filter(.iteration == sim_id) %>% 
  tidy_format_draws_time() %>% 
  select(name, time, value) %>% 
  mutate(time = time - min(time) + 1) %>% 
  rename(true_value = value)

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
tidy_prior_prop_variant_2 <- tidy_prop_variant_2(prior_generated_quantities)
tidy_posterior_prop_variant_2 <- tidy_prop_variant_2(posterior_generated_quantities)

posterior_predictive_score <- 
  posterior_predictive %>% 
  tidy_format_draws_time() %>% 
  select(sample = .draw, prediction = value, name, time) %>% 
  mutate(weeks_ahead = time - target_max_t) %>%
  filter(weeks_ahead > 0) %>% 
  left_join(dat_tidy) %>%
  mutate(max_t = target_max_t) %>%
  mutate(model = model_name) %>% 
  select(time, model, target_type = name, max_t, weeks_ahead, prediction, sample, true_value) %>%
  score_with_override(override = "continuous")


dir_create(path(results_dir, "tidy_prior_generated_quantities"))
dir_create(path(results_dir, "tidy_posterior_generated_quantities"))
dir_create(path(results_dir, "tidy_prior_predictive"))
dir_create(path(results_dir, "tidy_posterior_predictive"))
dir_create(path(results_dir, "tidy_prior_prop_variant_2"))
dir_create(path(results_dir, "tidy_posterior_prop_variant_2"))
dir_create(path(results_dir, "posterior_predictive_score"))

write_csv(tidy_prior_generated_quantities, file = path(results_dir, "tidy_prior_generated_quantities", str_c("tidy_prior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_generated_quantities, file = path(results_dir, "tidy_posterior_generated_quantities", str_c("tidy_posterior_generated_quantities_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_prior_predictive, file = path(results_dir, "tidy_prior_predictive", str_c("tidy_prior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_predictive, file = path(results_dir, "tidy_posterior_predictive", str_c("tidy_posterior_predictive_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_prior_prop_variant_2, file = path(results_dir, "tidy_prior_prop_variant_2", str_c("tidy_prior_prop_variant_2_max_t=", target_max_t), ext = "csv"))
write_csv(tidy_posterior_prop_variant_2, file = path(results_dir, "tidy_posterior_prop_variant_2", str_c("tidy_posterior_prop_variant_2_max_t=", target_max_t), ext = "csv"))
write_csv(posterior_predictive_score, file = path(results_dir, "posterior_predictive_score", str_c("posterior_predictive_score_max_t=", target_max_t), ext = "csv"))


