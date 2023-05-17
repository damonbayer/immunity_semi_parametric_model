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

model_table <- read_csv(path(simulation_dir, "model_table.csv"))

posterior_samples_summary_tbl <- 
  model_table %>%
  mutate(posterior_samples_summary_file = path(results_dir, "posterior_samples_csv", str_c("posterior_samples_fit_id=", fit_id), ext = "csv")) %>% 
  mutate(posterior_samples_summary = 
           future_map(posterior_samples_summary_file,
                    ~read_csv_as_draws(.x) %>%
                      select(-all_of(meta_columns)) %>% 
                      summarize_draws())) %>% 
  select(fit_id, posterior_samples_summary) %>% 
  unnest(posterior_samples_summary)
  
fit_ids_with_bad_convergence <- 
  posterior_samples_summary_tbl %>% 
  filter(rhat > 1.2) %>% 
  distinct(fit_id) %>% 
  pull(fit_id)

chain_to_exclude_tbl <- 
  tibble(fit_id = fit_ids_with_bad_convergence) %>% 
  mutate(posterior_samples_summary_file = path(str_c("results/simulation/sensitivity_analysis_cross_immunity_transmissibility/posterior_samples_csv/posterior_samples_fit_id=", fit_id), ext = "csv")) %>% 
  mutate(posterior_samples = future_map(posterior_samples_summary_file,
                                      ~read_csv_as_draws(.x) %>%
                                        select(-all_of(meta_columns)))) %>% 
  select(fit_id, posterior_samples) %>% 
  mutate(chain_to_exclude = map(posterior_samples, ~unique(.x$.chain))) %>% 
  unnest(chain_to_exclude) %>% 
  mutate(max_rhat = map2_dbl(posterior_samples, chain_to_exclude, ~{
    .x %>% 
      filter(.chain != .y) %>% 
      summarize_draws() %>% 
      pull(rhat) %>% 
      max()
  })) %>% 
  group_by(fit_id) %>% 
  filter(max_rhat == min(max_rhat)) %>% 
  ungroup() %>% 
  select(fit_id, chain_to_exclude, max_rhat)

write_csv(chain_to_exclude_tbl, path(simulation_dir, "chain_to_exclude", ext = "csv"))
          