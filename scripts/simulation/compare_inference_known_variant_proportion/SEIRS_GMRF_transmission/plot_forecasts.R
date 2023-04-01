source("src/immunity_semi_parametric_model.R")
library(fs)
context <- path("simulation", "compare_inference_known_variant_proportion")
model_name <- "SEIRS_GMRF_transmission"
results_dir <- path("results", context, model_name)
sim_id <- 1

dat_tidy <-
  path("data", context, "simulated_data", ext = "csv") %>%
  read_csv_as_draws() %>%
  filter(.iteration == sim_id) %>%
  tidy_format_draws_time() %>%
  select(name, time, value) %>%
  mutate(time = time - min(time) + 1)

posterior_predictive <-
  tibble(file_path = path(results_dir, "tidy_posterior_predictive") %>% dir_ls()) %>%
  mutate(max_t = file_path %>%
           path_file() %>%
           path_ext_remove() %>%
           str_extract("(?<=max_t=)\\d+") %>%
           as.integer()) %>%
  mutate(posterior_predictive = map(file_path, read_csv)) %>%
  select(-file_path) %>%
  unnest(posterior_predictive) %>%
  filter(weeks_ahead %in% c(0, 1, 4))


weeks_ahead_forecast_plot <-
  ggplot(mapping = aes(time, value)) +
  facet_grid(name ~ weeks_ahead, scales = "free_y",
             labeller = labeller(weeks_ahead = as_labeller(~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
                                 name = as_labeller(~str_remove_all(.x, "data_") %>% str_replace_all("_", " ") %>% str_to_title))) +
  geom_lineribbon(data = posterior_predictive,
                  mapping = aes(ymin =.lower, ymax = .upper), step = "mid") +
  geom_point(data = posterior_predictive %>%
               group_by(weeks_ahead) %>%
               summarize(min_time = min(time),
                         max_time = max(time)) %>%
               expand_grid(dat_tidy) %>%
               filter(time <= max_time)) +
  scale_y_continuous(name = "Value", labels = comma) +
  scale_x_continuous(name = "Time") +
  my_theme +
  ggtitle("SEIRS GMRF Transmission Forecasts")

save_plot(filename = path("figures", "SEIRS_GMRF_transmission_weeks_ahead_forecast_plot", ext = "pdf"),
          plot = weeks_ahead_forecast_plot,
          ncol = 3,
          nrow = 4)
