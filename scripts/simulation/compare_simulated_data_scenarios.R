source("src/immunity_semi_parametric_model.R")

library(fs)


all_dat <- 
  tibble(file_path =dir_ls(path("data", "simulation"), recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "simulated_data.csv") %>% 
  mutate(simulation = map_chr(path_split(file_path), ~pluck(.x, 3))) %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_predictive())) %>% 
  select(-file_path) %>% 
  unnest(dat)

all_generated_quantities <- 
  tibble(file_path =dir_ls(path("data", "simulation"), recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "true_generated_quantities.csv") %>% 
  mutate(simulation = map_chr(path_split(file_path), ~pluck(.x, 3))) %>% 
  mutate(dat = map(file_path, ~read_csv_as_draws(.x) %>% tidy_generated_quantities())) %>% 
  select(-file_path) %>% 
  unnest(dat)

all_dat %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, ymin = .lower, ymax = .upper, fill = simulation, group = .width)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5)

all_generated_quantities %>% 
  filter(!is.na(time)) %>% 
  distinct(simulation, time, name, value) %>% 
  mutate(simulation = factor(simulation)) %>% 
  ggplot(aes(time, value, color = simulation)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line(alpha = 0.5)

unique(all_generated_quantities$simulation)
all_generated_quantities %>% 
  filter(simulation == "compare_inference_unknown_variant_proportion_ihr_contact",
         name == "prop_variant_2") %>% 
  distinct(time, name, value) %>% 
  arrange(abs(value - 0.99))
