library(tidyverse)
source("src/immunity_semi_parametric_model.R")
library(fs)
context <- "ba1_immunity_contact_cdr"
data_dir <- path("data", "real_data", "ba1")
source(path(data_dir, "shared_constants", ext = "txt"))
initialization_values <- path(data_dir, "initialization_values", ext = "csv") %>% 
  read_csv()


model_table <- 
  initialization_values %>% 
  distinct(county_id = id, county) %>% 
  filter(county %in% c("Orange", "California")) %>% 
  expand_grid(max_t = 30:36,
              CDR_model = "gmrf") %>% 
  expand_grid(tribble(
    ~`R₀_model`, ~immunity_model,
    "seq-informed",  "seq-informed",
    "constant",  "seq-informed",
    "gmrf",      "constant",
    "constant",          "gmrf"
  )) %>% 
  mutate(fit_id = 1:n(), .before = 1)

write_csv(model_table, path("scripts", "real_data", context, "model_table", ext = "csv"))
