library(tidyverse)
source("src/immunity_semi_parametric_model.R")

experiment_name <- "ba1_immunity_contact_cdr"
context <- path("real_data", experiment_name)
results_dir <- path("results", context)
data_dir <- path("data", "real_data", "ba1")
simulation_dir <- path("scripts", context)

source(path(data_dir, "shared_constants.txt"))

date_time_0 <- ymd(date_time_0)
time_to_date <- function(time) date_time_0 + time  * 7

all_dat <- 
  full_join(
    read_csv(path(data_dir, "cdph_data.csv")),
    read_csv(path(data_dir, "seq_data.csv"))) %>% 
  mutate(date = if_else(is.na(date), end_date, date))
                     
                     
model_table <- read_csv(path(simulation_dir, "model_table.csv"))

figure_tbl <- 
  model_table %>% 
  group_by(county) %>% 
  summarize(min_max_t = min(max_t),
            max_max_t = max(max_t)) %>% 
  mutate(county_name = str_c(county, if_else(county == "California", "", " County")),
         first_forecast_target_date = time_to_date(min_max_t + 1),
         last_forecast_target_date = time_to_date(max_max_t + 4)) %>% 
  mutate(figure = pmap(list(county,
                         first_forecast_target_date,
                         last_forecast_target_date,
                         county_name),
                    ~{
                      all_dat %>% 
                        filter(county == ..1,
                               date < ..3) %>%
                        select(
                          date,
                          `New Cases` = cases,
                          `Hospital Occupancy` = hospitalized,
                          `ICU Occupancy` = icu,
                          `New Deaths` = deaths,
                          `BA.1 Sequences` = lineage_count,
                          `All Sequences` = total_count,
                          `Prop. BA.1 Sequences` = prop_lineage) %>% 
                        pivot_longer(-date) %>% 
                        drop_na() %>% 
                        mutate(name = fct_inorder(name)) %>% 
                        ggplot(aes(date, value)) +
                        facet_wrap(~name, scale = "free_y") +
                        annotate("rect", xmin = ..2, xmax = ..3, ymin = -Inf, ymax = Inf, alpha = 0.5) +
                        geom_line() +
                        geom_point() +
                        scale_y_continuous(name = "Value", labels = comma) +
                        scale_x_date(name = "Date", date_labels = "%b '%y", date_breaks = "2 months", guide = guide_axis(angle = 90)) +
                        ggtitle(glue("{..4} Omicron BA.1 Wave"))
  })) %>% 
  mutate(figure_dims = future_map(figure, gg_facet_dims)) %>%
  unnest_wider(figure_dims) %>%
  rename(n_row = ROW, n_col = COL) %>% 
  mutate(file_path = path(manuscript_figure_dir, glue("{county_name %>% str_to_lower() %>% str_replace_all(' ', '_')}_binned_data_plot"), ext = "pdf")) %>% 
  select(figure, n_row, n_col, file_path)

pwalk(as.list(figure_tbl), ~save_plot(filename = ..4, plot = ..1, ncol = ..3, nrow = ..2, base_asp = 1, base_height = 3))