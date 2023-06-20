library(gt)
tidy_duration_tbl %>% 
  mutate(compute_hrs = as.numeric(compute) / (60 * 60)) %>% 
  right_join(model_table) %>% 
  group_by(model_description) %>% 
  summarize(mean(compute_hrs), min(compute_hrs), max(compute_hrs)) %>% 
  xtable::xtable()



tidy_duration_tbl %>% 
  mutate(compute_hrs = as.numeric(compute) / (60 * 60)) %>% 
  right_join(model_table) %>% 
  ggplot(aes(max_t, compute_hrs, color = model_description)) +
  facet_wrap(~data_takeover_speed, scales = "free", labeller = labeller(data_takeover_speed = ~glue("{str_to_title(.x)} Takeover Data"))) +
  geom_line() +
  geom_point() +
  scale_x_continuous("Forecast Time") +
  scale_y_continuous("Compute Hours", labels = comma) +
  scale_color_discrete("Model", labels = label_parse()) +
  ggtitle("Computation Time for Simulation Study") +
  theme(legend.position = "bottom")
