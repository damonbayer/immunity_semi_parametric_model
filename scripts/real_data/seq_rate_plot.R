
seq_data <- read_csv("/Users/damon/Documents/immunity_semi_parametric_model/data/real_data/ba1/seq_data.csv")
full_cdph_dat <- read_csv("~/full_cdph_dat.csv")

tmp <- 
  seq_data %>% 
  summarize(min_date = min(date),
            max_date = max(date))
tmp$min_date
library(splines)
library(scales)
ca_county_seq_data %>% 
  distinct(date, county, total_count) %>% 
  left_join(full_cdph_dat) %>% 
  rename(seq_count = total_count) %>% 
  mutate(seq_rate = seq_count / cases) %>% 
  mutate(not_seq_count = cases - seq_count) %>% 
  # filter(county %in% c("Orange", "California")) %>%
  filter(date >= tmp$min_date,
         date <= tmp$max_date) %>% 
  nest(.by = county) %>% 
  mutate(date_summaries = map(data,
                              ~bind_cols(.x,
                                         glm(formula = cbind(seq_count, not_seq_count) ~ ns(date, 8), family = "binomial", data = .x) %>%
                                           augment(type.predict = "response")))) %>% 
  select(-data) %>% 
  unnest(date_summaries) %>% 
  ggplot(aes(date)) +
  facet_wrap(~county) +
  geom_point(mapping = aes(y = seq_rate)) +
  geom_line(mapping = aes(y = .fitted)) +
  scale_y_continuous(labels = percent)

