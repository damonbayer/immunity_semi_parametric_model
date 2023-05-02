library(tidyverse)
library(ckanr)
library(lubridate)
library(splines)
library(fs)
# Get CDPH Data -----------------------------------------------------------
time_interval_in_days <- 7

county_region_key <- read_csv(path("data", "real_data", "county_region_key.csv"))

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

ckanr_setup(url = "https://data.ca.gov")
ckan <- quiet(ckanr::src_ckan("https://data.ca.gov"))

# get resources
resources <- rbind(
  resource_search("name:covid-19", as = "table")$results,
  resource_search("name:hospitals by county", as = "table")$results
)


cases_deaths_url <- resources %>%
  filter(name == "Statewide COVID-19 Cases Deaths Tests") %>%
  pull(url)

hosp_url <- resources %>%
  filter(name == "Statewide Covid-19 Hospital County Data") %>%
  pull(url)

cases_raw <-
  read_csv(cases_deaths_url) %>%
  filter(area_type == "County") %>%
  rename(county = area) %>%
  filter(!(county %in% c("Out of state", "Unknown")))

county_id_pop <-
  cases_raw %>%
  distinct(county, population) %>%
  mutate(population = as.integer(population)) %>%
  bind_rows(
    .,
    left_join(., county_region_key, by = "county") %>%
      select(-county) %>%
      group_by(region) %>%
      summarize(
        population = sum(population),
        .groups = "drop"
      ) %>%
      rename(county = region),
    select(., -county) %>%
      summarize(population = sum(population)) %>%
      mutate(county = "California")
  ) %>%
  mutate(id = seq_len(n()), .before = 1)


top_counties <- 
  county_id_pop %>%
  head(-6) %>%
  filter(population > 1e6) %>% 
  pull(county)

cases <-
  cases_raw %>%
  select(county, date, cases, deaths, cumulative_deaths) %>%
  drop_na()

# Hospitalized should be hospitalized_covid_confirmed_patients
# ICU should be icu_covid_confirmed_patients
hosp <-
  read_csv(hosp_url) %>%
  select(county,
         date = todays_date,
         hospitalized = hospitalized_covid_confirmed_patients,
         icu = icu_covid_confirmed_patients
  )


full_cdph_dat <-
  full_join(cases, hosp) %>%
  # left_join(all_ca_seq_count_data) %>% 
  # left_join(all_ca_omicron_data) %>% 
  replace_na(list(
    hospitalized = 0,
    icu = 0
  )) %>% 
  filter(date >= global_start_date,
         date <= global_end_date) %>% 
  bind_rows(.,
            select(., -county) %>% 
              group_by(date) %>%
              summarize(across(everything(), sum), .groups = "drop") %>%
              mutate(county = "California"))

# dat <-
#   full_dat %>%
#   mutate(time = floor(as.numeric(date - min(date)) / time_interval_in_days)) %>%
#   group_by(county, time) %>%
#   summarize(
#     start_date = min(date),
#     end_date = max(date),
#     cases = sum(cases),
#     deaths = sum(deaths),
#     cumulative_deaths = last(cumulative_deaths),
#     hospitalized = last(hospitalized),
#     icu = last(icu),
#     total_seq = sum(total_seq),
#     omicron_seq = sum(omicron_seq),
#     .groups = "drop"
#   ) %>%
#   bind_rows(
#     .,
#     left_join(., county_region_key, by = "county") %>%
#       select(-county) %>%
#       group_by(region, start_date, end_date, time) %>%
#       summarize(across(everything(), sum), .groups = "drop") %>%
#       rename(county = region),
#     select(., -county) %>%
#       group_by(start_date, end_date, time) %>%
#       summarize(across(everything(), sum), .groups = "drop") %>%
#       mutate(county = "California")
#   )