library(tidyverse)
library(outbreakinfo)
library(glue)
library(ckanr)
library(fs)


# Setup -------------------------------------------------------------------
first_date_to_report <- ymd("2021-05-03")
last_date_to_report <- ymd("2022-05-01")

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

# Get GISAID Data ---------------------------------------------------------
getGenomicData_known_location <- function(query_url, location=NULL, cumulative=NULL, pangolin_lineage=NULL, mutations=NULL, ndays=NULL, frequency=NULL, subadmin=NULL, other_threshold=NULL, nday_threshold=NULL, other_exclude=NULL, logInfo=TRUE){
  
  genomic_url <- "https://api.outbreak.info/genomics/"
  
  q <- c()
  
  q <- c(q, paste0(query_url), "?")
  
  if(!is.null(location)){
    # location <- getLocationIdGenomic(location)
    if(length(location) == 0){
      cat(paste0("Could not find location ", location, "\n"))
      return(NULL)
    }
    q <- c(q, paste0("location_id=", location, "&"))
  }
  if(!is.null(cumulative)){
    if (!is.logical(cumulative)){
      stop("cumulative must be in Boolean format")
    }else{
      q <- c(q, paste0("cumulative=", tolower(cumulative)), "&")
    }
  }
  if(!is.null(subadmin)){
    if (!is.logical(subadmin)){
      stop("subadmin must be in Boolean format")
    }else{
      q <- c(q, paste0("subadmin=", tolower(subadmin)), "&")
    }
  }
  if(!is.null(pangolin_lineage)){
    q <- c(q, paste0("pangolin_lineage=", pangolin_lineage, "&"))
  }
  if(!is.null(mutations)){
    check_cond <- grepl("[A-Za-z0-9]+:[A-Za-z][0-9]+[A-Za-z]", mutations)
    if(!all(check_cond))
      warning(paste0("Mutations should be specified in the format gene:mutation, like \"S:E484K\". The following mutations are not in the specified format: ",  paste(mutations[!check_cond], collapse=", ")))
    mutations <- paste(mutations, collapse=" AND ")
    q <- c(q, paste0("mutations=", mutations, "&"))
  }
  if(!is.null(ndays)){
    q <- c(q, paste0("ndays=", ndays, "&"))
  }
  if(!is.null(frequency)){
    q <- c(q, paste0("frequency=", frequency, "&"))
  }
  if(!is.null(other_threshold)){
    q <- c(q, paste0("other_threshold=", other_threshold, "&"))
  }
  if(!is.null(nday_threshold)){
    q <- c(q, paste0("nday_threshold=", nday_threshold, "&"))
  }
  if(!is.null(other_exclude)){
    other_exclude <- paste(other_exclude, collapse=",")
    q <- c(q, paste0("other_exclude=", other_exclude, "&"))
  }
  
  q <- paste(q, sep="", collapse = "")
  q <- sub("&$", "", q)
  
  dataurl <- paste0(genomic_url, q)
  results <- getGenomicsResponse(dataurl, logInfo);
  
  if (length(results) > 1){
    hits <- rbind_pages(results)
  }else{
    hits <- data.frame(results)
  }
  if ("date" %in% colnames(hits)){
    hits$date=as.Date(hits$date, "%Y-%m-%d")
    hits <- hits[order(as.Date(hits$date, format = "%Y-%m-%d")),]
  }
  return(hits)
}

county_fips_key <-
  tibble(county = c(
    "Alameda", "Kings", "Placer",
    "Sierra", "Alpine", "Lake", "Plumas", "Siskiyou", "Amador", "Lassen",
    "Riverside", "Solano", "Butte", "Los Angeles", "Sacramento",
    "Sonoma", "Calaveras", "Madera", "San Benito", "Stanislaus",
    "Colusa", "Marin", "San Bernardino", "Sutter", "Contra Costa",
    "Mariposa", "San Diego", "Tehama", "Del Norte", "Mendocino",
    "San Francisco", "Trinity", "El Dorado", "Merced", "San Joaquin",
    "Tulare", "Fresno", "Modoc", "San Luis Obispo", "Tuolumne", "Glenn",
    "Mono", "San Mateo", "Ventura", "Humboldt", "Monterey", "Santa Barbara",
    "Yolo", "Imperial", "Napa", "Santa Clara", "Yuba", "Inyo", "Nevada",
    "Santa Cruz", "Kern", "Orange", "Shasta"
  ), fips_code = c(
    1, 31,
    61, 91, 3, 33, 63, 93, 5, 35, 65, 95, 7, 37, 67, 97, 9, 39, 69,
    99, 11, 41, 71, 101, 13, 43, 73, 103, 15, 45, 75, 105, 17, 47,
    77, 107, 19, 49, 79, 109, 21, 51, 81, 111, 23, 53, 83, 113, 25,
    55, 85, 115, 27, 57, 87, 29, 59, 89
  )) %>%
  arrange(fips_code) %>%
  mutate(location_string = glue("USA_US-CA_06{sprintf('%03d', fips_code)}"))



# Multiple lineage queries should be separated by the parameter ` OR `. The `lookupSublinage` function is a convenient wrapper to look up what sublineages are associated with Omicron (e.g. B.1.617.2 and all the AY sublineages) and collapse it into a string to be used in `getPrevalence`.
omicron_lineages_string <- lookupSublineages("Omicron", returnQueryString = TRUE)

all_ca_omicron_data <-
  county_fips_key %>%
  mutate(prevalence_by_location = map(
    location_string,
    ~ as_tibble(
      getGenomicData_known_location(query_url = "prevalence-by-location", pangolin_lineage = omicron_lineages_string, location = .x)
    )
  )) %>%
  select(county, prevalence_by_location) %>%
  unnest(prevalence_by_location) %>%
  select(county, date, total_seq = total_count, omicron_seq = lineage_count)

all_ca_seq_count_data <-
  county_fips_key %>%
  mutate(sequence_count = map(
    location_string,
    ~ as_tibble(
      getGenomicData_known_location(query_url = "sequence-count", pangolin_lineage = omicron_lineages_string, location = .x)
    )
  )) %>%
  select(county, sequence_count) %>%
  unnest(sequence_count) %>% 
  rename(total_seq = total_count)



# Combine Data ------------------------------------------------------------

full_dat <-
  full_join(cases, hosp) %>%
  left_join(all_ca_seq_count_data) %>% 
  left_join(all_ca_omicron_data) %>% 
  replace_na(list(
    hospitalized = 0,
    icu = 0,
    total_seq = 0,
    omicron_seq = 0
  )) %>% 
  filter(date >= first_date_to_report,
         date <= last_date_to_report)

dat <-
  full_dat %>%
  mutate(time = floor(as.numeric(date - min(date)) / time_interval_in_days)) %>%
  group_by(county, time) %>%
  summarize(
    start_date = min(date),
    end_date = max(date),
    cases = sum(cases),
    deaths = sum(deaths),
    cumulative_deaths = last(cumulative_deaths),
    hospitalized = last(hospitalized),
    icu = last(icu),
    total_seq = sum(total_seq),
    omicron_seq = sum(omicron_seq),
    .groups = "drop"
  ) %>%
  bind_rows(
    .,
    left_join(., county_region_key, by = "county") %>%
      select(-county) %>%
      group_by(region, start_date, end_date, time) %>%
      summarize(across(everything(), sum), .groups = "drop") %>%
      rename(county = region),
    select(., -county) %>%
      group_by(start_date, end_date, time) %>%
      summarize(across(everything(), sum), .groups = "drop") %>%
      mutate(county = "California")
  )

initialization_values <-
  dat %>%
  filter(time == 0) %>%
  select(county,
         H = hospitalized,
         ICU = icu,
         D = cumulative_deaths
  ) %>%
  right_join(county_id_pop, .) %>%
  mutate(across(c(H, ICU, D), ~ . + 1)) %>%
  mutate(remaining_population = population - (H + ICU + D))


# Save Data ---------------------------------------------------------------
write_csv(dat, file = path("data", "real_data", "ca_data", ext = "csv"))
write_csv(initialization_values, file = path("data", "real_data", "initialization_values_ca_data", ext = "csv"))


# Plot Data ---------------------------------------------------------------
plot_all_data_for_region <- function(target_region) {
  popsize <- 
    county_id_pop %>% 
    filter(county == target_region) %>% 
    pull(population)
  
  dat %>% 
    select(-cumulative_deaths) %>% 
    mutate(prop_omicron = omicron_seq / total_seq) %>% 
    filter(county == target_region) %>% 
    pivot_longer(-c(county, time, ends_with("date"))) %>% 
    ggplot(aes(end_date, value)) +
    facet_wrap(~name, scales = "free_y") +
    geom_point() +
    scale_y_continuous(labels = comma) +
    ggtitle(label = glue("Omicron Wave in {target_region}"),
            subtitle = glue("Population: {comma(popsize)}")) +
    cowplot::theme_cowplot()
}


dir_create(path("figures", "data", "real_data"))

figure_tbl <- 
  county_id_pop %>%
  arrange(desc(population)) %>% 
  mutate(figure = map(county, plot_all_data_for_region)) %>% 
  mutate(file_path = path("figures", "data", "real_data", glue("all_ca_omicron_data_{county %>% str_to_lower() %>% str_replace_all(' ', '_')}_plot"), ext = "pdf"))

walk2(figure_tbl$figure, figure_tbl$file_path, ~save_plot(filename = .y, plot = .x, ncol = 3, nrow = 3))

Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))

figure_tbl$file_path %>% 
  c(path("figures", "data", "real_data", "all_ca_omicron_data", ext = "pdf")) %>% 
  str_c(collapse = " ") %>% 
  system2("pdfunite", args = .)

walk(figure_tbl$file_path, file_delete)
