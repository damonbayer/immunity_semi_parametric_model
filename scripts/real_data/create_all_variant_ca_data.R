library(tidyverse)
library(ckanr)
library(lubridate)
library(splines)
library(fs)
library(broom)
library(cowplot)
library(glue)
library(outbreakinfo)
theme_set(cowplot::theme_minimal_grid())
n_forecast_weeks <- 21
# global_start_date <- ymd("2021-11-22")
global_start_date <- ymd("2021-01-01")

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
  mutate(location_string = glue("USA_US-CA_06{sprintf('%03d', fips_code)}")) %>% 
  add_row(tibble(county = "California", fips_code = NA, location_string = "USA_US-CA"), .before = 1) %>% 
  mutate(county_id = seq_along(county), .before = 1)



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


other_threshold <- 0.01
nday_threshold <- 1
ndays <- 1000
other_exclude <- NULL
cumulative <- F

# -------------------------------------------------------------------------
ca_county_seq_data_raw <-
  county_fips_key %>%
  mutate(prevalence_by_location = map(
    location_string,
    ~ as_tibble(
      getGenomicData_known_location(query_url = "prevalence-by-location-all-lineages",
                                    location = .x, other_threshold = other_threshold,
                                    nday_threshold = nday_threshold, ndays = ndays, other_exclude = other_exclude,
                                    cumulative = cumulative)
    )
  ))

ca_county_seq_data <-
  ca_county_seq_data_raw %>% 
  select(county, prevalence_by_location) %>% 
  unnest(prevalence_by_location) %>% 
  filter(lineage_count > 0) %>% 
  filter(date >= global_start_date) %>% 
  mutate(lineage_1 = str_extract(lineage, "^([^\\.]*\\.?)"),
         lineage_2 = str_extract(lineage, "^([^\\.]*\\.?[^\\.]*\\.?)")) %>% 
  mutate(across(matches("lineage_\\d+"), ~if_else(str_ends(.x, "\\."), str_sub(.x, end = -2), .x))) %>% 
  mutate(lineage_collapsed = if_else(lineage_1 %in% c("ay", "xbb", "be", "bf", "cq"), lineage_1, lineage_2)) %>% 
  mutate(lineage_collapsed = fct_collapse(lineage_collapsed,
                                          `ba.4/ba.5/bf` = c("ba.4", "ba.5", "bf"),
                                          `bq.1/xbb.1` = c("bq.1", "xbb"))) %>%
  group_by(date, county, lineage_collapsed) %>% 
  summarize(lineage_count = sum(lineage_count),
            .groups = "drop") %>% 
  complete(date, county, lineage_collapsed, fill = list(lineage_count = 0)) %>% 
  add_count(date, county, wt = lineage_count, name = "total_count") %>% 
  mutate(other_count = total_count - lineage_count,
         prop_lineage = lineage_count / total_count)


global_end_date <- max(ca_county_seq_data$date)
source(path("scripts", "real_data", "pull_cdph_data.R"))

top_lineages <- 
  ca_county_seq_data %>% 
  filter(county == "California") %>% 
  count(lineage_collapsed, wt = lineage_count, sort = T) %>% 
  head(6) %>% 
  pull(lineage_collapsed) %>% 
  as.character()

ca_seq_augmented <- 
  ca_county_seq_data %>% 
  filter(county == "California") %>% 
  filter(lineage_collapsed %in% top_lineages) %>%
  nest(.by = lineage_collapsed) %>% 
  mutate(date_summaries = map(data,
                              ~bind_cols(.x,
                                         glm(formula = cbind(lineage_count, other_count) ~ ns(date, 8), family = "binomial", data = .x) %>%
                                           augment(type.predict = "response")))) %>% 
  select(-data) %>% 
  unnest(date_summaries)

lineage_collapsed_in_order <- 
  ca_seq_augmented %>% 
  group_by(lineage_collapsed) %>% 
  summarize(peak_date = date[.fitted == max(.fitted)]) %>% 
  arrange(peak_date) %>% 
  pull(lineage_collapsed) %>% 
  as.character()

ca_lineage_prop_plot <- 
  ca_seq_augmented %>% 
  mutate(lineage_collapsed = fct_relevel(lineage_collapsed, lineage_collapsed_in_order)) %>% 
  ggplot(aes(x = date, color = lineage_collapsed)) +
  geom_point(mapping = aes(y = prop_lineage), alpha = 0.25) +
  geom_line(mapping = aes(y = .fitted)) +
  scale_y_continuous(name = "Lineage Proportion", labels = scales::percent) +
  scale_color_discrete(name = "Lineage Group", labels = str_to_upper) +
  ggtitle("Variant Proportions in California",
          subtitle = "Smoothing via logistic regression with spline") +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(global_start_date, global_end_date))


key_dates_by_wave <-
  ca_seq_augmented %>% 
  group_by(lineage_collapsed) %>% 
  mutate(started = as.logical(cumsum(.fitted > 0.01)),
         passed_5 = as.logical(cumsum(.fitted > 0.05)),
         peaked = as.logical(cumsum(.fitted >= max(.fitted) - 0.01))) %>% 
  mutate(date_type = case_when(
    started & !passed_5 & !peaked ~ "wave_start",
    started & passed_5 & !peaked ~ "wave_offset",
    started & passed_5 & peaked ~ "wave_peak",
  )) %>% 
  drop_na() %>% 
  group_by(lineage_collapsed, date_type) %>% 
  summarize(date = min(date), .groups = "drop") %>% 
  pivot_wider(names_from = date_type, values_from = date) %>% 
  arrange(wave_peak) %>% 
  mutate(start_date = lag(wave_start), # start including cdph  data
         start_seq = wave_start - 7, # start including seq data
         start_fitting = wave_start, # max_t for first fit
         stop_fitting = wave_peak + 7 # max_t for first fit
         ) %>% 
  slice(-1) %>% 
  mutate(start_date = start_date + as.numeric(start_fitting - start_date) %% 7,
         stop_fitting = stop_fitting - as.numeric(stop_fitting - start_fitting) %% 7 - 1) %>% 
  mutate(end_date = stop_fitting + 7 * n_forecast_weeks) %>% 
  filter(end_date <= global_end_date) %>% 
  replace_na(list(start_date = global_start_date))

# Make data by variant ----------------------------------------------------
create_wave_dataset <- function(target_lineage_collapsed) {
  tmp_key_dates <- 
    key_dates_by_wave %>% 
    filter(lineage_collapsed == target_lineage_collapsed) %>% 
    select(-lineage_collapsed) %>% 
    pivot_longer(everything()) %>% 
    deframe()
  
  tmp_key_times <- set_names(as.numeric(tmp_key_dates - tmp_key_dates["start_date"]) / time_interval_in_days, names(tmp_key_dates)) %>% sort()
  
  tmp_seq_data <- 
    ca_county_seq_data %>% 
    filter(lineage_collapsed == target_lineage_collapsed,
           date >= tmp_key_dates["start_seq"],
           date <= tmp_key_dates["end_date"]) %>% 
    arrange(county, date) %>% 
    mutate(time = as.numeric(date - tmp_key_dates["start_date"]) / time_interval_in_days) %>% 
    select(-lineage_collapsed) %>% 
    left_join(select(county_id_pop, id, county))
  
  # tmp_peak_hospitalization <- 
  #   full_cdph_dat %>% 
  #   filter(date >= tmp_key_dates["start_date"],
  #          date <= tmp_key_dates["end_date"]) %>% 
  #   mutate(time = floor(as.numeric(date - tmp_key_dates["start_date"]) / time_interval_in_days)) %>% 
  #   select(county, date, hospitalized) %>% 
  #   group_by(county) %>% 
  #   filter(date > tmp_key_dates["wave_start"]) %>% 
  #   filter(hospitalized == max(hospitalized[date > tmp_key_dates["wave_start"]])) %>% 
  #   slice(1) %>% 
  #   ungroup() %>% 
  #   left_join(select(county_id_pop, id, county)) %>% 
  #   select(id, county, date, hospitalized)
  # 
  
  tmp_cdph_data <- 
    full_cdph_dat %>% 
    filter(date >= tmp_key_dates["start_date"],
           date <= tmp_key_dates["end_date"]) %>% 
    mutate(time = floor(as.numeric(date - tmp_key_dates["start_date"]) / time_interval_in_days)) %>% 
    group_by(county, time) %>% 
    summarize(
      start_date = min(date),
      end_date = max(date),
      cases = sum(cases),
      deaths = sum(deaths),
      cumulative_deaths = last(cumulative_deaths),
      hospitalized = last(hospitalized),
      icu = last(icu),
      .groups = "drop"
    ) %>% 
    left_join(select(county_id_pop, id, county))
  
  tmp_shared_constants <- 
    c("time_look_for_second_wave" = unname(tmp_key_times["wave_offset"]),
      "logistic_growth_time_offset" = unname(tmp_key_times["wave_offset"]),
      "variant_2_import_time" = unname(tmp_key_times["start_seq"]),
      "date_time_0" = tmp_cdph_data %>% filter(time == 0) %>% pull(end_date) %>% pluck(1) %>% as.character() %>% str_c('"', ., '"')) %>% 
    str_c(names(.), ., sep = " = ")
  
  list(cdph_data = tmp_cdph_data,
       seq_data = tmp_seq_data,
       # peak_hospitalization = tmp_peak_hospitalization,
       shared_constants = tmp_shared_constants)
}

all_lineage_data_tbl <- 
  key_dates_by_wave %>% 
  select(lineage_collapsed) %>% 
  mutate(all_data = map(lineage_collapsed, create_wave_dataset)) %>% 
  unnest_wider(all_data)


all_lineage_figure_tbl <- 
  all_lineage_data_tbl %>% 
  mutate(combined_data = map2(cdph_data, seq_data, full_join)) %>% 
  select(lineage_collapsed, combined_data) %>% 
  unnest(combined_data) %>% 
  mutate(date = if_else(is.na(date), end_date, date)) %>% 
  group_by(lineage_collapsed, county) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(file_path = "xxx") %>% 
  mutate(figure = pmap(list(lineage_collapsed, county, data), ~{
      select(..3, -ends_with("_date"), -time) %>% 
      pivot_longer(-c(date)) %>% 
      drop_na() %>% 
      ggplot(aes(date, value)) +
      facet_wrap(~name, scale = "free_y") +
      geom_line() +
      geom_point() +
      scale_x_date(date_breaks = "1 month", date_labels = "%m/%y") +
      ggtitle(glue("{..1} - {..2}"))
  })) %>% 
  mutate(lineage_collapsed = lineage_collapsed %>%
           str_remove_all("\\.") %>% 
           str_replace_all("/", "_"),
         county = county %>%
           str_to_lower() %>% 
           str_replace_all(" ", "_")) %>% 
  mutate(file_path = path("figures", "data", "real_data", glue("{lineage_collapsed}_{county}"), ext = "pdf")) %>% 
  select(lineage_collapsed, file_path, figure)


# Save Data ---------------------------------------------------------------
path("data", "real_data", unique(all_lineage_figure_tbl$lineage_collapsed)) %>% dir_create()

all_lineage_data_tbl %>% 
  mutate(lineage_collapsed = lineage_collapsed %>%
           str_remove_all("\\.") %>% 
           str_replace_all("/", "_")) %>% 
  mutate(initialization_values = map(
    cdph_data,
    ~{filter(.x, time == 0) %>%
        select(county,
               cases,
               H = hospitalized,
               ICU = icu,
               D = cumulative_deaths
        ) %>%
        right_join(county_id_pop, .) %>%
        mutate(across(c(H, ICU, D), ~ . + 1)) %>%
        mutate(remaining_population = population - (H + ICU + D))})) %>% 
  mutate(cdph_data = map(cdph_data, ~filter(.x, time > 0))) %>% 
  pivot_longer(-lineage_collapsed) %>% 
  mutate(data_type = map_chr(value, typeof)) %>% 
  mutate(ext = case_when( data_type == "list" ~ "csv",
                          data_type == "character" ~ "txt")) %>% 
  mutate(file_path = path("data", "real_data", lineage_collapsed, glue("{name}.{ext}"))) %>% 
  select(ext, file_path, value) %>% 
  as.list() %>% 
  pwalk(~{
    if (..1 == "csv") {
      write_csv(x = ..3, file = ..2)
    } else if (..1 == "txt") {
      write_lines(x = ..3, file = ..2)
    }
  })

# Save figures ------------------------------------------------------------
path("figures", "data", "real_data") %>% dir_create()
pwalk(as.list(all_lineage_figure_tbl), ~save_plot(filename = ..2, plot = ..3, ncol = 3, nrow = 3, device = cairo_pdf))

if (Sys.info()['sysname'] == "Darwin") {
  Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/opt/homebrew/bin/", sep=":"))
}

all_lineage_figure_tbl %>% 
  group_by(lineage_collapsed) %>% 
  summarize(file_paths = str_c(file_path, collapse = " ")) %>% 
  mutate(combined_file_path = path("figures", "data", "real_data", lineage_collapsed, ext = "pdf")) %>% 
  mutate(command_args = str_c(file_paths, combined_file_path, sep = " ")) %>%
  pull(command_args) %>% 
  walk(~system2("pdfunite", args = .x))

file_delete(all_lineage_figure_tbl$file_path)

save_plot(filename = path("figures", "data", "real_data", "ca_lineage_prop_plot", ext = "pdf"), plot = ca_lineage_prop_plot, device = cairo_pdf, base_height = 6)
