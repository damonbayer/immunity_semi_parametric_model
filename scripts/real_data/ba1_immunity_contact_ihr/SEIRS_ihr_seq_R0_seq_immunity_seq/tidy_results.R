target_max_t <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 30, as.integer(commandArgs(trailingOnly=T[1])))
source("src/immunity_semi_parametric_model.R")
library(fs)
context <- path("real_data", "ba1_immunity_contact_ihr")
model_name <- "SEIRS_ihr_seq_R0_seq_immunity_seq"
data_dir <- path("data/real_data/ba1")
source(path("scripts", context, "tidy_results_code.R"))
