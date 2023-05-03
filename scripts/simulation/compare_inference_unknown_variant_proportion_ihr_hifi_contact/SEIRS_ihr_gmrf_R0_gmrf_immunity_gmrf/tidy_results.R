target_max_t <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 28, as.integer(commandArgs(trailingOnly=T[1])))
source("src/immunity_semi_parametric_model.R")
library(fs)
context <- path("simulation", "compare_inference_unknown_variant_proportion_ihr_hifi_contact")
model_name <- "SEIRS_ihr_gmrf_R0_gmrf_immunity_gmrf"

source(path("scripts", context, "tidy_results_code.R"))
