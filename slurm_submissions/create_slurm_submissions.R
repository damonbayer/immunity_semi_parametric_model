library(tidyverse)
library(fs)
slurm_dir <- "slurm_submissions/simulation/compare_inference_unknown_variant_proportion_ihr_hifi"

tibble(file_path = dir_ls(slurm_dir, recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "fit_models.sh") %>% 
  pull(file_path) %>% 
  str_c("sbatch ", .) %>% 
  cat(sep = "\n")


tibble(file_path = dir_ls(slurm_dir, recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "generate_predictive_and_generated_quantities.sh") %>% 
  pull(file_path) %>% 
  str_c("sbatch ", .) %>% 
  cat(sep = "\n")


tibble(file_path = dir_ls(slurm_dir, recurse = 1, type = "file")) %>% 
  filter(path_file(file_path) == "tidy_results.sh") %>% 
  pull(file_path) %>% 
  str_c("sbatch ", .) %>% 
  cat(sep = "\n")
