library(tidyverse)
library(fs)
slurm_dir <- "slurm_submissions/real_data/ba1_immunity_contact_cdr/"

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


tibble(file_path = dir_ls(slurm_dir, type = "file")) %>% 
  pull(file_path) %>% 
  str_c("sbatch ", .) %>% 
  cat(sep = "\n")
