args <- commandArgs(trailingOnly = TRUE)
setwd("/scratch/user/tongtong1993/Research/Analysis/singulomics")
.libPaths(c('/sw/hprc/sw/R_tamu/R_LIBS/4.1.2/GCC-11.2.0_OpenMPI-4.1.1', .libPaths()))
.libPaths(c('/scratch/user/tongtong1993/project/R/4.1', .libPaths()))

library(rlang)
library(tidyverse)
library(harmonicmeanp)

args[1] %>% as.numeric() -> se

list.files("./rda/output/Hepatocytes/Raw_data", pattern = "\\.csv", full.names = T) %>% 
#{.[grepl("(RNA|ATAC)\\.csv", .)]} %>%
{.[grepl("gene_activity\\.csv", .)]} %>%
gtools::mixedsort() -> files_
files_[se] -> raw_data

output_dir = "./rda/output/Hepatocytes/Cauchy"
output = sprintf("%s/%s", output_dir, gsub(".+/(.+)", "\\1", raw_data))
output = gsub("\\.csv", "_HMP.csv", output)

source("./Calculate_HMP.R")
res_ = cyclic_HMP(raw_data = raw_data)

write.csv(res_, file = output, row.names = F, quote = F)