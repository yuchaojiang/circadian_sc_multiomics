#res = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_deta/B_cells_seed_10_RNA.csv")

cauchy.test=function(p.vector){
  d <- length(p.vector)
  
  if(1 %in% p.vector){
    p.vector[p.vector==1]=max(0.9999, 1-1/length(p.vector))
  }
  
  Sd <- sum(tan((0.5-p.vector)*pi)/d)
  p.global <- pcauchy(Sd,lower.tail = F)
  p.global <- min(p.global,1)
  
  return(p.global)
}
#Either use csv (raw_data) as input or use exp_matrix ,gene and timepoints as input 
#The input exp matrix should be should have colnames in ZT0_REP1 format and gene name must be unique
cyclic_HMP = function(raw_data=NULL, exp_matrix=NULL, gene=NULL, timepoints=NULL, minper_ = 20, nCores_ = 5){
  library(tidyverse)
  
  # read file
  if (is_null(raw_data)){
    exp_matrix %>% as.data.frame() %>% 
      mutate(Gene = gene, .before = 1) -> df_
    #Generate tmp output files
    raw_data = sprintf("%s.csv", tempfile("HMP", tmpdir = "."))
    write.csv(df_, raw_data, row.names = FALSE, quote = FALSE)
  }else{
    print(raw_data)
    read.csv(raw_data, header = TRUE, stringsAsFactors = F) -> df_    
  }
  
  print(head(df_))
  
  #Run MetaCycle
  library(MetaCycle)
  print("Run MetaCycle")
  readLines(raw_data, n = 1) -> colnames_
  strsplit(colnames_, ",")[[1]] %>% 
    .[-1] %>% 
    gsub("ZT(\\d+)_.+", "\\1", .) %>% 
    as.numeric() -> timepoints
  print(timepoints)
  
  JTK.i=meta2d(infile = raw_data, filestyle = "csv",
               outputFile = FALSE, timepoints = timepoints, cycMethod = c("JTK"),
               maxper = 24, minper = minper_, combinePvalue = 'bonferroni', parallelize = T, nCores = nCores_)$meta
  
  #Run RAIN
#  print("Run RAIN")
#  library(rain)
#  table(timepoints) %>% unname() %>% as.integer() -> nr_
#  print(nr_)
  df_ %>% column_to_rownames("Gene") -> df_1
#  RAIN <- rain(t(df_1), period = 24, deltat = 4, peak.border = c(0.2, 0.8), method = 'independent', measure.sequence = nr_)
#  RAIN %>% rownames_to_column("CycID") -> RAIN
  
  #Run F24
#  print("Run F24")
#  source("~/Dropbox/singulomics/github/F24_BG.function.R")
#  source("~/Dropbox/singulomics/github/F24_BG.function_V2.R")
#  F24_resutls = F24(TPS = df_1, times = timepoints)
#  F24_resutls %>% as_tibble() %>% column_to_rownames("Probe") %>% rownames_to_column("CycID") -> F24_resutls
  
  #Run HR
  print("Run Harmonic Regression")
  library(HarmonicRegression)
  HARE <- harmonic.regression(inputts = t(as.matrix(df_1)), inputtime = timepoints, normalize = FALSE)
  
  HARE_results <- data.frame("CycID" = names(HARE$pvals), 
                             "p-value" = HARE$pvals, 
                             "q-value" = HARE$qvals, 
                             "phi" = HARE$pars[2], 
                             "amplitude" = HARE$pars[1])
  
  #Run BioCycle
#  print("Run BioCycle")
#  df_1 %>% 
#    "colnames<-"(., sprintf("ZT_%s_REP_%s", gsub("ZT(\\d+?)_REP(\\d+)", "\\1", colnames(.)) %>% as.numeric(), 
#                            gsub("ZT(\\d+?)_REP(\\d+)", "\\2", colnames(.)) %>% as.numeric())) %>% 
#    rownames_to_column("ID") -> Biocycle_input
#  
#  tmp_dir = tempfile(pattern = "BIOCYLCE", tmpdir = ".")
#  dir.create(tmp_dir)
#  write.table(Biocycle_input, file = sprintf("%s/Biocycle_input.tsv", tmp_dir), sep = "\t", quote = F, row.names = F)
#  
#  Biocycle_ = "~/Dropbox/singulomics/github/BioCycle_0.9.3/BioCycle.R"
#  input_file = sprintf("%s/Biocycle_input.tsv", tmp_dir)
#  command = sprintf("Rscript %s -i %s -o %s -s 18 -e 28", Biocycle_, input_file, tmp_dir)
#  system(command)
#  
#  read.table(sprintf("%s/results.tsv", tmp_dir), head = T, sep = "\t", stringsAsFactors = F) -> Biocycle_results
#  Biocycle_results %>% column_to_rownames("ID") %>% 
#    .[,1:7] %>%
#    "colnames<-"(., sprintf("BioCycle_%s", colnames(.))) %>% 
#    rownames_to_column("CycID") -> Biocycle_results
#  
#  unlink(tmp_dir, recursive = T)
  
  ####
#  list_ = list(MetaCycle = JTK.i, RAIN = RAIN, F24 = F24_resutls, HR = HARE_results)
#  list_ = list(MetaCycle = JTK.i, Biocycle = Biocycle_results, F24 = F24_resutls, HR = HARE_results)
#  list_ = list(MetaCycle = JTK.i, Biocycle = Biocycle_results, HR = HARE_results)
#  list_ = list(MetaCycle = JTK.i, F24 = F24_resutls, HR = HARE_results)
  list_ = list(MetaCycle = JTK.i, HR = HARE_results)
  
  print("Run HMP")
  library(harmonicmeanp)
  list_ %>% 
    map2(.x=.,.y=names(.),.f=function(x,y){
      method_ = y
      df_ = x
      df_ %>% as_tibble() %>% column_to_rownames("CycID") %>% 
        "colnames<-"(., sprintf("%s_%s", method_, colnames(.))) -> df_
      df_ %>% rownames_to_column("Gene") -> df_
    }) %>% purrr::reduce(., full_join, by = "Gene") -> df_
  
  df_ %>% column_to_rownames("Gene") %>% 
    dplyr::select(matches("pval|p\\.value|pvalue|pVal|P_VALUE")) %>%
    dplyr::select(-matches("meta2d")) %>%
    apply(., 1, function(x){hmp.stat(p = x, w = NULL)}) %>% 
    as.data.frame() %>% "colnames<-"(., "HMP") %>% 
    rownames_to_column("Gene") %>% 
    mutate(HMP_BH.Q = p.adjust(HMP, method = "BH")) %>%
    mutate(HMP_Bonferroni = p.adjust(HMP, method = "bonferroni")) -> df_HMP
  
  print("Run Cauchy")
  df_ %>% column_to_rownames("Gene") %>% 
    dplyr::select(matches("pval|p\\.value|pvalue|pVal|P_VALUE")) %>%
    dplyr::select(-matches("meta2d")) %>%
    apply(., 1, function(x){cauchy.test(x)}) %>% 
    as.data.frame() %>% "colnames<-"(., "cauchy_p") %>% 
    rownames_to_column("Gene") %>% 
    mutate(cauchy_BH.Q = p.adjust(cauchy_p, method = "BH")) %>%
    mutate(cauchy_Bonferroni = p.adjust(cauchy_p, method = "bonferroni")) -> df_cauchy
  
  left_join(x = df_HMP, y = df_, by = "Gene") -> df_
  left_join(x = df_cauchy, y = df_, by = "Gene") -> df_
  
  if (!is_null(exp_matrix)){
    file.remove(raw_data)
    print(sprintf("Removed: %s", raw_data))    
  }
  print("Done")
  return(df_)
}

recal_cauchy_p_and_hmp = function(df_){
  df_ %>% dplyr::select(-matches("HMP|cauchy")) -> df_
  
  print("Run HMP")
  library(harmonicmeanp)
  
  df_ %>% column_to_rownames("Gene") %>% 
    dplyr::select(matches("pval|p\\.value|pvalue|pVal|P_VALUE")) %>%
    dplyr::select(-matches("meta2d")) %>%
    apply(., 1, function(x){hmp.stat(p = x, w = NULL)}) %>% 
    as.data.frame() %>% "colnames<-"(., "HMP") %>% 
    rownames_to_column("Gene") %>% 
    mutate(HMP_BH.Q = p.adjust(HMP, method = "BH")) %>%
    mutate(HMP_Bonferroni = p.adjust(HMP, method = "bonferroni")) -> df_HMP
  
  print("Run Cauchy")
  df_ %>% column_to_rownames("Gene") %>% 
    dplyr::select(matches("pval|p\\.value|pvalue|pVal|P_VALUE")) %>%
    dplyr::select(-matches("meta2d")) %>%
    apply(., 1, function(x){cauchy.test(x)}) %>% 
    as.data.frame() %>% "colnames<-"(., "cauchy_p") %>% 
    rownames_to_column("Gene") %>% 
    mutate(cauchy_BH.Q = p.adjust(cauchy_p, method = "BH")) %>%
    mutate(cauchy_Bonferroni = p.adjust(cauchy_p, method = "bonferroni")) -> df_cauchy
  
  left_join(x = df_HMP, y = df_, by = "Gene") -> df_
  left_join(x = df_cauchy, y = df_, by = "Gene") -> df_
  
  return(df_)
}
