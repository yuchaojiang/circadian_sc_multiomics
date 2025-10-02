setwd("~/Dropbox/singulomics/")

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # mm10
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(biovizBase)
library(SingleCellExperiment)
library(mbkmeans)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(MetaCycle)
library(future)
library(furrr)

#source('github/runJTK.R')


# 1. Load in the integrated and QC'ed data ----
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells

sc <- sc[,hepatocytes_cells]
resolution=0.5
sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=resolution)
sc=sc[,!grepl('KO',sc$group)] # Remove the knockout to look at circadian rhythmicity
sc$group=droplevels(sc$group)

sc$ZT=as.numeric(gsub('ZT','',sc$ZT)) # Change ZT to be numeric time
sc=sc[,sc$celltype=='Hepatocytes'] # Keep only hepatocytes

sc$cc_cluster=NULL
sc$cc_clusters=sc$seurat_clusters

# Minimum number of cells per timepoint in each cluster: 50
sc=sc[,sc$cc_clusters %in% names(which(apply(table(sc$cc_clusters, sc$group), 1, min)>50))]
sc$cc_clusters=droplevels(sc$cc_clusters)
sc$cluster=sc$cc_clusters

genes=rownames(sc@assays$RNA)
cyc.genes=c('Arntl','Per1','Per2','Npas2','Dbp','Cry1','Clock','Nr1d1','Rorc')

# 2. Generate normalized RNA matrix (gene x cells)
sc@assays$RNA@counts %>% 
  as.data.frame() %>%
  mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) -> RNA_norm

# 3. Generate Mean RNA matrix (gene x metacells) ----
sc@meta.data$cluster %>% 
  as.character() %>%
  as.numeric() %>%
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
#  .[1:2] %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(cluster==cluster_) -> df_
    seq(2,22,4) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT==ZT_) -> df_
        print(nrow(df_))
        rownames(df_) -> cellnames
        RNA_norm[ ,cellnames] -> df_
        rowMeans(df_) %>% 
          as.data.frame() %>% 
          "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, cluster_+1)) -> df_
        df_ %>% rownames_to_column("Gene") -> df_
      }) %>% purrr::reduce(., left_join, by="Gene") -> df_
  }) %>% purrr::reduce(., full_join, by="Gene") %>% 
  {.[, c("Gene", colnames(.)[-1] %>% gtools::mixedsort())]} -> RNA_Mean

# 4. Generate RNA non-zero mean matrix (gene x metacells) ----
sc@meta.data$cluster %>% 
  as.character() %>%
  as.numeric() %>%
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
  #  .[1:2] %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(cluster==cluster_) -> df_
    seq(2,22,4) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT==ZT_) -> df_
        print(nrow(df_))
        rownames(df_) -> cellnames
        RNA_norm[ ,cellnames] -> df_
        apply(df_, 1, function(x){
          x[x!=0] -> x
          mean(x) -> x
          if (is.na(x)){
            x = 0
          }else{
            x = x
          }
          return(x)
          }) %>% 
#        rowMeans(df_) %>% 
          as.data.frame() %>% 
          "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, cluster_+1)) -> df_
        df_ %>% rownames_to_column("Gene") -> df_
      }) %>% purrr::reduce(., left_join, by="Gene") -> df_
  }) %>% purrr::reduce(., full_join, by="Gene") %>% 
  {.[, c("Gene", colnames(.)[-1] %>% gtools::mixedsort())]} -> RNA_Nonzero_Mean

# 5. Generate RNA non-zero proportion matrix (gene x metacells) ----
sc@meta.data$cluster %>% 
  as.character() %>%
  as.numeric() %>%
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
#    .[1:2] %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(cluster==cluster_) -> df_
    seq(2,22,4) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT==ZT_) -> df_
        print(nrow(df_))
        rownames(df_) -> cellnames
        RNA_norm[ ,cellnames] -> df_
        df_ %>% 
          apply(., 1, function(x){
            zero_cell = sum(x==0)
            non_zero_cell = sum(x!=0)
            nono_zero_prop = non_zero_cell/(zero_cell+non_zero_cell)
          }) %>% as.data.frame() %>% 
          "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, cluster_+1)) %>% 
          rownames_to_column("Gene") -> df_
      }) %>% purrr::reduce(., left_join, by="Gene") -> df_
  }) %>% purrr::reduce(., full_join, by="Gene") %>% 
  {.[, c("Gene", colnames(.)[-1] %>% gtools::mixedsort())]} -> RNA_NonZeroProp

# 6. Generate RNA gini matrix (gene x metacells) ----
sc@meta.data$cluster %>% 
  as.character() %>%
  as.numeric() %>%
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
  #    .[1:2] %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(cluster==cluster_) -> df_
    seq(2,22,4) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT==ZT_) -> df_
        print(nrow(df_))
        rownames(df_) -> cellnames
        RNA_norm[ ,cellnames] -> df_
        df_ %>% 
          apply(., 1, function(x){
            DescTools::Gini(x) -> x
            if (is.na(x)){
              x = 1
            }else{
              x = x
            }
            return(x)
          }) %>% as.data.frame() %>% 
          "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, cluster_+1)) %>% 
          rownames_to_column("Gene") -> df_
      }) %>% purrr::reduce(., left_join, by="Gene") -> df_
  }) %>% purrr::reduce(., full_join, by="Gene") %>% 
  {.[, c("Gene", colnames(.)[-1] %>% gtools::mixedsort())]} -> RNA_gini

#7. Generate RNA cv matrix (gene x metacells) ----
sc@meta.data$cluster %>% 
  as.character() %>%
  as.numeric() %>%
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
  #    .[1:2] %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(cluster==cluster_) -> df_
    seq(2,22,4) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT==ZT_) -> df_
        print(nrow(df_))
        rownames(df_) -> cellnames
        RNA_norm[ ,cellnames] -> df_
        df_ %>% 
          apply(., 1, function(x){
            x = sd(x)/mean(x)
            if (is.na(x)){
              x = 0 
            }else{
              x = x
            }
            return(x)
          }) %>% as.data.frame() %>% 
          "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, cluster_+1)) %>% 
          rownames_to_column("Gene") -> df_
      }) %>% purrr::reduce(., left_join, by="Gene") -> df_
  }) %>% purrr::reduce(., full_join, by="Gene") %>% 
  {.[, c("Gene", colnames(.)[-1] %>% gtools::mixedsort())]} -> RNA_cv

#Save matrix
write.csv(RNA_Mean, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_Mean.csv", row.names = F, quote = F)
write.csv(RNA_NonZeroProp, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_NonZeroProp.csv", row.names = F, quote = F)
write.csv(RNA_Nonzero_Mean, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_Nonzero_Mean.csv", row.names = F, quote = F)
write.csv(RNA_gini, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_gini.csv", row.names = F, quote = F)
write.csv(RNA_cv, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_cv.csv", row.names = F, quote = F)
####

####
radian_to_phase = function(radian){
  phase = (radian/(2*pi))*24
  return(phase)
}
####

#8. Calculate rhtyhmicity ----
source("~/Dropbox/singulomics/github/Calculate_HMP.R")
res_RNA_Mean = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_Mean.csv", minper_ = 24)
res_RNA_Mean %>% 
  dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
  dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase) %>% 
  recal_cauchy_p_and_hmp(.) -> res_RNA_Mean
write.csv(res_RNA_Mean, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Mean.csv", row.names = F, quote = F, col.names = T)

res_RNA_NonZeroProp = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_NonZeroProp.csv", minper_ = 24)
res_RNA_NonZeroProp %>% 
  dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
  dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase) %>% 
  recal_cauchy_p_and_hmp(.) -> res_RNA_NonZeroProp
write.csv(res_RNA_NonZeroProp, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_NonZeroProp.csv", row.names = F, quote = F, col.names = T)

res_RNA_Nonzero_Mean = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_Nonzero_Mean.csv", minper_ = 24)
res_RNA_Nonzero_Mean %>% 
  dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
  dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase) %>% 
  recal_cauchy_p_and_hmp(.) -> res_RNA_Nonzero_Mean
write.csv(res_RNA_Nonzero_Mean, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Nonzero_Mean.csv", row.names = F, quote = F, col.names = T)

res_RNA_gini = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_gini.csv", minper_ = 24)
res_RNA_gini %>% 
  dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
  dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase) %>% 
  recal_cauchy_p_and_hmp(.) -> res_RNA_gini
write.csv(res_RNA_gini, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_gini.csv", row.names = F, quote = F, col.names = T)

res_RNA_cv = cyclic_HMP(raw_data = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_cv.csv", minper_ = 24)
res_RNA_cv %>% 
  dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
  dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase) %>% 
  recal_cauchy_p_and_hmp(.) -> res_RNA_cv
write.csv(res_RNA_cv, "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_cv.csv", row.names = F, quote = F, col.names = T)

#9. Plot RNA mean and RNA non zero proportion of clock genes ----
se_genes = c("Arntl", "Per1", "Per2", "Npas2", "Dbp", "Cry1", "Clock", "Nr1d1", "Rorc")
se_genes %>% 
  map(function(x){
    res_RNA_Mean %>% dplyr::filter(Gene == x)
  }) %>% do.call(rbind, .) %>% .[,1:4]
se_genes %>% 
  map(function(x){
    res_RNA_NonZeroProp %>% dplyr::filter(Gene == x)
  }) %>% do.call(rbind, .) %>% .[,1:4]

list_df = list()
se_genes %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    RNA_Mean %>% dplyr::filter(Gene == gene_) %>% 
      pivot_longer(cols = -Gene, names_to = "ZT", values_to = "Mean") %>% 
      mutate(REP = gsub(".+_(REP\\d+)", "\\1", ZT)) %>% 
      mutate(ZT = gsub("ZT(\\d+)_REP\\d+", "\\1", ZT) %>% as.numeric()) -> df_
    
    df_ %>% mutate(linetype_ = case_when(
      REP == "REP8" ~ "B_dashed",
      TRUE ~ "A_solid"
    )) -> df_
    list_df[[gene_]] <<- df_
    
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean, group = REP, linetype = linetype_)) + 
      geom_line(color = "#56B4E9") + 
      geom_point(color = "#56B4E9") -> p
    
#    res_RNA_Mean %>% filter(Gene == gene_) %>% .$cauchy_Bonferroni %>% sprintf("%.2e", .) -> p_val
    res_RNA_Mean %>% dplyr::filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .) -> p_val
#    p + ggtitle(sprintf("%s\nBonferroni=%s", gene_, p_val)) -> p
    p + ggtitle(sprintf("%s\np.adj=%s", gene_, p_val)) -> p
    p + theme(legend.position = "none") -> p
    p + scale_x_continuous(breaks = seq(2, 22, 4)) -> p
  }) -> p_list_Mean

patchwork::wrap_plots(p_list_Mean, ncol = 3) -> p_mean

list_df_1 = list()
se_genes %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    RNA_NonZeroProp %>% dplyr::filter(Gene == gene_) %>% 
      pivot_longer(cols = -Gene, names_to = "ZT", values_to = "Nonzero_Prop") %>% 
      mutate(REP = gsub(".+_(REP\\d+)", "\\1", ZT)) %>% 
      mutate(ZT = gsub("ZT(\\d+)_REP\\d+", "\\1", ZT) %>% as.numeric()) -> df_
    
    df_ %>% mutate(linetype_ = case_when(
      REP == "REP8" ~ "B_dashed",
      TRUE ~ "A_solid"
    )) -> df_
    list_df_1[[gene_]] <<- df_
    
    df_ %>% 
      ggplot(aes(x = ZT, y = Nonzero_Prop, group = REP, linetype = linetype_)) + 
      geom_line(color = "#009E73") + 
      geom_point(color = "#009E73") -> p
    
#    res_RNA_NonZeroProp %>% filter(Gene == gene_) %>% .$cauchy_Bonferroni %>% sprintf("%.2e", .) -> p_val
    res_RNA_NonZeroProp %>% dplyr::filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .) -> p_val
#    p + ggtitle(sprintf("%s\nBonferroni=%s", gene_, p_val)) -> p
    p + ggtitle(sprintf("%s\np.adj=%s", gene_, p_val)) -> p
    p + theme(legend.position = "none") -> p
    p + scale_x_continuous(breaks = seq(2, 22, 4)) -> p
  }) -> p_list_NonZeroProp

patchwork::wrap_plots(p_list_NonZeroProp, ncol = 3) -> p_non_zero

se_genes %>% 
  "names<-"(., .) %>% 
#  .[1] %>% 
  map(function(x){
    gene_ = x
    list_df[[gene_]] -> df_
    list_df_1[[gene_]] -> df_1
    
    df_ %>% group_by(Gene, ZT) %>% 
      group_map(function(x,y){
        df_ = x
        data.frame(Gene = y$Gene, ZT = y$ZT, Mean = mean(df_$Mean), Mean_SD = sd(df_$Mean))
      }, .keep = T) %>% 
      do.call(rbind, .) -> df_
    
    df_1 %>% group_by(Gene, ZT) %>% 
      group_map(function(x,y){
        df_ = x
        data.frame(Gene = y$Gene, ZT = y$ZT, Nonzero_Prop = mean(df_$Nonzero_Prop), Nonzero_SD = sd(df_$Nonzero_Prop))
      }, .keep = T) %>% 
      do.call(rbind, .) -> df_1
    
    left_join(df_, df_1, by = c("Gene", "ZT")) -> df_
    
    scaleFactor <- max(df_$Mean) / max(df_$Nonzero_Prop)
    res_RNA_Mean %>% dplyr::filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .) -> p_val_mean
    res_RNA_NonZeroProp %>% dplyr::filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .) -> p_val_non_zero
#    
    ggplot(df_, aes(x = ZT)) + 
      geom_line(aes(y = Mean), color = "#56B4E9") + 
      geom_ribbon(aes(ymin = Mean - Mean_SD, ymax = Mean + Mean_SD), fill = "#56B4E9", alpha = 0.3) + 
#      geom_point(aes(y = Mean), color = "#56B4E9") + 
      geom_line(aes(y = Nonzero_Prop*scaleFactor), color = "#009E73") + 
      geom_ribbon(aes(ymin = Nonzero_Prop*scaleFactor - Nonzero_SD*scaleFactor, ymax = Nonzero_Prop*scaleFactor + Nonzero_SD*scaleFactor), fill = "#009E73", alpha = 0.3) + 
#      geom_point(aes(y = Nonzero_Prop*scaleFactor), color = "#009E73") +
      scale_y_continuous(name="Mean", sec.axis=sec_axis(~./scaleFactor, name="Nonzero_Prop")) + 
      scale_x_continuous(breaks = seq(2, 22, 4)) +
      theme(legend.position = "none") + 
##      ggtitle(gene_) -> p
      ggtitle(sprintf("%s\n(Mean) p.adj=%s\n(Nonzero prop) p.adj=%s", gene_, p_val_mean, p_val_non_zero)) + 
      theme(
        plot.title = element_text(size = 10)  # Adjust the font size here
      ) -> p
  }) -> p_list
wrap_plots(p_list, ncol = 3) -> p #Fig_3B

####

# 10. Plot venndiagram (RNA Mean vs RNA Non zero prop) ----
read.csv(file = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_Mean.csv", header = T) -> RNA_mean
read.csv(file = "~/Dropbox/singulomics/github_rda/output/Beyond_Mean/RNA_NonZeroProp.csv", header = T) -> Nonzero_prop_mean

source("~/Dropbox/singulomics/github/Calculate_HMP.R")
radian_to_phase = function(radian){
  phase = (radian/(2*pi))*24
  return(phase)
}

list_df = 
  list(
    Mean = RNA_mean, 
    Nonzero_prop = Nonzero_prop_mean
  )

list_df %>% 
  purrr::map(function(df_){
    df_ %>% pivot_longer(-Gene, names_to = "ZT", values_to = "Expression") %>% dplyr::mutate(ZT = gsub("ZT(\\d+)_.+", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% group_by(Gene, ZT) %>% summarise(Mean_Expression = mean(Expression)) %>% ungroup() %>% as.data.frame() -> df_mean
    df_mean %>% group_by(Gene) %>% summarise(Max = max(Mean_Expression), Min = min(Mean_Expression)) %>% ungroup() %>% as.data.frame() -> df_range
    df_range %>% dplyr::mutate(fc = Max/Min) %>% dplyr::mutate(log2_fc = log(fc, 2)) -> df_final
  }) -> list_df_fc

#Some fc show Nan or Inf ,Add pseudo count (non zero min) to and re-calculate fc and log2_fc
list_df_fc %>% 
  purrr::map(function(df_){
    df_$Min %>% .[.!=0] %>% min() -> non_zero_min
    df_[["fc"]] = (df_[["Max"]]+non_zero_min)/(df_[["Min"]]+non_zero_min)
    df_[["log2_fc"]] = log(df_[["fc"]], 2)
    return(df_)
  }) -> list_df_fc

list_df %>% 
  purrr::map(function(df_){
    write.csv(df_, "~/Downloads/temp_mean.csv", row.names = F, quote = F, col.names = T)  
    res_Mean = cyclic_HMP(raw_data = "~/Downloads/temp_mean.csv", minper_ = 24)
    res_Mean %>% 
      dplyr::mutate(F24_Phase = radian_to_phase(HR_phi)) %>% 
      dplyr::select(Gene, MetaCycle_JTK_pvalue, MetaCycle_JTK_BH.Q, HR_p.value, HR_q.value, MetaCycle_JTK_amplitude, F24_Phase, MetaCycle_meta2d_Base, MetaCycle_meta2d_AMP, MetaCycle_meta2d_rAMP) %>% 
      recal_cauchy_p_and_hmp(.) -> res_Mean
  }) -> list_res_df

list_res_df %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, assay_){
    list_df_fc[[assay_]] -> df_fc
    left_join(df_, df_fc, by="Gene")
  }) -> res_list

list(
  Mean = res_list$Mean %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01, !is.na(MetaCycle_meta2d_AMP), fc > 1.6) %>% .$Gene,
  Nonzero_prop = res_list$Nonzero_prop %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01, !is.na(MetaCycle_meta2d_AMP), fc > 1.6) %>% .$Gene #5322 
) -> genes_list
ggvenn::ggvenn(genes_list) -> p_ggvenn_new #Fig_3C
####

# 11. Plot boxplot (Mean RNA expression and non zero proportion) ----
list(
Overlapped = intersect(genes_list$Mean, genes_list$Nonzero_Prop),
Mean_specific = setdiff(genes_list$Mean, genes_list$Nonzero_Prop),
Nonzero_specific = setdiff(genes_list$Nonzero_Prop, genes_list$Mean)
) -> genes_list

c("p_val", "log10_p_val") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    type_ = x
    if (type_ == "log10_p_val"){
      genes_list %>% 
        map2(.x=.,.y=names(.),.f=function(x,y){
          gene_ = x
          group_ = y %>% gsub("_specific", "", .)
          RNA_Mean %>% column_to_rownames("Gene") %>% .[gene_,] %>% as.vector() -> mean_
          purrr::reduce(mean_, c) -> mean_
          mean_ %>% 
            #      {.+1} %>% 
            log(., 10) %>% as.data.frame() %>% 
            "colnames<-"(., "expression") %>% mutate(group = group_) -> df_
        }) %>% do.call(rbind, .) -> df_
      df_ %>% dplyr::filter(is.finite(expression)) -> df_ 
    }else{
      genes_list %>% 
        map2(.x=.,.y=names(.),.f=function(x,y){
          gene_ = x
          group_ = y %>% gsub("_specific", "", .)
          RNA_Mean %>% column_to_rownames("Gene") %>% .[gene_,] %>% as.vector() -> mean_
          purrr::reduce(mean_, c) -> mean_
          mean_ %>% 
            #      {.+1} %>% 
            as.data.frame() %>% 
            "colnames<-"(., "expression") %>% mutate(group = group_) -> df_
        }) %>% do.call(rbind, .) -> df_
      df_ %>% dplyr::filter(is.finite(expression)) -> df_ 
    }
    return(df_)
  }) -> list_tmp

list_tmp$log10_p_val %>% 
  mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
  ggplot(aes(x = group, y = expression, fill = group)) + 
  geom_boxplot(width = 0.5, outlier.shape = NA) + 
  theme_classic() + 
  ggpubr::stat_compare_means(comparisons = list(c("Mean", "Overlapped"), 
                                                c("Overlapped", "Nonzero"), 
                                                c("Mean", "Nonzero")), method = "t.test") + 
  theme(legend.position = "none") + 
  xlab(NULL) + 
  ylab("log10(exprssion)") -> p_1

library(ggbreak)
list_tmp$p_val %>% 
  mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
  ggplot(aes(x = group, y = expression, fill = group)) + 
  geom_boxplot(width = 0.5) + 
  scale_y_break(c(1, 50), scales = 0.5) + 
  theme_classic() + 
  ggpubr::stat_compare_means(comparisons = list(c("Mean", "Overlapped"), 
                                                c("Overlapped", "Nonzero"), 
                                                c("Mean", "Nonzero")), method = "t.test") + 
  theme(legend.position = "none") + 
  xlab(NULL) + 
  ylab("exprssion") -> p_1_1

genes_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = x
    group_ = y %>% gsub("_specific", "", .)
    RNA_NonZeroProp %>% column_to_rownames("Gene") %>% .[gene_,] %>% as.vector() -> mean_
    purrr::reduce(mean_, c) -> mean_
    mean_ %>% 
      as.data.frame() %>% 
      "colnames<-"(., "expression") %>% mutate(group = group_) -> df_
  }) %>% do.call(rbind, .) -> df_
df_ %>% dplyr::filter(is.finite(expression)) -> df_

df_ %>% 
  mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
  ggplot(aes(x = group, y = expression, fill = group)) + 
  geom_boxplot(width = 0.5) + 
  theme_classic() + 
  ggpubr::stat_compare_means(comparisons = list(c("Mean", "Overlapped"), 
                                                c("Overlapped", "Nonzero"), 
                                                c("Mean", "Nonzero")), method = "t.test") + 
  theme(legend.position = "none") + 
  xlab(NULL) + 
  ylab("Non zero proportion") -> p_2

patchwork::wrap_plots(p_2, p_1, ncol = 2)
patchwork::wrap_plots(p_2, p_1_1, ncol = 2) #Fig_3D

read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Nonzero_Mean.csv", header = T, stringsAsFactors = F) -> res_RNA_Nonzero_Mean
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_gini.csv", header = T, stringsAsFactors = F) -> res_RNA_gini
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_cv.csv", header = T, stringsAsFactors = F) -> res_RNA_cv

library(ggupset)
list(
  Mean = res_list$Mean %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01, !is.na(MetaCycle_meta2d_AMP), fc > 1.6) %>% .$Gene, 
  Nonzero_Prop = res_list$Nonzero_prop %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01, !is.na(MetaCycle_meta2d_AMP), fc > 1.6) %>% .$Gene,
  Nonzero_mean = res_RNA_Nonzero_Mean %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene,
  cv = res_RNA_cv %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.05) %>% .$Gene,
  gini = res_RNA_gini %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>% 
  as.data.frame() %>% 
  UpSetR::upset(sets = c("Mean", "cv", "Nonzero_Prop", "gini", "Nonzero_mean"), keep.order = T, order.by = "freq") -> p_upset #Fig_S16A

# 12. Plot correlation between pval, phase, and amp between RNA mean and RNA non zero proportion ----
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Mean.csv", header = T, stringsAsFactors = F) -> res_RNA_Mean
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_NonZeroProp.csv", header = T, stringsAsFactors = F) -> res_RNA_NonZeroProp

c("cauchy_p", "cauchy_BH.Q") %>% 
  .[2] %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    method_ = x
    list(Mean = res_RNA_Mean, Nonzero_Prop = res_RNA_NonZeroProp) %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        type_ = y
        df_ = x
        -log(x[[method_]], 10) -> p_val_
        data.frame(Gene = df_$Gene, p_val = p_val_) %>% 
          "colnames<-"(., c("Gene", type_))
      }) %>% purrr::reduce(., left_join, by = "Gene") -> df_
    head(df_)
    
    df_ %>% 
      ggplot(., aes(x = Mean, y = Nonzero_Prop)) + 
#      geom_point(size = 0.1) + 
      geom_hex(bins = 100) + 
      scale_fill_continuous(trans = 'log', name = "Count (log scale)") +
#      ggtitle(sprintf("%s", method_)) + 
      xlab(sprintf("-log10(%s):mean", "p.adj")) + 
      ylab(sprintf("-log10(%s):nonzero prop.", "p.adj")) + 
      theme_classic() -> p
    
    xintercept_ = df_$Mean %>% median(., na.rm = T)
    yintercept_ = df_$Nonzero_Prop %>% median(., na.rm = T)
    
    p + geom_vline(xintercept = xintercept_, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = yintercept_, linetype = "dashed", color = "red") -> p
    
    df_ %>% drop_na() -> df_
    cor_ = cor(df_$Mean, df_$Nonzero_Prop, method = "pearson")
    p + 
      annotate(geom = "text", x = 3, y = 13, label = sprintf("r = %.2f", cor_)) -> p
    
  }) %>% .[[1]] -> p_pval_corr #Fig_3E

list(Mean = res_RNA_Mean, Nonzero_Prop = res_RNA_NonZeroProp) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    type_ = y
    df_ = x
    df_ %>% dplyr::select(Gene, F24_Phase) %>% 
      "colnames<-"(., c("Gene", type_))
  }) %>% purrr::reduce(., left_join, by = "Gene") %>% 
  dplyr::mutate(Nonzero_Prop = case_when(
#    Mean - Nonzero_Prop > 12 ~ Nonzero_Prop+24,
#    Mean - Nonzero_Prop < -12 ~ Nonzero_Prop-24,
    Nonzero_Prop - Mean > 12 ~ Nonzero_Prop-24,
    TRUE ~ Nonzero_Prop
  )) %>% 
  dplyr::mutate(Mean = case_when(
    Nonzero_Prop - Mean < -12 ~ Mean - 24, 
    TRUE ~ Mean
  )) %>% 
  list(phaes = .) %>% 
  map(function(x){
    df_ = x
    df_ %>% drop_na() -> df_
    cor(df_$Mean, df_$Nonzero_Prop, method = "pearson") -> cor_
    df_ %>% 
      ggplot(aes(x = Mean, y = Nonzero_Prop)) + 
      geom_hex(bins = 100) -> p
    p + 
      theme_classic() +
      ylab("Adjusted phase: nonzero prop.") + 
      xlab("Adjusted phase: mean") -> p
    p + 
      annotate(geom = "text", x = 1, y = 24, label = sprintf("r = %.2f", cor_)) -> p
    p + 
      geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") -> p
    p + scale_x_continuous(breaks = seq(-10, 20, 5)) -> p
    p + scale_y_continuous(breaks = seq(-10, 20, 5)) -> p
  }) -> p_phase_hex  #Fig_3E

list(Mean = res_RNA_Mean, Nonzero_Prop = res_RNA_NonZeroProp) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    type_ = y
    df_ = x
    df_ %>% dplyr::select(Gene, MetaCycle_JTK_amplitude) %>% 
      mutate(MetaCycle_JTK_amplitude = log(MetaCycle_JTK_amplitude)) %>%
      "colnames<-"(., c("Gene", type_))
  }) %>% purrr::reduce(., left_join, by = "Gene") %>% 
  dplyr::filter(is.finite(Mean)&is.finite(Nonzero_Prop)) %>% 
  list(amp = .) %>% 
  map(function(x){
    df_ = x
    df_ %>% drop_na() -> df_
    cor(df_$Mean, df_$Nonzero_Prop, method = "pearson") -> cor_
    df_ %>% 
      ggplot(aes(x = Mean, y = Nonzero_Prop)) + 
      geom_hex(bins = 100) -> p
    p + 
      theme_classic() +
      ylab("log(amplitude): nonzero prop.") + 
      xlab("log(amplitude): mean") -> p
    p + 
      annotate(geom = "text", x = -11, y = -2, label = sprintf("r = %.2f", cor_)) -> p
  }) -> p_amp_hex #Fig_3E

# 13. smFISH (Arntl & Per1) ----
library(mclust)
library(tidyverse)
list_tmp = list()
list_tmp_1 = list()
list(Arntl = 1, Per1 = 2) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = y 
    channel_ = x
    c("_ZT0_", "_ZT03_", "_ZT06_", "_ZT09_", "_ZT12_", "_ZT15_", "_ZT18_", "_ZT21_") %>% 
      "names<-"(.,.) %>% 
#      .[3] %>% 
      purrr::map(function(x){
        ZT_ = x
        list.files("~/Dropbox/singulomics/RNA_FISH/Bmal1_Per1_FISH/Table", pattern = "csv", full.names = T) %>% 
          {.[grepl(ZT_,.)]} -> files_
        print(files_)
        ZT_ = gsub("_", "", ZT_)
        i = 0
        files_ %>% 
#          .[1] %>% 
          purrr::map(function(x){
            i <<- i +1
            read.csv(x, header = T, stringsAsFactors = F) -> df_
            df_ %>% dplyr::filter(Channel == channel_) -> df_
            IntDen = df_$IntDen
            df_$NucIdx %>% table() %>% {.[names(.) != 0]} %>% sum() -> nuc_fish_signal
            
            model <- Mclust(IntDen, G = 2)
            list_tmp_1[[gene_]][[ZT_]][[i]] <<- model
            plot(model, what = "classification") 
            clusters <- model$classification
            unique(clusters) %>% sort() %>% rev() %>% .[1] -> last_clust
            print(last_clust)
            (clusters %in% last_clust) -> se_
            list_tmp[[gene_]][[ZT_]][[i]] <<- df_ %>% dplyr::mutate(cluster = case_when(se_ ~ 2, TRUE ~ 1))
            
            print(mean(df_[se_,] %>% .$IntDen %>% mean()))
            print(mean(df_[!se_,] %>% .$IntDen %>% mean()))
            
            df_[se_,] -> df_
            total_nuc = df_$Tot_nuc %>% unique()
            
            df_$NucIdx %>% table() %>% {.[names(.) != 0]} %>% length() -> nonzero_nuc_num
#            df_$NucIdx %>% table() %>% {.[names(.) != 0]} %>% sum() -> nuc_fish_signal
            data.frame(Gene = gene_, ZT = ZT_, nonzero_nuc = nonzero_nuc_num, total_nuc = total_nuc, nuc_fish_signal = nuc_fish_signal)
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .) -> df_
    df_ %>% dplyr::mutate(nonzero_prop = nonzero_nuc/total_nuc, avg_nuc_fish_signal = nuc_fish_signal/total_nuc) -> df_
  }) -> list_summary

list_summary %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = y
    df_ = x
    df_$ZT = gsub("ZT(\\d+)", "\\1", df_$ZT) %>% as.integer()
    df_$avg_nuc_fish_signal = log(df_$avg_nuc_fish_signal, 2)
    df_ %>% 
      group_by(Gene, ZT) %>% 
      summarise(sd_nonzero_prop = sd(nonzero_prop), mean_nonzero_prop = mean(nonzero_prop), 
                sd_nuc_fish_signal = sd(avg_nuc_fish_signal), mean_nuc_fish_signal = mean(avg_nuc_fish_signal)) -> df_
    df_ %>% 
      dplyr::bind_rows(df_ %>% dplyr::filter(ZT == 0) %>% dplyr::mutate(ZT = 24)) -> df_
    
    # Minimum and maximum of mean_nonzero_prop
    min_prop <- min(df_$mean_nonzero_prop)
    max_prop <- max(df_$mean_nonzero_prop)
    
    # Minimum and maximum of mean_nuc_fish_signal
    min_nuc <- min(df_$mean_nuc_fish_signal)
    max_nuc <- max(df_$mean_nuc_fish_signal)
    
    # Calculate 'a'
    a <- (max_nuc - min_nuc) / (max_prop - min_prop)
    
    # Calculate 'b'
    b <- min_nuc - a * min_prop
    
    df_$mean_nonzero_prop_scaled <- a * df_$mean_nonzero_prop + b
    df_$sd_nonzero_prop_scaled <- a * df_$sd_nonzero_prop
    
    df_
    
    df_ %>% 
      ggplot(aes(x = ZT, y = mean_nuc_fish_signal, group = Gene)) + 
      geom_line(color = "#56B4E9") + 
      geom_line(aes(y = mean_nonzero_prop_scaled), color = "#009E73") +
      scale_y_continuous(
        name = "log2[nuclear mRNA signal]",
        # Left y-axis label
        sec.axis = sec_axis(
          trans = ~ (. - b) / a,   # Inverse transformation
          name = "Nonzero prop."
        )
      ) + 
      geom_ribbon(aes(ymin = mean_nuc_fish_signal - sd_nuc_fish_signal, ymax = mean_nuc_fish_signal + sd_nuc_fish_signal), fill = "#56B4E9", alpha = 0.2) +
      geom_ribbon(aes(ymin = mean_nonzero_prop_scaled - sd_nonzero_prop_scaled, ymax = mean_nonzero_prop_scaled + sd_nonzero_prop_scaled), fill = "#009E73", alpha = 0.2) +
      ggtitle(gene_)
  }) %>% patchwork::wrap_plots(., ncol = 2) -> p_fish #Fig_3F

list_tmp %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = y
    x %>%
      purrr::map2(.x=.,.y=names(.),.f=function(x,y){
        ZT_ = y
        list_ = x
        1:length(list_) %>% 
          "names<-"(., sprintf("image_%s", .)) %>% 
          purrr::map(function(i){
            
            list_tmp_1[[gene_]][[ZT_]][[i]] -> model
            model$parameters$mean -> mean_
            model$parameters$variance$sigmasq -> variance_
            model$parameters$pro -> proportions_
            
            x_min <- min(mean_) - 3 * sqrt(max(variance_))
            x_max <- max(mean_) + 3 * sqrt(max(variance_))
            x_values <- seq(x_min, x_max, length.out = 1000)
            
            density_data <- data.frame(x = x_values)
            for (k in 1:length(mean_)) {
              density_data[[paste0("Component", k)]] <- proportions_[k] * dnorm(x_values, mean = mean_[k], sd = sqrt(variance_[k]))
            }
            density_data$TotalDensity <- rowSums(density_data[ , paste0("Component", 1:length(mean_))])
            density_data
            
            density_long <- gather(density_data, key = "Component", value = "Density", -x)
            density_long %>% dplyr::filter(Component != "TotalDensity") -> density_long
            density_long %>% dplyr::filter(x >= 0) -> density_long
            
            list_[[i]] -> df_
            df_$cluster = factor(df_$cluster)
            
            df_ %>% 
              ggplot(aes(x = IntDen)) + 
#              geom_density() +
              geom_histogram(aes(y = ..density..), bins = 50, fill = "white", alpha = 0.3, position = "dodge", color = "black") + 
              geom_line(data = density_long, aes(x = x, y = Density, color = Component), size = 1) + 
              ggtitle(sprintf("%s %s", gene_, ZT_)) + 
              theme_minimal()
          })
      })
  }) -> p_list_5

c("Arntl", "Per1") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(x){
    gene_ = x
    c("ZT0", "ZT03", "ZT06", "ZT09", "ZT12", "ZT15", "ZT18", "ZT21") %>% 
      purrr::map(function(x){
        ZT_ = x
        p_list_5[[gene_]][[ZT_]][["image_1"]] -> p
      }) -> p_list
    patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p
  }) -> p_list_6

p_list_6$Arntl #Fig_S28A
p_list_6$Per1 #Fig_S28B

# 14. Burst frequency (Nr1d1, Cry1, Bmal1) ----
setwd("~/Dropbox/Research/singulomics/")
# Burst frequency and size from Phillips et al. pulled by Jay Suh
burst.freq=read.csv('bursting/smFISH_burst_freq.csv')
colnames(burst.freq)[1]='Time'

library(ggplot2)
p1=ggplot(burst.freq, aes(x=Time))+
  geom_ribbon(aes(x=Time, ymin=Nr1d1_pred_lower, ymax=Nr1d1_pred_upper), fill=adjustcolor( "#CC6666", alpha.f = 0.2)) +
  geom_line(aes(y=Nr1d1_pred_median, colour='#CC6666'))+
  geom_ribbon(aes(x=Time, ymin=Bmal1_pred_lower, ymax=Bmal1_pred_upper), fill=adjustcolor( "#9999CC", alpha.f = 0.2)) +
  geom_line(aes(y=Bmal1_pred_median, colour='#9999CC'))+
  geom_ribbon(aes(x=Time, ymin=Cry1_pred_lower, ymax=Cry1_pred_upper), fill=adjustcolor( "#66CC99", alpha.f = 0.2)) +
  geom_line(aes(y=Cry1_pred_median, colour='#66CC99'))

burst.size=read.csv('bursting/smFISH_burst_size.csv')
burst.size=as.matrix(burst.size[,-1])

hist(burst.size, xlim=c(0,20), 100)
burst.size.toplot=reshape2::melt(burst.size)

p2=ggplot(burst.size.toplot, aes(x=value, fill=Var2))+geom_histogram(bins=100)


mean=as.matrix(read.csv('bursting/smFISH_mean_mRNA_cell.csv'))
colnames(mean)[1]='Time'

p3=ggplot(mean, aes(x=Time))+
  geom_line(aes(y=Nr1d1, colour='#CC6666'))+
  geom_line(aes(y=Bmal1, colour='#9999CC'))+
  geom_line(aes(y=Cry1, colour='#66CC99'))

Nr1d1.dot=read.csv('bursting/smFISH_mean_Nr1d1_dots.csv')
Cry1.dot=read.csv('bursting/smFISH_mean_Cry1_dots.csv')
Bmal1.dot=read.csv('bursting/smFISH_mean_Bmal1_dots.csv')

p3=p3+geom_point(data=Nr1d1.dot, aes(x=time, y=avg, color='#CC6666'), size=1, alpha=0.9)+
  geom_point(data=Bmal1.dot, aes(x=time, y=avg, color='#9999CC'), size=1, alpha=0.9)+
  geom_point(data=Cry1.dot, aes(x=time, y=avg, color='#66CC99'), size=1, alpha=0.9)



library(patchwork)
p=p3+p1+p2
p
ggsave(filename='bursting/bursting_mean.pdf', width=13, height=3.5, plot=p) #Fig_S15
####

# 16. Pair-wise correlation (Nonzero_mean, gini, Nonzero_Prop, cv, and mean) ----
setwd("~/Dropbox/Research/singulomics/")

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # mm10
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(biovizBase)
library(SingleCellExperiment)
library(mbkmeans)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(MetaCycle)
library(future)
library(furrr)
library(DescTools)

source('github/runJTK.R')


# Load in the integrated and QC'ed data
readRDS("github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells

colnames(sc) %>% length()
length(hepatocytes_cells)

sc <- sc[,hepatocytes_cells]
sc$celltype %>% table()

DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "multicca.umap", repel = TRUE)
DefaultAssay(sc) <- "SCT"

# Choose resolution: do not need too many metacells, especially because we need to characterize the distributions of cells within the metacells
resolution=0.5
sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=resolution)
DimPlot(sc, group.by = "seurat_clusters", label = TRUE, reduction = "multicca.umap", repel = TRUE)
sc@meta.data$seurat_clusters %>% unique() %>% length()

sc=sc[,!grepl('KO',sc$group)] # Remove the knockout to look at circadian rhythmicity
sc$group=droplevels(sc$group)
table(sc$group)

sc$ZT=as.numeric(gsub('ZT','',sc$ZT)) # Change ZT to be numeric time
sc=sc[,sc$celltype=='Hepatocytes'] # Keep only hepatocytes
table(sc$ZT)
table(sc$celltype)

sc$cc_cluster=NULL
sc$cc_clusters=sc$seurat_clusters
table(sc$cc_clusters)
table(sc$cc_clusters, sc$group)
dim(sc)

# Minimum number of cells per timepoint in each cluster: 50
sc=sc[,sc$cc_clusters %in% names(which(apply(table(sc$cc_clusters, sc$group), 1, min)>50))]
sc$cc_clusters=droplevels(sc$cc_clusters)
table(sc$cc_clusters, sc$group)
length(table(sc$cc_clusters)) # Total number of replicates after filtering
sc$cluster=sc$cc_clusters

# We will only focus on RNA for now
genes=rownames(sc@assays$RNA)
cyc.genes=c('Arntl','Per1','Per2','Npas2','Dbp','Cry1','Clock','Nr1d1','Rorc')
cyc.genes %in% genes
# Nr1d1 (peak ~ZT06)
# Dbp (peak ~ZT06)
# Arntl (peak ~ZT22)
# Cry1 (peak ~ZT18/22)
# Ciart (peak ~ZT06/10)
# Rorc (peak ~ZT18)

output.rna.list=vector(mode = "list", length=length(genes))
names(output.rna.list)=genes
i=which(genes=='Arntl')
for(i in 1:length(genes)){
  if(i%%100==0) cat(i,'\t')
  gene.name=genes[i]
  gene.count=sc@assays$RNA@counts[gene.name,] # gene expression vector
  gene.count=gene.count/sc$nCount_RNA*median(sc$nCount_RNA) # Adjust for library size
  
  cell.meta=c('celltype','ZT','cluster')
  
  # output=aggregate(gene.count, by=sc@meta.data[,cell.meta],
  #           FUN=function(x) c(mean=mean(x),
  #                             nonzero.prop=mean(x>0),
  #                             nonzero.mean=mean(x[x>0])))
  
  output=aggregate(gene.count, by=sc@meta.data[,cell.meta],
            FUN=function(x) c(mean=mean(x),
                              sd=sd(x),
                              gini=Gini(x),
                              nonzero.prop=mean(x>0),
                              nonzero.mean=mean(x[x>0])))
  #   burst.freq=moment_fun_kon(x)
  
  output=do.call(data.frame,output)
  colnames(output)=gsub('x.','', colnames(output))
  
  output$cv=(output$sd)/(output$mean)
  #output$var=(output$sd)^2 # Use sd instead
  #output$sd.resid=lm(output$sd~output$mean)$residuals # No signals
  output[,-(1:3)]=round(output[,-(1:3)],4)
  
  output.rna.list[[i]]=output
}
output.rna.list[[1]] # 6 timepoints x number of replicates

for(FUN in colnames(output.rna.list[[1]])[-(1:3)]){
  cat('Running JTK on', FUN,'...\n')
  eval(parse(text=paste0('JTK.rna.',FUN,'=runJTK.FUN(output.rna.list, group.by=\'celltype\', cell.meta=cell.meta, FUN=\'',FUN,'\')[[1]]')))
}

save.image(file='~/Dropbox/Research/singulomics/github_rda/hepatocytes_proportion_v2.rda')
####

####
setwd("~/Dropbox/Research/singulomics/")

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # mm10
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(biovizBase)
library(SingleCellExperiment)
library(mbkmeans)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(MetaCycle)
library(future)
library(furrr)
library(DescTools)


load('~/Dropbox/Research/singulomics/github_rda/hepatocytes_proportion_v2.rda')
cyc.genes2=cyc.genes

cyc.genes=read.csv('bursting/NIHMS299852-supplement-Supplemental_tables.csv')
cyc.genes=cyc.genes$Gene.Symbol[cyc.genes$BH.Q< 0.01]
cyc.genes=unique(cyc.genes)
cyc.genes=cyc.genes[cyc.genes %in% rownames(sc)]
length(cyc.genes)
#cyc.genes=cyc.genes[1:500]

toplot=output.rna.list[[cyc.genes[1]]][,-(1:3)]
for(i in cyc.genes[-1]){
  toplot=rbind(toplot, output.rna.list[[i]][,-(1:3)])
}

toplot=toplot[,c('mean','nonzero.prop','nonzero.mean','gini','cv')]

panel.cor <- function(x,y,digits = 2,prefix = "",cex.cor = NULL,...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  filter = !is.na(y) &
    !is.na(x) &
    !is.nan(x) & !is.nan(y) & !is.infinite(x) & !is.infinite(y)
  r <- cor(x[filter], y[filter], method='spearman')
  r.temp = abs(r)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, 'r = ', txt)
  cex.cor <- 1
  text(0.5, 0.5, txt, cex = cex.cor * r.temp^(1/4) * 1.5)
}

cyc.genes=cyc.genes2
rm(cyc.genes2)

pdf(file='bursting/pairwise.pdf', width=8, height=8)
pairs(toplot, lower.panel = function(x,y){smoothScatter(x,y,add=TRUE)}, upper.panel = panel.cor, main='Pairwise Spearman correlation') #Fig_S16B
dev.off()

####

# 17. Validation mouse SCN ----

#Generate the expression matrix
list.files("~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN", pattern = "\\.csv", full.names = T) %>% 
  {.[grep("CT", .)]} %>% 
  gtools::mixedsort() -> files_

ZT_list = list()
files_ %>% 
  {
    files_ = .
    seq(14,58,4) %>% sprintf("CT%s", .) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        CT_ = x
        files_ %>% {.[grepl(CT_, .)]} -> file_
        print(file_)
        read.csv(file_, header = T, stringsAsFactors = F) -> dat_
        colnames(dat_) ->> ZT_list[[CT_]]
        dat_ %>% rownames_to_column("Gene") -> dat_
      }) %>% purrr::reduce(., full_join, by = "Gene") %>% 
      column_to_rownames("Gene") %>% 
      dplyr::mutate(across(everything(), ~ replace_na(., 0))) -> dat_
    dat_
  } -> dat_

save(dat_, file = "~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/exprssion_matrix.rda")

#Generate the seurat object 
sc = Seurat::CreateSeuratObject(counts = dat_, project = "GSE117295_RAW_Mouse_SCN", assay = "RNA")
##remove cells with duplicate barcodes
colnames(sc) %>% {.[grepl("\\.x|\\.y", .)]} -> duplicated_barcode
barcode_kept = setdiff(colnames(sc), duplicated_barcode)
sc = sc[,barcode_kept]

#Generate celltype metadata
read.table("~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/cluster_vs_barcode/all_cell_cluster_cell_barcode.tsv", 
           header =T, sep = "\t", stringsAsFactors = F) %>% as_tibble() -> barcode_df
ZT_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    CT_ = y
    barcode_ = x
    barcode_ %>% 
      as_tibble() %>% 
      "colnames<-"(., "barcode") %>% 
      dplyr::mutate(CT = CT_)
  }) %>% do.call(rbind, .) %>% 
  left_join(x = barcode_df, y = ., by = "barcode") -> meta_df

save(meta_df, file = "~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/metadata.rda")

##further filter sc
sc@meta.data %>% 
  rownames_to_column("barcode") %>% 
  left_join(x =. , y = meta_df, by = "barcode") %>% 
  .$barcode %>% table() %>% {.[. > 1]} %>% names() -> duplicated_barcode
sc = sc[,setdiff(colnames(sc), duplicated_barcode)]

sc@meta.data %>% 
  rownames_to_column("barcode") %>% 
  left_join(x =. , y = meta_df, by = "barcode") %>% 
  column_to_rownames("barcode") -> sc@meta.data

sc@meta.data %>% dplyr::filter(is.na(cell_type)) %>% rownames() -> non_annotaed_barcode
sc = sc[,setdiff(colnames(sc), non_annotaed_barcode)]
save(sc, file = "~/Dropbox/singulomics/github_rda/GSE117295_RAW_Mouse_SCN/seurat_object.rda")

sc <- SCTransform(sc, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
p1 <- DimPlot(sc, reduction = "umap.rna", group.by='cell_type', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA")  #Fig_S17A

sc@meta.data %>% dplyr::mutate(celltype = case_when(
  grepl("Neurons", cell_type) ~ "Neurons",
  TRUE ~ cell_type
)) -> sc@meta.data
sc$celltype %>% table()

# Generate RNA expression matrix
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) -> RNA_norm

# 5 pseudo bulk (random sampling)
c("Neurons1", "Ependymocytes") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  map(function(x){
    celltype_ = x
    print(celltype_)
    sc@meta.data %>% dplyr::filter(cell_type == celltype_) -> meta_df
    seq(14,58,4) %>% sprintf("CT%s", .) %>% 
      "names<-"(.,.) %>% 
#      .[1:2] %>% 
      map(function(x){
        set.seed(123)
        CT_ = x
        meta_df %>% dplyr::filter(CT == CT_) %>% 
          dplyr::mutate(group = sample(1:5, nrow(.), replace = T)) -> meta_df
        c(1:5) %>% 
          map(function(x){
            i = x
            meta_df %>% dplyr::filter(group == i) -> meta_df_
            barcode = rownames(meta_df_)
            RNA_norm[,barcode] -> RNA_norm
            rowMeans(RNA_norm) %>% 
              as.data.frame() %>% 
              "colnames<-"(., sprintf("%s_REP%s", CT_, i)) -> RNA_nrom
          }) %>% do.call(cbind, .) -> RNA_norm
#        head(RNA_norm)
      }) %>% do.call(cbind, .) %>% 
      "colnames<-"(., gsub("^CT\\d+\\.CT", "ZT", colnames(.))) -> RNA_norm
#    head(RNA_norm)
  }) -> list_RNA_norm

c("Neurons1", "Ependymocytes") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  map(function(x){
    celltype_ = x
    print(celltype_)
    sc@meta.data %>% dplyr::filter(cell_type == celltype_) -> meta_df
    seq(14,58,4) %>% sprintf("CT%s", .) %>% 
      "names<-"(.,.) %>% 
      #      .[1:2] %>% 
      map(function(x){
        set.seed(123)
        CT_ = x
        meta_df %>% dplyr::filter(CT == CT_) %>% 
          dplyr::mutate(group = sample(1:5, nrow(.), replace = T)) -> meta_df
        c(1:5) %>% 
          map(function(x){
            i = x
            meta_df %>% dplyr::filter(group == i) -> meta_df_
            barcode = rownames(meta_df_)
            RNA_norm[,barcode] -> RNA_norm
#            rowMeans(RNA_norm) %>% 
            apply(RNA_norm, 1, function(x){sum(x != 0)/length(x)}) %>% 
              as.data.frame() %>% 
              "colnames<-"(., sprintf("%s_REP%s", CT_, i)) -> RNA_non_zero
          }) %>% do.call(cbind, .) -> RNA_non_zero
        #        head(RNA_norm)
      }) %>% do.call(cbind, .) %>% 
      "colnames<-"(., gsub("^CT\\d+\\.CT", "ZT", colnames(.))) -> RNA_non_zero
    #    head(RNA_norm)
  }) -> list_non_zero_prop

# Calculate gene rhythmicity (Mean and non zero proportion)
source("~/Dropbox/singulomics/github/Calculate_HMP.R")
list_RNA_norm %>% 
  map(function(x){
    df_ = x 
    gene_ = rownames(df_)
    df_ = df_ %>% as_tibble()
    timepoints_ = colnames(df_) %>% gsub("ZT(\\d+)_.+", "\\1", .) %>% as.integer()
    res_ = cyclic_HMP(timepoints = timepoints, gene = gene_, exp_matrix = df_)
  }) -> list_RNA_norm_pval
list_RNA_norm_pval %>% 
  map(function(x){
    #    class(x)
    df_ = x
    df_ %>% dplyr::select(!matches("LS|meta2d")) -> df_
    recal_cauchy_p_and_hmp(df_) -> df_
  }) -> list_RNA_norm_pval

list_non_zero_prop %>% 
  map(function(x){
    df_ = x 
    gene_ = rownames(df_)
    df_ = df_ %>% as_tibble()
    timepoints_ = colnames(df_) %>% gsub("ZT(\\d+)_.+", "\\1", .) %>% as.integer()
    res_ = cyclic_HMP(timepoints = timepoints, gene = gene_, exp_matrix = df_)
  }) -> list_non_zero_prop_pval
list_non_zero_prop_pval %>% 
  map(function(x){
#    class(x)
    df_ = x
    df_ %>% dplyr::select(!matches("LS|meta2d")) -> df_
    recal_cauchy_p_and_hmp(df_) -> df_
  }) -> list_non_zero_prop_pval

# Genome-wide comparision (pvalues, phase and amplitude)
#5 Genomewide comparison 
data.frame(
  Mean = list_RNA_norm_pval$Neurons1$cauchy_BH.Q,
  Nonzero = list_non_zero_prop_pval$Neurons1$cauchy_BH.Q
) %>% dplyr::filter((!is.na(Mean))|(!is.na(Nonzero))) %>% 
  {-log(., 10)} %>% 
  {
    df_ = .
    cor_ = cor(df_$Mean, df_$Nonzero, method = "pearson")
    df_ %>% 
    ggplot(aes(x = Mean, y = Nonzero)) + 
#      geom_point(size = 0.1) + 
      geom_hex(bins = 100) + 
      scale_fill_continuous(trans = 'log', name = "Count (log scale)") +
      annotate("text", x = 1, y = 10, label = sprintf("cor = %.2f", cor_)) + 
      xlab("-log10(p.adj):mean") + 
      ylab("-log10(p.adj):nonzero prop.") -> p
    
    xintercept_ = df_$Mean %>% median(., na.rm = T)
    yintercept_ = df_$Nonzero %>% median(., na.rm = T)
    
    p + geom_vline(xintercept = xintercept_, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = yintercept_, linetype = "dashed", color = "red") -> p
    p + theme_classic()
  } -> p1

data.frame(
  Mean = list_RNA_norm_pval$Neurons1$HR_phi,
  Nonzero = list_non_zero_prop_pval$Neurons1$HR_phi
) %>% dplyr::filter((!is.na(Mean))|(!is.na(Nonzero))) %>% 
  dplyr::mutate(across(everything(), function(x){(x/(2*pi))*24})) %>% 
  dplyr::mutate(Nonzero = case_when(
#    Mean - Nonzero > 12 ~ Nonzero+24,
#    Mean - Nonzero < -12 ~ Nonzero-24,
    Nonzero - Mean > 12 ~ Nonzero - 24,
    TRUE ~ Nonzero
  )) %>% 
  dplyr::mutate(Mean = case_when(
    Nonzero - Mean < -12 ~ Mean - 24, 
    TRUE ~ Mean
  )) %>% 
  {
    df_ = .
    print(nrow(df_))
    cor_ = cor(df_$Mean, df_$Nonzero, method = "pearson")
    df_ %>%
      ggplot(aes(x = Mean, y = Nonzero)) + 
      geom_hex(bins = 100) + 
      scale_fill_continuous(trans = 'log', name = "Count (log scale)") +
      geom_abline(intercept = 0, slope = 1, color = "red") + 
      annotate("text", x = 2, y = 24, label = sprintf("cor = %.2f", cor_)) + 
      ylab("Adjusted phase: nonzero prop.") + 
      xlab("Adjusted phase: mean") + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(-10, 20, 5)) +
      scale_y_continuous(breaks = seq(-10, 20, 5))
  } -> p2

data.frame(
  Mean = list_RNA_norm_pval$Neurons1$MetaCycle_JTK_amplitude,
  Nonzero = list_non_zero_prop_pval$Neurons1$MetaCycle_JTK_amplitude
) %>% dplyr::filter((!is.na(Mean))|(!is.na(Nonzero))) %>% 
  {log(.)} %>% 
  dplyr::filter((is.finite(Mean))&(is.finite(Nonzero))) %>% 
  {
    df_ = .
#    df_ = df_ %>% dplyr::filter(Nonzero > -30) -> df_
    df_ = df_ %>% dplyr::filter(Nonzero > quantile(Nonzero, 0.01)) -> df_
    print(nrow(df_))
    cor_ = cor(df_$Mean, df_$Nonzero, method = "pearson")
    cor_
    df_ %>%
      ggplot(aes(x = Mean, y = Nonzero)) + 
      geom_hex(bins = 100) + 
#      scale_fill_continuous(trans = 'log', name = "Count (log scale)") +
#      geom_abline(intercept = 0, slope = 1, color = "red") + 
      annotate("text", x = -12, y = -2.5, label = sprintf("cor = %.2f", cor_)) + 
      ylab("log(amplitude): nonzero prop.") + 
      xlab("log(amplitude): mean") + 
      theme_classic()
  } -> p3

patchwork::wrap_plots(p1,p2,p3, ncol = 3) #Fig_3E

#Venn diagram of cyclic genes (mean vs nonzero)
list(
  Mean = list_RNA_norm_pval$Neurons1 %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene,
  Nonzero = list_non_zero_prop_pval$Neurons1 %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene
) %>% 
  ggvenn::ggvenn() #Fig_S17C

list(
  Mean = list_RNA_norm_pval$Ependymocytes %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene,
  Nonzero = list_non_zero_prop_pval$Ependymocytes %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene
) %>% 
  ggvenn::ggvenn()

#7. boxplot of mean and nonzero prop
list(
  Mean = list_RNA_norm_pval$Neurons1 %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene,
  Nonzero = list_non_zero_prop_pval$Neurons1 %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% .$Gene
) %>% 
  ggvenn::list_to_data_frame() %>% 
  {
    df_ = .
    Mean = df_ %>% dplyr::filter((Mean==TRUE)&(Nonzero==FALSE)) %>% .$key
    Nonzero = df_ %>% dplyr::filter((Mean==FALSE)&(Nonzero==TRUE)) %>% .$key
    Overlapped = df_ %>% dplyr::filter((Mean==TRUE)&(Nonzero==TRUE)) %>% .$key
    
    list(Mean = Mean, Nonzero = Nonzero, Overlapped = Overlapped)
  } -> list_genes_venn

c("Neurons1", "Ependymocytes") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    celltype_ = x
    list_val = list()
    list_genes_venn %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        gene_ = x
        group_ = y
        #    list_non_zero_prop$Neurons1[gene_, ] -> df_
        list_non_zero_prop[[celltype_]][gene_, ] -> df_
        unlist(df_) -> df_
        df_ ->> list_val[[group_]]
        df_ %>% 
          as_tibble() %>% 
          dplyr::mutate(group = group_)
      }) %>% do.call(rbind, .) %>% 
      {
        df_ = .
        df_ %>% dplyr::filter(group != "Overlapped") -> df_
        
        df_ %>% dplyr::filter(group == "Mean") %>% .$value %>% quantile(., 0.75) %>% .[[1]] %>% {. + 1} -> low_break
        df_ %>% dplyr::filter(group == "Mean") %>% .$value %>% quantile(., 0.99) %>% .[[1]] -> up_break
        print(low_break)
        
        df_ %>% 
          #    dplyr::mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
          dplyr::mutate(group = factor(group, levels = c("Mean", "Nonzero"))) %>% 
          ggplot(aes(x = group, y = value, fill = group)) + 
          geom_boxplot(width = 0.4) + 
          scale_fill_manual(values = c(Mean="#56B4E9", Nonzero="#009E73")) + 
          theme_classic() + 
          xlab(NULL) + 
          ylab("Non zero proportion") + 
          theme(legend.position = "none") + 
          scale_y_break(breaks = c(0.4,0.5), scale = 0.3) + 
#          scale_y_break(breaks = c(low_break ,up_break), scale = 0.3) + 
          ggtitle(celltype_)
      } -> p1
    print(t.test(list_val$Mean, list_val$Nonzero))
    return(p1)
  }) -> p1
#p1$Neurons1+p1$Ependymocytes
p1$Neurons1 #Fig_S17D


c("Neurons1", "Ependymocytes") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    celltype_ = x
list_val = list()
list_genes_venn %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = x
    group_ = y
    list_RNA_norm[[celltype_]][gene_, ] -> df_
    unlist(df_) -> df_
    df_ ->> list_val[[group_]]
    df_ %>% 
      as_tibble() %>% 
      dplyr::mutate(group = group_)
  }) %>% do.call(rbind, .) %>% 
  {
    df_ = .
    print(median(df_$value))
    df_ %>% dplyr::filter(group != "Overlapped") -> df_
    df_ %>% dplyr::filter(group == "Mean") %>% .$value %>% quantile(., 0.75) %>% .[[1]] %>% {. + 1} -> low_break
    df_ %>% dplyr::filter(group == "Mean") %>% .$value %>% quantile(., 0.99) %>% .[[1]] -> up_break
    sprintf("low_break: %f, up_break: %f", low_break, up_break) %>% print()
    df_ %>% 
#      dplyr::mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
      dplyr::mutate(group = factor(group, levels = c("Mean", "Nonzero"))) %>% 
      ggplot(aes(x = group, y = value, fill = group)) + 
      geom_boxplot(width = 0.4) + 
      scale_fill_manual(values = c(Mean="#56B4E9", Nonzero="#009E73")) + 
      theme_classic() + 
      xlab(NULL) + 
      ylab("Normalized expression") + 
      theme(legend.position = "none") + 
#      scale_y_break(breaks = c(0.3,9), scale = 0.3) + 
      scale_y_break(breaks = c(0.25, up_break), scale = 0.3) + 
      ggtitle(celltype_)
  } -> p2
print(t.test(list_val$Mean, list_val$Nonzero))
return(p2)
}) -> p2
#p2$Neurons1+p2$Ependymocytes
p2$Neurons1 #Fig_S17D

# 18. Validation drosophila brain ----
library(tidyverse)
library(Seurat)

#Load data
list.files("~/Dropbox/singulomics/github_rda/GSE157504_RAW", pattern = "csv", full.names = T) %>% 
  map(function(x){
    read.csv(file = x, header = T, sep = ",", stringsAsFactors = F) -> df_
  }) -> df_list

purrr::reduce(df_list, full_join, by = "X") -> df_rna_all
df_rna_all %>% column_to_rownames(var = "X") %>% 
  as.data.frame() -> df_rna_all
rm(df_list)
gc()

df_rna_all %>% replace(is.na(.), 0) -> df_rna_all

CreateSeuratObject(counts = df_rna_all, project = "drosophila_scRNA", assay = "RNA") -> scRNA_obj
rm(df_rna_all)

sc = scRNA_obj
rm(scRNA_obj)

colnames(sc)

sc@meta.data %>% mutate(AR = gsub(pattern = ".+(AR\\d+?)_.+", replacement = "\\1", x = rownames(.))) -> sc@meta.data
sc@meta.data %>% rownames() %>% gsub(pattern = "^.+?_.+?_((.+?)_.+?)_.+", replace = "\\2", x = .) -> sc@meta.data$experimental_design
sc@meta.data %>% rownames() %>% gsub(pattern = "^.+?_.+?_((.+?)_(.+?))_.+", replace = "\\3", x = .) -> sc@meta.data$time
sc@meta.data %>% rownames() %>% gsub(pattern = "^X(.+?)_.+", replace = "\\1", x = .) -> sc@meta.data$date
sc@meta.data <- sc@meta.data %>% mutate(orig.ident = "scRNA_drosophila")

f.entropy <- function( m.in ){
  m.in <- apply( m.in, 2, function( x ){ x / sum(x, na.rm=T ); } );
  v.entropy <- apply( m.in, 2, function(x){ y <- unname(x); y <- y[ y > 0 ]; sum( -y * log(y) ); } );
  return( v.entropy );
}

sc@meta.data$entropy = f.entropy(sc@assays$RNA@counts)

sc_QC <- subset(
  x = sc,
  subset = nCount_RNA < 75000 &
    nCount_RNA > 6000 &
    nFeature_RNA < 6000 &
    nFeature_RNA > 1000 &
    entropy > 5.5
)

sc = sc_QC
rm(sc_QC)

sc@meta.data$log2nCount_RNA <- log2(sc@meta.data$nCount_RNA+1)
sc@meta.data$log2nFeature_RNA <- log2(sc@meta.data$nFeature_RNA+1)

library(gtools)
#structure containing information about data
l.info <- list();

#compute genes not to be used in clustering
l.info[["geneSets"]][["mt"]] <- v.mt.genes <- grep("^mt:", rownames(sc@assays$RNA@data), value=T, ignore.case=T)
l.info[["geneSets"]][["rpls"]] <- v.ribo.genes <- grep("^rp[ls]", rownames(sc@assays$RNA@data), value=T, ignore.case=T)
l.info[["geneSets"]][["rRNA"]] <- v.rRNA.genes <- grep("rRNA", rownames(sc@assays$RNA@data), value=T, ignore.case=T)
l.info[["geneSets"]][["tRNA"]] <- v.tRNA.genes <- grep("tRNA", rownames(sc@assays$RNA@data), value=T, ignore.case=T)
l.info[["geneSets"]][["ERCC"]] <- v.ercc.genes <- grep("^ERCC", rownames(sc@assays$RNA@data), value=T, ignore.case=F)
l.info[["geneSets"]][["exclude"]] <- v.genes.exclude <- mixedsort(unique(c(v.mt.genes, v.ribo.genes, v.rRNA.genes, v.tRNA.genes, v.ercc.genes, "EGFP")))

cat( "length( v.genes.exclude ) = ", length( v.genes.exclude ), "\n" )
#length( v.genes.exclude ) =  472 

write(v.genes.exclude, file="~/Dropbox/singulomics/github_rda/GSE157504_RAW/gene.exclude.txt")

sc@meta.data %>% mutate(experiment = sprintf("%s_%s", experimental_design, time)) -> sc@meta.data

v.experiments <- unique(sc@meta.data$experiment)
cat( "v.experiments:", v.experiments, "\n" );

for (experiment in v.experiments){
  l.info[["nCells_condition"]][[experiment]] <- sum(sc@meta.data$experiment == experiment)
}
print( l.info )

#Data integration
l.info$geneSets$exclude -> v.var.genes.exclude

sc = PercentageFeatureSet(sc, pattern = "^mt:", col.name = "percent.mito")

sc_splited <- SplitObject(sc, split.by = "experiment")

for (i in 1:length(sc_splited)){
  sc_splited[[i]] <- SCTransform(sc_splited[[i]], 
                                 vars.to.regress = c("date", "nCount_RNA", "nFeature_RNA", "percent.mito"),
                                 assay="RNA", verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 3000)
  sc_splited[[i]]@assays$SCT@var.features <- setdiff(sc_splited[[i]]@assays$SCT@var.features, v.var.genes.exclude)
}

v.common <- rownames(sc_splited[[1]]$SCT@data)
for (i in 1:length(sc_splited)){
  if (i == 1){
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(v.common))) 
  }else{
    v.common <- intersect(v.common, rownames(sc_splited[[i]]$SCT@data))
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(rownames(sc_splited[[i]]$SCT@data))))
  }
} 
length(v.common)
cat("length(v.common)=", length(v.common), "\n")

v.inte.features <- vector()
if (1){                                     # find variable features/genes common to all
  v.var.common <- sc_splited[[1]]@assays$SCT@var.features
  print(sprintf("condition %s number of SCT var.features: %s", 1, length(v.var.common)))
  for (i in 2:length(sc_splited)){
    v.var.common <- c(v.var.common, sc_splited[[i]]@assays$SCT@var.features)
    print(sprintf("condition %s number of SCT var.features: %s", i, length(sc_splited[[i]]@assays$SCT@var.features)))
  }
  length(unique(v.var.common))
  v.var.common <- names(which(table(v.var.common) == length(sc_splited)))
  # v.var.common <- names( which( table( v.var.common ) == 8 ) ); # Change this just for integration LD and DD
  
  cat("length(v.var.common)=", length(v.var.common), "\n")
  v.var.common <- setdiff(v.var.common, v.var.genes.exclude)
  length(v.var.common)
  print("After excluding v.var.genes.exclude:")
  cat("length(v.var.common)=", length(v.var.common), "\n")
  v.inte.features <- v.var.common
} else {
  v.inte.features <- SelectIntegrationFeatures(object.list = sc_splited, nfeatures = 3000)
}

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
k.weight.use = 100
sc_splited.anchors <- FindIntegrationAnchors(object.list = sc_splited, 
                                             dims = 1:nDims, 
                                             assay=rep( "SCT", length(sc_splited)),
                                             anchor.features=v.inte.features, 
                                             k.anchor=k.anchor.use,
                                             k.filter=k.filter.use,
                                             k.score=k.score.use, verbose=T )

sc_integrated <- IntegrateData(anchorset = sc_splited.anchors, dims = 1:nDims, 
                               features.to.integrate=v.common, k.weight = k.weight.use, verbose = T)
VariableFeatures(sc_integrated, assay="integrated") <- v.inte.features

DefaultAssay(sc_integrated) <- "integrated"

sc_integrated <- ScaleData(sc_integrated, verbose=TRUE, assay="integrated", features=rownames(sc_integrated$integrated@data))
sc_integrated <- RunPCA(sc_integrated, assay="integrated", npcs = nDims, verbose = FALSE, reduction.name = "rna_integrated_pca", 
                        reduction.key = "rnaintegratedpc_"); cat( "RunPCA done\n" )
sc_integrated <- RunTSNE(sc_integrated, reduction = "rna_integrated_pca",
                         reduction.name = "rna_integrated_tsne", 
                         reduction.key = "rna_integrated_tSNE_", 
                         dims = 1:nDims ); cat( "RunTSNE done\n" )
sc_integrated <- FindNeighbors(sc_integrated, assay="integrated", reduction = "rna_integrated_pca", dims = 1:nDims, force.recalc=T)
sc_integrated <- FindClusters(sc_integrated, assay="integrated", resolution = 1.0 )

rm(sc_splited.anchors, sc_splited)
gc()

DimPlot(sc_integrated, group.by = "seurat_clusters", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq") -> p1
DimPlot(sc_integrated, group.by = "experiment", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq") -> p2
p1|p2

sc_integrated$seurat_clusters -> sc_integrated$integrated_cluster

#Cluster annotation 
DefaultAssay(sc_integrated) <- "SCT"
DimPlot(sc_integrated, group.by = "seurat_clusters", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq") -> p1


#Use new annotation
read.csv("~/Dropbox/singulomics/github_rda/clock_neurons_annotation.csv", head = T) -> df_anno

#filter DD cells
sc_integrated_tmp = sc_integrated
sc_integrated_tmp@meta.data %>% filter(grepl("LD", experiment)) %>% .$experiment %>% table()
sc_integrated_tmp@meta.data %>% filter(grepl("LD", experiment)) %>% rownames() %>% {sc_integrated_tmp[,.]} -> sc_integrated_tmp
sc_integrated_tmp@meta.data -> meta.data
sc_integrated_tmp@meta.data %>% rownames_to_column("X") %>% 
  {
    sc_meta = .
    df_anno %>% dplyr::mutate(X = sprintf("X%s", X) %>% toupper()) %>% 
      distinct() -> df_anno
    
    intersect(sc_meta$X, df_anno$X) %>% length() %>% print()
    
    left_join(sc_meta, df_anno, by = "X") %>% 
      column_to_rownames("X")
  } -> sc_integrated_tmp@meta.data
sc_integrated_tmp@meta.data %>% dplyr::mutate(celltype = gsub("\\d+:(.+)", "\\1", Idents)) -> sc_integrated_tmp@meta.data
DimPlot(sc_integrated_tmp, group.by = "Idents", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq")
DimPlot(sc_integrated_tmp, group.by = "celltype", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq")
rm(sc_integrated_tmp)

df_anno$X %>% sprintf("X%s", .) %>% toupper() -> cell_sele

sc_integrated[, cell_sele] -> sc_integrated_1
DefaultAssay(sc_integrated_1) <- "integrated"
sc_integrated_1 <- FindNeighbors(sc_integrated_1, assay="integrated", reduction = "rna_integrated_pca", dims = 1:nDims, force.recalc=T)
sc_integrated_1 <- FindClusters(sc_integrated_1, resolution = 1.0, graph.name = "integrated_snn")

list(
sc_integrated_1@meta.data %>% rownames_to_column("cellname"), 
df_anno %>% mutate(cellname = toupper(sprintf("X%s", X)), .after = 1) %>% 
  .[,c(2,5,6)]
) %>% purrr::reduce(., left_join, by = "cellname") %>% 
  column_to_rownames("cellname") -> sc_integrated_1@meta.data

DimPlot(sc_integrated_1, group.by = "seurat_clusters", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + ggtitle("scRNA-seq")
DimPlot(sc_integrated_1, group.by = "Idents", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + 
  theme(legend.position = "none") -> p1

saveRDS(sc_integrated_1, "~/Dropbox/singulomics/github_rda/Drosophia_scRNAseq_sc_obj.rds")
rm(list=ls())

library(tidyverse)
library(Seurat)

readRDS("~/Dropbox/singulomics/github_rda/Drosophia_scRNAseq_sc_obj.rds") -> sc

sc@meta.data %>% dplyr::mutate(celltype = gsub(".+:(.+)", "\\1", Idents)) -> sc@meta.data
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "rna_integrated_tsne", repel = TRUE) + 
  theme(legend.position = "none") -> p_tsne
p_tsne + theme_classic() + theme(legend.position = "none") + 
  ggtitle(NULL) -> p_tsne #Supp Fig 16A

sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) %>% 
  as.matrix() -> RNA_norm

sc$Idents %>% unique() %>%
  gtools::mixedsort() %>% 
  "names<-"(., .) %>% 
  map(function(x){
    print(x)
    sc@meta.data %>% 
      rownames_to_column("cellname") %>%
      filter(Idents == x) -> meta_
    meta_$experiment %>% unique() %>% 
      gtools::mixedsort() %>% "names<-"(.,.) -> exp_
    exp_ %>%
      map(function(x){
        meta_ %>% filter(experiment == x) %>% .$cellname -> cell_sele
        ZT_ = gsub(".+(ZT.+)", "\\1_REP%s", x) %>% {sprintf(., 1:length(cell_sele))}
        length(cell_sele)
        RNA_norm[,cell_sele] %>%
          as.data.frame() %>%
          "colnames<-"(., ZT_) -> RNA_norm
      }) %>% purrr::reduce(., cbind) -> RNA_norm
    RNA_norm %>% rownames_to_column("Gene") -> RNA_norm
  }) ->RNA_norm_list

##Circadian genes detection (mean and non zero proportion)
RNA_norm_list %>% 
#  .[1] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    cluster_ = y
    cluster_ %>% gsub("(\\d+):.+", "\\1", .) -> cluster_
    df_ = x %>% column_to_rownames("Gene")
    c("ZT02", "ZT06", "ZT10", "ZT14", "ZT18", "ZT22") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::select(contains(ZT_)) -> df_
        df_ %>% rowMeans() %>% 
          as.data.frame() %>% 
          "colnames<-"(., sprintf("%s_REP%s", ZT_, cluster_)) -> df_
      }) %>% purrr::reduce(., cbind) -> df_
  }) %>% purrr::reduce(., cbind) %>% 
  {
    df_ = .
    colnames(df_) %>% gtools::mixedsort() -> colnames_
    df_[, colnames_] -> df_
    df_ %>% rownames_to_column("Gene")
  } -> df_cyclic_neuron_mean

names(RNA_norm_list) %>% 
  {.[grepl("DN1p", .)]} %>% 
  map(function(x){
    cluster_ = x
    RNA_norm_list[[cluster_]] %>% column_to_rownames("Gene") -> df_
    c("ZT02", "ZT06", "ZT10", "ZT14", "ZT18", "ZT22") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::select(contains(ZT_)) -> df_
        df_ %>% rowMeans() %>% 
          as.data.frame() %>% 
          "colnames<-"(., sprintf("%s_REP%s", ZT_, cluster_)) -> df_
      }) %>% purrr::reduce(., cbind) -> df_
  }) %>% purrr::reduce(., cbind) %>% 
  {
    df_ = .
    colnames(df_) %>% gtools::mixedsort() -> colnames_
    df_[, colnames_] -> df_
    df_ %>% rownames_to_column("Gene")
  } %>% "colnames<-"(., gsub("(.+):.+", "\\1", colnames(.))) -> df_cyclic_neuron_mean_DN1p

RNA_norm_list %>% 
#    .[1] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    cluster_ = y
    cluster_ %>% gsub("(\\d+):.+", "\\1", .) -> cluster_
    df_ = x %>% column_to_rownames("Gene")
    c("ZT02", "ZT06", "ZT10", "ZT14", "ZT18", "ZT22") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::select(contains(ZT_)) -> df_
        n_total_cell = ncol(df_)
        n_non_zero_cell = apply(df_, 1, function(x){sum(x!=0)})
        n_non_zero_prop = n_non_zero_cell/n_total_cell
        n_non_zero_prop %>% as.data.frame() %>% 
          "colnames<-"(., sprintf("%s_REP%s", ZT_, cluster_)) -> df_
      }) %>% purrr::reduce(., cbind) -> df_
  }) %>% purrr::reduce(., cbind) %>% 
  {
    df_ = .
    colnames(df_) %>% gtools::mixedsort() -> colnames_
    df_[, colnames_] -> df_
    df_ %>% rownames_to_column("Gene")
  } -> df_cyclic_neuron_non_zero_prop

names(RNA_norm_list) %>% 
  {.[grepl("DN1p", .)]} %>% 
  map(function(x){
    cluster_ = x
    RNA_norm_list[[cluster_]] %>% column_to_rownames("Gene") -> df_
    c("ZT02", "ZT06", "ZT10", "ZT14", "ZT18", "ZT22") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::select(contains(ZT_)) -> df_
        df_ %>% 
#          rowMeans() %>% 
          apply(., 1, function(x){sum(x!=0)/length(x)}) %>% 
          as.data.frame() %>% 
          "colnames<-"(., sprintf("%s_REP%s", ZT_, cluster_)) -> df_
      }) %>% purrr::reduce(., cbind) -> df_
  }) %>% purrr::reduce(., cbind) %>% 
  {
    df_ = .
    colnames(df_) %>% gtools::mixedsort() -> colnames_
    df_[, colnames_] -> df_
    df_ %>% rownames_to_column("Gene")
  } %>% "colnames<-"(., gsub("(.+):.+", "\\1", colnames(.))) -> df_cyclic_neuron_non_zero_prop_DN1p


source("~/Dropbox/singulomics/github/Calculate_HMP.R")

list(df_ = df_cyclic_neuron_mean) %>% 
  map(function(x){
    df_ = x
    df_$Gene -> gene_
    df_ %>% dplyr::select(-Gene) -> exp_matrix_
    colnames(exp_matrix_) %>% gsub("ZT(.+?)_.+", "\\1", .) %>% as.numeric() -> timepoints_
    cyclic_HMP(exp_matrix=exp_matrix_, gene = gene_, timepoints = timepoints_) -> res_
    res_ %>% dplyr::select(!matches("LS|meta2d", ignore.case = FALSE)) %>% recal_cauchy_p_and_hmp(.) -> res_
  }) -> res_mean
res_mean$df_ %>% dplyr::mutate(F24_phase = (HR_phi/(2*pi))*24) -> res_mean$df_

list(df_ = df_cyclic_neuron_non_zero_prop) %>% 
  map(function(x){
    df_ = x
    df_$Gene -> gene_
    df_ %>% dplyr::select(-Gene) -> exp_matrix_
    colnames(exp_matrix_) %>% gsub("ZT(.+?)_.+", "\\1", .) %>% as.numeric() -> timepoints_
    cyclic_HMP(exp_matrix=exp_matrix_, gene = gene_, timepoints = timepoints_) -> res_
    res_ %>% dplyr::select(!matches("LS|meta2d", ignore.case = FALSE)) %>% recal_cauchy_p_and_hmp(.) -> res_
  }) -> res_non_zero_prop
res_non_zero_prop$df_ %>% dplyr::mutate(F24_phase = (HR_phi/(2*pi))*24) -> res_non_zero_prop$df_

list(df_ = df_cyclic_neuron_mean_DN1p) %>% 
  map(function(x){
    df_ = x
    df_$Gene -> gene_
    df_ %>% dplyr::select(-Gene) -> exp_matrix_
    colnames(exp_matrix_) %>% gsub("ZT(.+?)_.+", "\\1", .) %>% as.numeric() -> timepoints_
    cyclic_HMP(exp_matrix=exp_matrix_, gene = gene_, timepoints = timepoints_) -> res_
    res_ %>% dplyr::select(!matches("LS|meta2d", ignore.case = FALSE)) %>% recal_cauchy_p_and_hmp(.) -> res_
  }) -> res_mean_DN1p
res_mean_DN1p$df_ %>% dplyr::mutate(F24_phase = (HR_phi/(2*pi))*24) -> res_mean_DN1p$df_

list(df_ = df_cyclic_neuron_non_zero_prop) %>% 
  map(function(x){
    df_ = x
    df_$Gene -> gene_
    df_ %>% dplyr::select(-Gene) -> exp_matrix_
    colnames(exp_matrix_) %>% gsub("ZT(.+?)_.+", "\\1", .) %>% as.numeric() -> timepoints_
    cyclic_HMP(exp_matrix=exp_matrix_, gene = gene_, timepoints = timepoints_) -> res_
    res_ %>% dplyr::select(!matches("LS|meta2d", ignore.case = FALSE)) %>% recal_cauchy_p_and_hmp(.) -> res_
  }) -> res_non_zero_prop_DN1p
res_non_zero_prop_DN1p$df_ %>% dplyr::mutate(F24_phase = (HR_phi/(2*pi))*24) -> res_non_zero_prop_DN1p$df_

#Plot circadian genes
library(patchwork)
list(all_cyclic_neuron = df_cyclic_neuron_mean, DN1p_cyclic_neuron = df_cyclic_neuron_mean_DN1p) -> df_cycli_neuron_mean
list(all_cyclic_neuron = df_cyclic_neuron_non_zero_prop, DN1p_cyclic_neuron = df_cyclic_neuron_non_zero_prop_DN1p) -> df_cyclic_neuron_non_zero_prop
res_mean = list(all_cyclic_neuron = res_mean, DN1p_cyclic_neuron = res_mean_DN1p)
res_non_zero_prop = list(all_cyclic_neuron = res_non_zero_prop, DN1p_cyclic_neuron = res_non_zero_prop_DN1p)

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    c("tim", "Clk", "per", "cry", "DIP-beta", "vri") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        gene_ = x
        
        df_cyclic_neuron_mean[[neuron_type_]] %>% 
          filter(Gene == gene_) %>% 
          pivot_longer(cols = -Gene, names_to = "timepoint", values_to = "expression") %>% 
          mutate(cluster = gsub(".+_REP(\\d+)", "\\1", timepoint)) %>% 
          mutate(timepoint = gsub("ZT(\\d+)_.+", "\\1", timepoint) %>% as.numeric()) -> df_mean
        
        df_cyclic_neuron_non_zero_prop[[neuron_type_]] %>% 
          filter(Gene == gene_) %>% 
          pivot_longer(cols = -Gene, names_to = "timepoint", values_to = "expression") %>% 
          mutate(cluster = gsub(".+_REP(\\d+)", "\\1", timepoint)) %>% 
          mutate(timepoint = gsub("ZT(\\d+)_.+", "\\1", timepoint) %>% as.numeric()) -> df_non_zero_prop
        
        mean_p_adj = res_mean[[neuron_type_]]$df_ %>% filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
        df_mean %>% 
          ggplot(aes(x = timepoint, y = expression, group = timepoint)) + 
          geom_boxplot(fill = "lightblue") + 
          theme_classic() + 
          ylab("Gene expression") + 
          xlab("ZT") + 
          scale_x_continuous(breaks = seq(2,22,4)) + 
          ggtitle(sprintf("%s\np.adj=%s", gene_, mean_p_adj)) -> p_mean
        
        non_zero_proportion_p_adj = res_non_zero_prop[[neuron_type_]]$df_ %>% filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
        df_non_zero_prop %>% 
          ggplot(aes(x = timepoint, y = expression, group = timepoint)) + 
          geom_boxplot(fill = "lightgreen") + 
          scale_y_continuous(limits = c(0,1)) + 
          theme_classic() + 
          ylab("Non zero proportion") + 
          xlab("ZT") + 
          scale_x_continuous(breaks = seq(2,22,4)) + 
          ggtitle(sprintf("%s\np.adj=%s", gene_, non_zero_proportion_p_adj)) -> p_non_zero_prop
        
        wrap_plots(p_mean, p_non_zero_prop, ncol = 2) -> p
        
      }) -> p_list
  }) -> p_list

wrap_plots(p_list$all_cyclic_neuron, ncol = 1)
wrap_plots(p_list$DN1p_cyclic_neuron, ncol = 1)
wrap_plots(p_list, ncol = 2) -> p_cyclic_genes
wrap_plots(p_tsne, p_cyclic_genes, ncol = 2, widths = c(0.8, 1)) -> p_r1

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    c("tim", "Clk", "per", "cry", "DIP-beta", "vri") %>% 
      "names<-"(.,.) %>% 
      #  .[1] %>% 
      map(function(x){
        gene_ = x
        gene_expr_mat = p_list[[neuron_type_]][[gene_]][[1]][["data"]] %>% group_by(Gene, timepoint) %>% 
          summarise(mean_expression = mean(expression), sd = sd(expression)) %>% 
          dplyr::mutate(group = "Mean")
        non_zero_mat = p_list[[neuron_type_]][[gene_]][[2]][["data"]] %>% group_by(Gene, timepoint) %>%
          summarise(mean_expression = mean(expression), sd = sd(expression)) %>% 
          dplyr::mutate(group = "Nonzero")
        scale_factor = max(gene_expr_mat$mean_expression)/max(non_zero_mat$mean_expression)
        non_zero_mat$mean_expression = non_zero_mat$mean_expression*scale_factor
        non_zero_mat$sd = non_zero_mat$sd*scale_factor
        rbind(gene_expr_mat, non_zero_mat) -> df_
        df_ %>% 
          ggplot(aes(x = timepoint, y = mean_expression, group = group, fill = group, color = group)) + 
          geom_line() + 
          geom_ribbon(aes(ymin = mean_expression - sd, ymax = mean_expression + sd), alpha = 0.2, color = NA) +
          scale_y_continuous(name = "Mean", sec.axis = sec_axis(~./scale_factor, name = "Nonzero_Prop")) + 
          scale_color_manual(values = c(Mean="#56B4E9", Nonzero="#009E73")) +
          scale_fill_manual(values = c(Mean="#56B4E9", Nonzero="#009E73")) + 
          xlab("ZT") -> p
        
        mean_p_adj = res_mean[[neuron_type_]]$df_ %>% filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
        non_zero_proportion_p_adj = res_non_zero_prop[[neuron_type_]]$df_ %>% filter(Gene == gene_) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
        p + ggtitle(sprintf("%s\n(Mean) p.adj=%s\n(Nonzero) p.adj=%s", gene_, mean_p_adj, non_zero_proportion_p_adj))
      }) -> p_list_1
  }) -> p_list_1 #Supp Fig 16B


c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    list(
      Mean = res_mean[[neuron_type_]]$df_ %>% filter(cauchy_BH.Q < 0.05) %>% .$Gene, 
      `Nonzero Prop.` = res_non_zero_prop[[neuron_type_]]$df_ %>% filter(cauchy_BH.Q < 0.05) %>% .$Gene
    ) %>% ggvenn::ggvenn(., fill_color = c("lightblue", "lightgreen")) -> p_venn
  }) -> p_venn
p_venn$all_cyclic_neuron #Supp Fig 16C
#p_venn$DN1p_cyclic_neuron


## P adjust correlation (Mean vs Non zero prop.)
c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    list(
      Mean = res_mean[[neuron_type_]]$df_$cauchy_BH.Q,
      Non_zero = res_non_zero_prop[[neuron_type_]]$df_$cauchy_BH.Q
    ) %>% as.data.frame() %>% 
      -log(., 10) %>% 
      drop_na() %>% 
      {
        df_ = .
        x_intercept = median(df_$Mean)
        y_intercept = median(df_$Non_zero)
        cor_ = cor(df_$Mean, df_$Non_zero, method = "pearson") %>% round(., 2)    
        
        df_ %>% 
          ggplot(aes(x = Mean, y = Non_zero)) + 
          #      geom_point(size = 0.3, alpha = 0.5) + 
          geom_hex(bins = 100) + 
          scale_fill_gradient(trans = 'log', name = "Count (log scale)") +
          geom_vline(xintercept = x_intercept, linetype = "dashed", color = "red") +
          geom_hline(yintercept = y_intercept, linetype = "dashed", color = "red") + 
          theme_classic() + 
          annotate(geom = "text", x = 2.5, y = 11.5, label = sprintf("r = %s", cor_)) + 
          xlab("-log10(p.adj):mean") + 
          ylab("-log10(p.adj):non zero prop.")
      } -> p_val_correlation
  }) -> p_val_correlation
p_val_correlation$all_cyclic_neuron #Fig 3E
#p_val_correlation$DN1p_cyclic_neuron


# Phase correlation (Mean vs Non zero prop.)
c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    
    list(
      Mean = res_mean[[neuron_type_]]$df_$F24_phase,
      Non_zero = res_non_zero_prop[[neuron_type_]]$df_$F24_phase
    ) %>% as.data.frame() %>% 
      drop_na() %>% 
      dplyr::mutate(Non_zero = case_when(
#        Mean - Non_zero > 12 ~ Non_zero+24,
#        Mean - Non_zero < -12 ~ Non_zero-24,
        Non_zero - Mean > 12 ~ Non_zero - 24,
        TRUE ~ Non_zero
      )) %>% 
      dplyr::mutate(Mean = case_when(
        Non_zero - Mean < -12 ~ Mean - 24, 
        TRUE ~ Mean
      )) %>% 
      {
        df_ = .
        print(nrow(df_))
        cor_ = cor(df_$Mean, df_$Non_zero, method = "pearson")
        df_ %>%
          ggplot(aes(x = Mean, y = Non_zero)) + 
          geom_hex(bins = 100) + 
          scale_fill_gradient(trans = 'log', name = "Count (log scale)") +
          geom_abline(intercept = 0, slope = 1, color = "red") + 
          annotate("text", x = 3, y = 23, label = sprintf("cor = %.2f", cor_)) + 
          ylab("Adjusted phase: nonzero prop.") + 
          xlab("Adjusted phase: mean") + 
          theme_classic() + 
          scale_x_continuous(breaks = seq(-10,20,5)) + 
          scale_y_continuous(breaks = seq(-10,20,5))
      } -> p_phase_correlation
  }) -> p_phase_correlation
p_phase_correlation$all_cyclic_neuron #Fig 3E
#p_phase_correlation$DN1p_cyclic_neuron

# Amplitude correlation (Mean vs Non zero prop.)
c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    list(
      Mean = res_mean[[neuron_type_]]$df_$MetaCycle_JTK_amplitude,
      Non_zero = res_non_zero_prop[[neuron_type_]]$df_$MetaCycle_JTK_amplitude
    ) %>% as.data.frame() %>% 
      drop_na() %>% 
      log(.) %>% 
      filter(is.finite(Mean)&is.finite(Non_zero)) %>% 
      {
        df_ = .
        #    df_ %>% dplyr::filter(Non_zero > -20) -> df_
        df_ %>% dplyr::filter(Non_zero > quantile(Non_zero, 0.01)) -> df_
        #    df_ %>% dplyr::filter(Non_zero < quantile(Non_zero, 0.95)) -> df_
        print(nrow(df_))
        cor_ = cor(df_$Mean, df_$Non_zero, method = "pearson")
        df_ %>%
          ggplot(aes(x = Mean, y = Non_zero)) + 
          geom_hex(bins = 100) + 
          #      scale_fill_gradient(trans = 'log', low = "lightblue", high = "darkblue", name = "Count (log scale)") +
          #      geom_abline(intercept = 0, slope = 1, color = "red") + 
          annotate("text", x = -8, y = -2, label = sprintf("cor = %.2f", cor_)) + 
          ylab("log(amplitude): nonzero prop.") + 
          xlab("log(amplitude): mean") + 
          theme_classic() 
        #      scale_x_continuous(limits = c(-10, 10))
      } -> p_amp_correlation
  }) -> p_amp_correlation
p_amp_correlation$all_cyclic_neuron #Fig 3E
#p_amp_correlation$DN1p_cyclic_neuron

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    list(
      Mean = res_mean[[neuron_type_]]$df_ %>% filter(cauchy_BH.Q < 0.05) %>% .$Gene, 
      Nonzero = res_non_zero_prop[[neuron_type_]]$df_ %>% filter(cauchy_BH.Q < 0.05) %>% .$Gene
    ) %>% ggvenn::list_to_data_frame() %>% 
      list(df_ = .) %>% 
      map(function(x){
        df_ = x
        df_ %>% filter(Mean&Nonzero) %>% .$key -> Overlapped
        df_ %>% filter(Mean&!Nonzero) %>% .$key -> Mean
        df_ %>% filter(!Mean&Nonzero) %>% .$key -> Nonzero
        list(Overlapped = Overlapped, Mean = Mean, Nonzero = Nonzero) -> list_
      }) %>% .[[1]] -> gene_set
  }) -> gene_set

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    gene_set[[neuron_type_]] %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        group_ = y
        gene_ = x
        df_cyclic_neuron_mean[[neuron_type_]] %>% column_to_rownames("Gene") %>% .[gene_, ] %>% as.vector() %>% 
          purrr::reduce(., c) %>% 
          as.data.frame() %>% 
          "colnames<-"(., c("expression")) %>% 
          mutate(group = group_)
      }) %>% do.call(rbind, .) -> df_
  }) -> df_

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    df_[[neuron_type_]] %>% 
      dplyr::filter(group %in% c("Mean", "Nonzero")) -> df_
    df_ %>% dplyr::filter(group == "Mean") %>% .$expression %>% quantile(., 0.75) %>% .[[1]] %>% {. + 1} %>% as.integer() -> up_break
    df_ %>% dplyr::filter(group == "Mean") %>% .$expression %>% quantile(., 0.99)%>% .[[1]] %>% {. + 1} %>% as.integer() -> up_break_max
    print(c(up_break, up_break_max))
#      mutate(group = factor(group, levels = c("Mean", "Overlapped", "Nonzero"))) %>% 
    df_ %>% 
      mutate(group = factor(group, levels = c("Mean", "Nonzero"))) %>% 
      ggplot(aes(x = group, y = expression, fill = group)) + 
      geom_boxplot(width = 0.45) + 
      scale_y_break(c(up_break, up_break_max), scales = 0.5) + 
      theme_classic() + 
      scale_fill_manual(values = c("Mean" = "#56B4E9", "Nonzero" = "#009E73")) +
#      ggpubr::stat_compare_means(comparisons = list(c("Mean", "Overlapped"), 
#                                                    c("Overlapped", "Nonzero"), 
#                                                    c("Mean", "Nonzero")), method = "t.test") + 
      theme(legend.position = "none") + 
      xlab(NULL) + 
      ylab("Normalized expression") -> p_1
  })-> p_1
#p_1$all_cyclic_neuron+p_1$DN1p_cyclic_neuron
p_1$all_cyclic_neuron #Supp Fig 16D

c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    gene_set[[neuron_type_]] %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        group_ = y
        gene_ = x
        df_cyclic_neuron_non_zero_prop[[neuron_type_]] %>% column_to_rownames("Gene") %>% .[gene_, ] %>% as.vector() %>% 
          purrr::reduce(., c) %>% 
          as.data.frame() %>% 
          "colnames<-"(., c("expression")) %>% 
          mutate(group = group_)
      }) %>% do.call(rbind, .) %>% 
      filter(is.finite(expression)) -> df_
  }) -> df_


c("all_cyclic_neuron", "DN1p_cyclic_neuron") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    neuron_type_ = x
    df_[[neuron_type_]] %>% 
      dplyr::filter(group %in% c("Mean", "Nonzero")) -> df_
    df_ %>% 
      mutate(group = factor(group, levels = c("Mean", "Nonzero"))) %>% 
      ggplot(aes(x = group, y = expression, fill = group)) + 
      geom_boxplot(width = 0.5) + 
      scale_fill_manual(values = c("Mean" = "#56B4E9", "Nonzero" = "#009E73")) +
      theme_classic() + 
      theme(legend.position = "none") + 
      xlab(NULL) + 
      ylab("Non zero proportion") -> p_2
  }) -> p_2
#p_2$all_cyclic_neuron+p_2$DN1p_cyclic_neuron
p_2$all_cyclic_neuron #Supp Fig 16D