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
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Mean.csv", header = T, stringsAsFactors = F) -> res_RNA_Mean
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_NonZeroProp.csv", header = T, stringsAsFactors = F) -> res_RNA_NonZeroProp
list(Mean = res_RNA_Mean, Nonzero_Prop = res_RNA_NonZeroProp) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    type_ = y
    if (type_ == "Mean"){
      x %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene -> genes
    }else{
      x %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene -> genes
    }
  }) -> genes_list
ggvenn::ggvenn(genes_list) -> p_venn #Fig_3C
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