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
ggsave(filename='bursting/bursting_mean.pdf', width=13, height=3.5, plot=p)
####

# 15. Upset plot of Nonzero_mean, gini, Nonzero_Prop, cv, and mean ----
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Mean.csv", header = T, stringsAsFactors = F) -> res_RNA_Mean
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_NonZeroProp.csv", header = T, stringsAsFactors = F) -> res_RNA_NonZeroProp
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_Nonzero_Mean.csv", header = T, stringsAsFactors = F) -> res_RNA_Nonzero_Mean
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_gini.csv", header = T, stringsAsFactors = F) -> res_RNA_gini
read.csv("~/Dropbox/singulomics/github_rda/output/Beyond_Mean/res_RNA_cv.csv", header = T, stringsAsFactors = F) -> res_RNA_cv
list(Mean = res_RNA_Mean, 
     Nonzero_Prop = res_RNA_NonZeroProp,
     Nonzero_mean = res_RNA_Nonzero_Mean
) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    type_ = y
    if (type_ == "Mean"){
      x %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene -> genes
    }else if (type_ == "Nonzero_Prop") {
      x %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene -> genes
    }else{
      x %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene -> genes 
    }
  }) -> genes_list
ggvenn::ggvenn(genes_list) -> p_venn

library(ggupset)
list(
  Mean = res_RNA_Mean %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene, 
  Nonzero_Prop = res_RNA_NonZeroProp %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene,
  Nonzero_mean = res_RNA_Nonzero_Mean %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene,
  cv = res_RNA_cv %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.05) %>% .$Gene,
  gini = res_RNA_gini %>% dplyr::filter(MetaCycle_JTK_BH.Q < 0.01) %>% .$Gene
) %>% 
  ggvenn::list_to_data_frame() %>% 
  dplyr::mutate(across(!matches("key"), function(x){x %>% as.integer()})) %>% 
  as.data.frame() %>% 
  UpSetR::upset(sets = c("Mean", "cv", "Nonzero_Prop", "gini", "Nonzero_mean"), keep.order = T, order.by = "freq") -> p_upset #Supp_Fig_14B

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
pairs(toplot, lower.panel = function(x,y){smoothScatter(x,y,add=TRUE)}, upper.panel = panel.cor, main='Pairwise Spearman correlation')
dev.off()

####