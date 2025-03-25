setwd("~/Library/CloudStorage/Dropbox/singulomics/github")

library(tidyverse)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(patchwork)

# 1. Read sc of hepatocytes with infered pseudo time ----

load("~/Dropbox/singulomics/github_rda/trajectory_analysis.rda")
rm(ATres, heatdata, p_list, RNA_norm, sce);gc()

sc$pseudotime=1-sc$pseudotime

table(sc$ZT)
table(sc$group)
sc$group=droplevels(sc$group)

sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))

# Generate normalized RNA matrix
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) %>% 
  as.matrix() -> RNA_norm

# Generate normalized ATAC activity matrix
sc@assays[["ATAC"]]@fragments[[1]]@path <- "~/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
DefaultAssay(sc) <- "ATAC"
gene.activities <- GeneActivity(sc, extend.upstream = 2000, biotypes = NULL)
colSums(gene.activities) %>% median() -> gene_activity_median
gene.activities %>%
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*gene_activity_median})) %>% 
  as.matrix() -> gene_activity_norm
rm(gene_activity_median, gene.activities);gc()

# 2. Genrate normalized RNA matrix and ATAC activity matrix (genes x metacells) with different resolutions ----
readRDS("~/Dropbox/singulomics/github_rda/output/Celltype_specific/intersect_genes_high_exp.rds") -> highly_expressed_genes
c("rna", "gene_activity") %>% 
  "names<-"(.,.) %>% 
  .[2] %>%
  map(function(x){
    seq_ = x
    seq(0.5, 20, 0.5) %>%
      "names<-"(., sprintf("res_%s", .)) %>% 
      .[29:40] %>% 
      map(function(x){
        res_ = x
        print(res_)
        sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                           resolution=res_)
        
        pseudotime.median=aggregate(sc$pseudotime, by=list(sc$seurat_clusters), FUN=median)
        pseudotime.median=cbind(pseudotime.median,rank(pseudotime.median[,2])-1)
        new.cluster.id=rep(NA, ncol(sc))
        
        for(i in 1:nrow(pseudotime.median)){
          new.cluster.id[which(sc$seurat_clusters==pseudotime.median[i,1])]=pseudotime.median[i,3]
        }
        sc$seurat_clusters=factor(as.numeric(new.cluster.id))
        rm(pseudotime.median); rm(new.cluster.id)
        
        sc$cluster=sc$seurat_clusters
        sc$cluster %>% unique() %>% as.character() %>% as.numeric() %>% gtools::mixedsort() %>% 
          map(function(x){
            cluster_ = x
            print(cluster_)
            sprintf("ZT%s", c("02", "06", "10", "14", "18", "22")) %>%
              "names<-"(.,.) %>% 
              map(function(x){
                ZT_ = x
                print(ZT_)
                sc@meta.data %>% filter(cluster == cluster_, ZT == ZT_) -> meta_data
                meta_data$pseudotime %>% median() -> pseudotime_median
                rownames(meta_data) -> cellnames_
                if(length(cellnames_) == 0){
                  df_ = data.frame(Gene = character(), Mean_expression = numeric(), 
                                   ZT = character(), cluster = numeric(), pseudotime = numeric())
                }else{
                  if(length(cellnames_) > 1){
                    
                    if(seq_ == "rna"){
                      RNA_norm[highly_expressed_genes, cellnames_] %>% rowMeans() -> mean_                       
                    }else{
                      gene_activity_norm[highly_expressed_genes, cellnames_] %>% rowMeans() -> mean_
                    }

                  }else{
                    
                    if(seq_ == "rna"){
                      RNA_norm[highly_expressed_genes, cellnames_] -> mean_                                   
                    }else{
                      gene_activity_norm[highly_expressed_genes, cellnames_] -> mean_
                    }

                  }
                  mean_ %>% 
                    as.data.frame() %>% 
                    "colnames<-"(., "Mean_expression") %>% 
                    rownames_to_column("Gene") %>%
                    mutate(ZT = ZT_, cluster = cluster_, pseudotime = pseudotime_median) -> df_              
                }
                return(df_)
              }) %>% do.call(rbind, .)
          }) %>% do.call(rbind, .) %>% 
          "rownames<-"(., 1:nrow(.)) %>% 
          dplyr::arrange(Gene) -> df_
        
        i = 0
        highly_expressed_genes %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            i <<- i + 1
            print(i)
            gene_ = x
            print(x)
            df_ %>% filter(Gene == gene_) -> df_
          }) -> output.list
        rm(df_);gc()
        saveRDS(output.list, sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_%s_%s.rds", res_, seq_))
        #    return(output.rna.list)
      }) -> list_tmp
  }) -> list_tmp

rm(list_tmp)
gc()

# 3. Calculate spatial and circadian genes along the liver lobule ----
source('dynamic.R')
#list.output = list()
c("rna", "gene_activity") %>% 
#  .[1] %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    seq_ = x
    print(seq_)
    list.files("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data", pattern = sprintf("%s\\.rds", seq_), full.names = T) %>% 
      gtools::mixedsort() -> files_
    files_ %>% 
#      .[1] %>% 
      "names<-"(., gsub("/.+/(res_.+?)_.+", "\\1", .)) %>% 
#      .[4] %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        res_ = y
        file_ = x
        print(file_)
        readRDS(file_) -> list_
#        list.output[[seq_]][[res_]] <<- list_
        highly_expressed_genes %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            gene_ = x
            list_[[gene_]] -> df_
            y = df_$Mean_expression
            t = df_$ZT %>% gsub("ZT(.+)", "\\1", .) %>% as.numeric()
            p = df_$pseudotime
            
            lrt.mod=runModel(y,p,t,k = 2) # polynomial with 2 df
            lrt.mod=lrt(lrt.mod)
            lrt.mod$lrt.pval %>% as.data.frame() %>% t() %>% 
              as.data.frame() %>% 
              mutate(Gene = gene_, .before = 1) -> result_
          }) %>% do.call(rbind, .) -> result_
        
        ####Cauchy combination
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
        #####
        
        alpha = 0.01
        result_ %>% as_tibble() %>% column_to_rownames("Gene") -> result_
        Gene = rownames(result_)
        Circadian_p = result_$H1vsH0
        Transient_p = result_$H2vsH0
        Circadian_p.adj = p.adjust(Circadian_p, method = "BH")
        Transient_p.adj = p.adjust(Transient_p, method = "BH")
        data.frame(Gene, Circadian_p, Transient_p, Circadian_p.adj, Transient_p.adj, 
                   H3vsH2 = result_$H3vsH2, H3vsH1 = result_$H3vsH1) -> result_
        result_[is.na(result_)] = 1
        print(any(is.na(result_)))
        

        result_ %>% 
          dplyr::rowwise() %>% 
          mutate(Circadian_transient_p = case_when(
          Circadian_p.adj<=alpha & Transient_p.adj<=alpha ~ cauchy.test(c(Circadian_p, Transient_p)),
          Circadian_p.adj<=alpha & Transient_p.adj>alpha ~ H3vsH1,
          Circadian_p.adj>alpha & Transient_p.adj<=alpha ~ H3vsH2,
          TRUE ~ NA
        )) -> result_
#        head(result_)
        result_ %>% mutate(Circadian_transient_p.adj = p.adjust(Circadian_transient_p, method = "BH")) -> result_
        saveRDS(result_, sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/%s.%s.rds", res_, seq_)) 
        return(result_) 
      })
  }) -> list_p_value 

# 4. Resolution tuning ----
ggroc_list = list()
density_list = list()
c("rna", "gene_activity") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    seq_ = x
    print(seq_)
    sprintf("res_%s", seq(0.5,20,0.5)) %>% 
#      c("res_2.5") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        res_ = x
        print(res_)
        list.files("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results", pattern = sprintf("%s\\.%s\\.rds", res_, seq_), full.names = T) -> file_
        readRDS(file_) -> df_
        n_metacell = readRDS(sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/%s_rna.rds", res_))[[1]]$cluster %>% unique() %>% length()
        
        ###density plot###
        df_[,c("Gene", "Circadian_p.adj")] %>% mutate(resolution = gsub("res_(.+)", "\\1", res_) %>% as.numeric()) ->> density_list[[seq_]][[res_]]
        ###
        
        c(0.05, 0.01) %>% 
          "names<-"(., sprintf("p%s", .)) %>% 
          map(function(x){
            threshold_ = x
            threshold_name_ = sprintf("p_%s", threshold_)
            n_circadian = df_ %>% filter(Circadian_p.adj<=threshold_) %>% nrow()
            n_transient = df_ %>% filter(Transient_p.adj<=threshold_) %>% nrow()
            n_both = df_ %>% filter(Circadian_transient_p.adj<=threshold_) %>% nrow()
            data.frame(resolution = res_, n_circadian, n_transient, n_both, n_metacell, threshold = threshold_) -> df_1
            df_1 %>% mutate(resolution = gsub("res_(.+)", "\\1", resolution) %>% as.numeric()) -> df_1
            
            ##caluclate auc###
            pos.genes=readRDS('~/Dropbox/singulomics/rda/pos.genes.RDS')
            neg.genes=readRDS('~/Dropbox/singulomics/rda/neg.genes.RDS')
            
            df_ %>% "rownames<-"(., .$Gene) %>% .[pos.genes, c("Gene", "Circadian_p.adj")] %>% 
              mutate(group = "pos.genes") %>% drop_na() -> df_pos
            df_ %>% "rownames<-"(., .$Gene) %>% .[neg.genes, c("Gene", "Circadian_p.adj")] %>% 
              mutate(group = "neg.genes") %>% drop_na() -> df_neg
            print("Done")
            
            rbind(df_pos, df_neg) %>% 
              mutate(response = ifelse(group == "pos.genes", 0, 1)) -> df_roc
            pROC::roc(df_roc$response, df_roc$Circadian_p.adj) -> roc_obj
#            ggroc_list[[seq_]][[res_]] <<- roc_obj
            ggroc_list[[threshold_name_]][[seq_]][[res_]] <<- roc_obj
            auc_ = roc_obj$auc %>% as.numeric()
            df_1 %>% mutate(auc = auc_) -> df_1
            df_1 %>% mutate(assay = seq_) -> df_1
            ###End ###
            
            return(df_1)
          }) %>% do.call(rbind, .) -> df_
        gc()
        return(df_)
      }) %>% do.call(rbind, .) -> df_summary
  }) %>% do.call(rbind, .) -> df_summary
df_summary %>% dplyr::arrange(resolution) -> df_summary

ggroc_list$p_0.01 %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seq_ = y
    pROC::ggroc(x)$data %>%
      as.data.frame() %>% 
      mutate(res = gsub("res_(.+)", "\\1", name) %>% as.numeric()) %>% 
      ggplot(aes(x = 1-specificity, y = sensitivity, color = res, group = res)) + 
      geom_line() + 
      scale_color_gradient(low = "blue", high = "red") + 
      theme_classic() + 
      ggtitle(sprintf("%s: H1 model. 12435 genes in total", seq_))-> p_roc
  }) -> p_roc #Supp_Fig8B

c(0.05, 0.01) %>% 
  "names<-"(., sprintf("p_%s", .)) %>% 
  map(function(x){
    threshold_ = x
    df_summary %>% dplyr::filter(threshold == threshold_) -> df_summary
    df_summary[, c("assay", "n_metacell", "auc")] %>% 
      distinct() %>% 
      mutate(assay = factor(assay, levels = c("rna", "gene_activity"))) %>% 
      ggplot(aes(x = n_metacell, y = auc, color = assay)) + 
      geom_point() +
      geom_smooth(span = 0.15) + 
      scale_color_manual(values = c("blue", "red")) +
      geom_vline(xintercept = 25, linetype = "dashed") + 
      theme_classic() -> p_auc
  }) -> p_auc #Supp_Fig8C

df_summary %>% 
  filter(threshold == 0.01) %>% 
  mutate(assay = factor(assay, levels = c("rna", "gene_activity"))) %>% 
  ggplot(aes(x = n_metacell, y = n_circadian, color = assay)) + 
           geom_point() +
           geom_smooth() + 
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = 25, linetype = "dashed") + 
           theme_classic() -> p_n_circadian #Supp_Fig8D

df_summary[, c("resolution", "n_metacell")] %>% 
  distinct() %>%
  ggplot(aes(x = resolution, y = n_metacell)) +
  geom_point() +
  geom_smooth(color = "black") + 
  geom_vline(xintercept = 2.5, linetype = "dashed") +
  geom_hline(yintercept = 25, linetype = "dashed") +
  theme_classic() -> p_res_nmetacell #Supp_Fig8A

density_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seq_ = y
    x %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        res_ = y
        x %>% mutate(assay = seq_)
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>% 
  as_tibble() %>% 
  dplyr::mutate(assay = factor(assay, levels = c("rna", "gene_activity"))) %>% 
  ggplot(aes(x = Circadian_p.adj, color = resolution, group = resolution)) + 
  scale_color_gradient(low = "blue", high = "red") +
  geom_density() + 
  theme_classic() + 
  facet_wrap(~assay, scales = "free", nrow = 1) -> p_density #Supp_Fig8E

library(patchwork)
p_roc + p_auc + p_n_circadian + p_res_nmetacell + p_density

p_res_nmetacell + annotate("text", x = 7, y = 120, label = "resolution=2.5\nn_metacell=25") -> p_res_nmetacell_1
p_roc$rna +
  theme(legend.position = c(0.55, 0.3), legend.direction = "horizontal", legend.box = "horizontal", legend.background = element_rect(fill = NA)) +
 ggtitle(NULL) + labs(color = NULL) + 
  annotate("text", label = "Resolution", x = 0.55, y = 0.6, size = 3) + 
  theme(legend.key.size = unit(0.4, 'cm')) -> p_roc_rna
p_roc$gene_activity + 
  theme(legend.position = c(0.55, 0.3), legend.direction = "horizontal", legend.box = "horizontal", legend.background = element_rect(fill = NA)) +
  ggtitle(NULL) + labs(color = NULL) + 
  annotate("text", label = "Resolution", x = 0.55, y = 0.6, size = 3) + 
  theme(legend.key.size = unit(0.4, 'cm')) -> p_roc_gene_activity
p_auc$p_0.01 + theme(legend.direction = "horizontal", legend.box = "horizontal", legend.position = c(0.5,0.2), legend.background = element_rect(fill = NA)) -> p_auc_1
p_n_circadian + theme(legend.position = c(0.5,0.2), legend.direction = "horizontal", legend.box = "horizontal", legend.background = element_rect(fill = NA)) -> p_n_circadian_1
p_density + theme(legend.position = c(0.8, 0.9), legend.direction = "horizontal", legend.box = "horizontal", legend.background = element_rect(fill = NA), legend.key.size = unit(0.4, "cm")) -> p_density_1

(p_res_nmetacell_1 + (p_roc_rna / p_roc_gene_activity))/(p_auc_1+p_n_circadian_1)/p_density_1 -> p_resolution_tuning
####

# 5. Hepatocytes spatial and circadian analysis ----
setwd("~/Library/CloudStorage/Dropbox/singulomics/github")

library(tidyverse)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(patchwork)

load("~/Dropbox/singulomics/github_rda/trajectory_analysis.rda")
rm(ATres, heatdata, p_list, RNA_norm, sce);gc()

sc$pseudotime=1-sc$pseudotime

table(sc$ZT)
table(sc$group)
sc$group=droplevels(sc$group)

sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))

FeaturePlot(sc, features = 'pseudotime', reduction = "multicca.umap") + 
  annotate("text", x = -5.5, y = -7, label = "CV") + 
  annotate("text", x = 3, y = 4, label = "PN") +
  ggtitle(NULL) + 
  theme(legend.position = c(0.05, 0.8), legend.direction = "horizontal", legend.box = "horizontal", legend.text = element_text(size = 8), legend.background = element_rect(fill=NA)) + 
  annotate("text", label = "Pseudotime", x = -4.1, y = 3.75, size = 3) -> p_featureplot_pseudotime

VlnPlot(sc, features='pseudotime', group.by='ZT', pt.size = 0) + 
  theme(legend.position = "none") + 
  ggtitle(NULL) + 
  scale_fill_manual(values = rep("grey", 6)) + 
  xlab(NULL) + 
  ylab("pseudotime") -> p_violinplot_ZT #Supp_Fig8G

####VlnPlot hepatocyte subtype
res_ = 2.5
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=res_)

pseudotime.median=aggregate(sc$pseudotime, by=list(sc$seurat_clusters), FUN=median)
pseudotime.median=cbind(pseudotime.median,rank(pseudotime.median[,2])-1)
new.cluster.id=rep(NA, ncol(sc))

for(i in 1:nrow(pseudotime.median)){
  new.cluster.id[which(sc$seurat_clusters==pseudotime.median[i,1])]=pseudotime.median[i,3]
}
sc$seurat_clusters=factor(as.numeric(new.cluster.id))
rm(pseudotime.median); rm(new.cluster.id)

sc$cluster=sc$seurat_clusters

VlnPlot(sc, features='pseudotime', group.by='cluster', pt.size = 0) + 
  theme(legend.position = "none") + 
  ggtitle(NULL) + 
  xlab("Hepatocytes Subtype") + 
  ylab("pseudotime") + 
  theme(axis.text.x = element_text(angle=90, size = 8)) -> p_violinplot_cluster #Supp_Fig8F

# Plot Hepatocytes pseudotime
FeaturePlot(sc, features = "pseudotime", reduction = "multicca.umap")

# Plot cyp2f2 and cyp2e1 expression in hepatoytes umap
DefaultAssay(sc) <- 'SCT'
c("Cyp2f2", "Cyp2e1") %>% 
  "names<-"(.,.) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = y
    FeaturePlot(sc, features = x, reduction = "multicca.umap", cols = c("lightgrey", "#0072B2")) + 
      annotate("text", x = -5.5, y = -7, label = "CV") + 
      annotate("text", x = 3, y = 4, label = "PN") +
      ggtitle(x) + 
#      theme(legend.position = c(0.05, 0.8), legend.direction = "horizontal", legend.box = "horizontal", legend.text = element_text(size = 8), legend.background = element_rect(fill=NA)) + 
#      annotate("text", label = "Gene expression", x = -4.1, y = 3.75, size = 3) + 
      ggtitle(gene_) + 
      theme_classic()
  }) -> p_transient_rna #Fig_4C

DefaultAssay(sc) <- 'gene_activity'
c("Cyp2f2", "Cyp2e1") %>% 
  "names<-"(.,.) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = y
    FeaturePlot(sc, features = x, reduction = "multicca.umap", cols = c("lightgrey", "#EA3C53")) + 
      annotate("text", x = -5.5, y = -7, label = "CV") + 
      annotate("text", x = 3, y = 4, label = "PN") +
      ggtitle(x) + 
#      theme(legend.position = c(0.05, 0.8), legend.direction = "horizontal", legend.box = "horizontal", legend.text = element_text(size = 8), legend.background = element_rect(fill=NA)) + 
#      annotate("text", label = "Gene activity", x = -4.1, y = 3.75, size = 3) + 
      ggtitle(gene_) + 
      theme_classic()
  }) -> p_transient_gene_activity #Fig_4C

# Plot circadian genes, spatial genes, and spatial+circadian genes (Fig 4D)
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.rna.rds") -> rna.p.value
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.gene_activity.rds") -> atac.p.value
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_rna.rds") -> output.rna.list
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_gene_activity.rds") -> output.atac.list

c("Arntl", "Glul", "Cyp2e1", "Gsta3", "Acly", "Pck1", "Elovl3", "Arg1") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.rna.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    rna.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nCircadian p.adj=%s\nTransient p.adj=%s\nCircadian+Transient p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p

c("Arntl", "Glul", "Cyp2e1", "Gsta3", "Acly", "Pck1", "Elovl3", "Arg1") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.atac.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    atac.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nCircadian p.adj=%s\nTransient p.adj=%s\nCircadian+Transient p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
    p + scale_color_gradient(low = "grey", high = "red") -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p

# FISH validation (Acly, Pck1, Elovl3)
##Validation by smRNA Fish ----
library(R.matlab)
list_1 = readMat("~/Dropbox/singulomics/RNA_FISH/Archive/new_smFISH_val_genes.mat")
list_2 = readMat("~/Dropbox/singulomics/RNA_FISH/Archive/smFISH_validation2.mat")

list(
  Acly = list(
    ZT0 = list_2$Acly.table.ZT0[[1]][[1]],
    ZT12 = list_2$Acly.table.ZT12[[1]][[1]]
  ),
  Pck1 = list(
    ZT0 = list_2$Pck1.table.ZT0[[1]][[1]],
    ZT12 = list_2$Pck1.table.ZT12[[1]][[1]]
  ), 
  Elovl3 = list(
    ZT0 = list_2$Elovl3.table.ZT0[[1]][[1]],
    ZT12 = list_2$Elovl3.table.ZT12[[1]][[1]]
  ), 
  Arg1 = list(
    ZT0 = list_1$table.Arg1.ZT0[[1]][[1]],
    ZT12 = list_1$table.Arg1.ZT12[[1]][[1]]
  )    
) -> list_

names(list_) %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    c("ZT0", "ZT12") %>% 
      map(function(x){
        ZT_ = x
        list_[[gene_]][[ZT_]] -> expr_mat
        log(expr_mat, 2) -> expr_mat
        colMeans(expr_mat) -> expr_
        apply(expr_mat, 2, sd) -> sd_
        data.frame(mean_expr = expr_, sd = sd_) %>% 
          dplyr::mutate(dist = seq(0,1,length.out=30)) %>% 
          dplyr::mutate(group = ZT_) -> df_1
      }) %>% do.call(rbind, .) -> df_
    
    df_ %>% 
      ggplot(aes(x = dist, y = mean_expr, color = group)) + 
#      ggplot(aes(x = dist, y = mean_expr, linetype = group)) + 
      geom_smooth(se = FALSE) + 
      scale_color_manual(values = c(ZT0 = "#E08B00", ZT12 = "#CF78FF")) + 
#      geom_ribbon(aes(ymin = mean_expr - sd, ymax = mean_expr + sd), alpha = 0.2, color = NA, fill = "grey") -> p
      geom_ribbon(aes(ymin = mean_expr - sd, ymax = mean_expr + sd, fill = group), alpha = 0.2, color = NA) + 
    scale_fill_manual(values = c(ZT0 = "#E08B00", ZT12 = "#CF78FF")) -> p
    p + 
      theme_classic() + 
      ggtitle(gene_) + 
      ylab("mRNA concentration [log2]") + 
      xlab("Distance from central vein (µm)") -> p
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p

# Plot Actb and Agxt (Fig 4F)
c("Actb", "Agxt") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.rna.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    rna.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nCircadian p.adj=%s\nTransient p.adj=%s\nCircadian+Transient p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 2, guides = "collect") -> p

c("Actb", "Agxt") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.atac.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    atac.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nCircadian p.adj=%s\nTransient p.adj=%s\nCircadian+Transient p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
    p + scale_color_gradient(low = "grey", high = "red") -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 2, guides = "collect") -> p

# FISH validation agxt (Fig 4F)
list_QC = list()
c("Agxt") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(x){
    gene_ = x
    c("_ZT00_", "_ZT04_", "_ZT08_", "_ZT12_", "_ZT16_", "_ZT20_") %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      purrr::map(function(x){
        ZT_ = x
        dir_ = list.dirs("~/Dropbox/singulomics/RNA_FISH", full.names = T, recursive = F) %>% {.[grepl(gene_, .)]}
        list.files(dir_, pattern = "dots_measurements\\.csv", full.names = T, recursive = T) %>% 
          {.[grepl(ZT_,.)]} -> files_
        i_ = 0
        files_ %>% 
#          .[1] %>% 
#          .[1:3] %>% 
          purrr::map(function(x){
            ZT_ = gsub("_", "", ZT_)
            i_ + 1 ->> i_
            print(x)
            read.csv(x, header = T, stringsAsFactors = F) -> df_
            df_ %>% dplyr::filter(NucIdx != 0) -> df_
            df_ %>% dplyr::mutate(dist_ratio = Distance_PC/Distance_PP) -> df_
            df_$ZT = ZT_
            
            df_ %>% 
              dplyr::mutate(group = case_when(
                (ZT != "ZT16")&(Distance_PC < 3000)&(Distance_PP < 3000) ~ "selected", 
                (ZT == "ZT16")&(Distance_PC < 1500)&(Distance_PP < 1500) ~ "selected", 
#                Distance_PC < 2000 ~ "selected", 
                TRUE ~ "non-selected"
                )) -> df_
            
            df_ %>% 
              ggplot(aes(x = XM, y = YM, color = group)) + 
              geom_point(size = 1) -> p0
            df_ %>% dplyr::filter(group == "selected") %>% 
              ggplot(aes(x = Distance_PC)) +
              geom_histogram(bins = 200) -> p1
            patchwork::wrap_plots(p0, p1, ncol = 2) -> p
            list_QC[[ZT_]][[i_]] <<- p
            
            total_nuc = df_$Tot_nuc %>% unique()
            df_ %>% dplyr::filter(group == "selected") -> df_
            df_$Distance_PC_scaled = scales::rescale(df_$Distance_PC, to = c(0,1))
            
            seq(0,1, length.out = 30) -> interval_
            2:length(interval_) %>% 
              purrr::map(function(i){
                df_ %>% dplyr::filter(Distance_PC_scaled >= interval_[i-1] & Distance_PC_scaled < interval_[i]) -> df_
                df_ %>% dplyr::arrange(dist_ratio)
                RNA_concentration = nrow(df_)/total_nuc
                data.frame(Gene = gene_, ZT = ZT_, interval = interval_[i], RNA_concentration = RNA_concentration, file_idx = i_)
              }) %>% do.call(rbind, .)
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .) -> df_
    df_ %>% 
      dplyr::filter(!((ZT == "ZT00")&(file_idx %in% c(1,5)))) %>% 
      dplyr::filter(!((ZT == "ZT04")&(file_idx %in% c(4)))) %>% 
      dplyr::filter(!((ZT == "ZT12")&(file_idx %in% c(4)))) %>% 
      dplyr::filter(!((ZT == "ZT16")&(file_idx %in% c(1,5)))) %>% 
      dplyr::filter(!((ZT == "ZT16")&(file_idx %in% c(1,5)))) ->> df_1
    
    df_1 %>% 
      group_by(Gene, ZT, interval) %>% 
      summarise(sd_RNA_concentration = sd(RNA_concentration), mean_RNA_concentration = mean(RNA_concentration)) -> df_
  }) -> list_summary_1

list_summary_1$Agxt %>% 
  dplyr::mutate(log2_mean = log(mean_RNA_concentration, 2), cv = sd_RNA_concentration/mean_RNA_concentration) %>% 
  dplyr::mutate(log2_sd = sqrt(log(1+cv^2, 2))) %>% 
  dplyr::mutate(log2_sd = runmed(log2_sd, k = 9, endrule = "median"), log2_mean = runmed(log2_mean, k = 9, endrule = "median")) %>% 
  ggplot2::ggplot(aes(x = interval, y = log2_mean, group = ZT, color = ZT, fill = ZT)) + 
  geom_ribbon(aes(ymin = log2_mean - log2_sd, ymax = log2_mean + log2_sd), alpha = 0.1, color = NA) +
  geom_line(linewidth = 1) + 
#  scale_y_continuous(trans = "log2", labels = scales::number_format(accuracy = 0.01)) + 
#  coord_cartesian(ylim = c(0.02, 2.3)) + 
  xlab("Destance from CV") + 
  ylab("log2[RNA concentration]") -> p1_1

# Plot clock genes (Supp Fig 18)
c("Clock", "Npas2", "Bhlhe40", "Bhlhe41", "Dbp", "Nfil3", "Nr1d1", "Rorc", "Cry1", "Ciart", "Per1", "Per2") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.rna.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    rna.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nC p.adj=%s\nS p.adj=%s\nS&C p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
    p + ylab("scRNA-seq expression")
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p

c("Clock", "Npas2", "Bhlhe40", "Bhlhe41", "Dbp", "Nfil3", "Nr1d1", "Rorc", "Cry1", "Ciart", "Per1", "Per2") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    output.atac.list[[gene_]] %>% 
      dplyr::mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    atac.p.value %>% filter(Gene == gene_) -> df_
    circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
    transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
    both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    p + 
      ggtitle(sprintf("%s\nCircadian p.adj=%s\nTransient p.adj=%s\nCircadian+Transient p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p
    p + scale_color_gradient(low = "grey", high = "red") -> p
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 4, guides = "collect") -> p

# Plot venn diagram circadian genes, spatial genes and spatial + circadian genes

readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.rna.rds") -> rna.p.value
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.gene_activity.rds") -> atac.p.value
library(VennDiagram)
library(grid)
run_venn = function(){
  list(p_val = rna.p.value) %>% 
    map(function(x){
      df_ = x
      list(
        `Circadian` = df_ %>% filter(Circadian_p.adj < 0.01) %>% .$Gene,
        `Transient` = df_ %>% filter(Transient_p.adj < 0.01) %>% .$Gene,
        `Circadian+Transient` = df_ %>% filter(Circadian_transient_p.adj < 0.01) %>% .$Gene
      )
    }) %>% .[[1]] %>% 
    venn.diagram(., fill = c("#E69F00", "#56B4E9", "#CC79A7"), 
                 alpha = c(0.8, 0.8, 0.8), lwd =0, filename=NULL,
                 disable.logging = TRUE) %>% 
    grobTree()
}
p_venn_ = run_venn()
p_venn_ = ggplot() + annotation_custom(p_venn_) # RNA expression venn diagram (Supp Fig 17)

run_venn = function(){
  list(p_val = atac.p.value) %>% 
    map(function(x){
      df_ = x
      list(
        `Circadian` = df_ %>% filter(Circadian_p.adj < 0.01) %>% .$Gene,
        `Transient` = df_ %>% filter(Transient_p.adj < 0.01) %>% .$Gene,
        `Circadian+Transient` = df_ %>% filter(Circadian_transient_p.adj < 0.01) %>% .$Gene
      )
    }) %>% .[[1]] %>% 
    venn.diagram(., fill = c("#E69F00", "#56B4E9", "#CC79A7"), 
                 alpha = c(0.8, 0.8, 0.8), lwd =0, filename=NULL,
                 disable.logging = TRUE) %>% 
    grobTree()
}
p_venn_ = run_venn()
p_venn_ = ggplot() + annotation_custom(p_venn_) # ATAC activity venn diagram (Supp Fig 20)

# Plot Top 10 genes for each group
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.rna.rds") -> rna.p.value
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_rna.rds") -> output.rna.list

rna.p.value %>% dplyr::filter((Circadian_transient_p.adj < 0.01)|(Transient_p.adj < 0.01)) %>% .$Gene %>% 
  map(function(x){
    gene_ = x 
    rna.p.value %>% filter(Gene == gene_) %>% .$Circadian_transient_p.adj -> p.adj
    output.rna.list[[gene_]] -> df_
    cor.test(x = df_$pseudotime, y = df_$Mean_expression, method = "pearson") -> pearson_
    cor.test(x = df_$pseudotime, y = df_$Mean_expression, method = "spearman") -> spearman_
    data.frame(gene = gene_, 
               pearson_r = pearson_$estimate, pearson_p = pearson_$estimate, 
               spearman_r = spearman_$estimate, spearman_p = spearman_$estimate, 
               ZR_p.adj = p.adj) -> df_
  }) %>% do.call(rbind, .) %>% 
  mutate(zonation = case_when(
    pearson_r > 0 ~ "PN",
    TRUE ~ "CV"
  )) -> cor_df_all

list(p_val = rna.p.value) %>% 
  map(function(x){
    df_ = x
    list(
      `Circadian` = df_ %>% filter(Circadian_p.adj < 0.01) %>% .$Gene,
      `Transient` = df_ %>% filter(Transient_p.adj < 0.01) %>% .$Gene,
      `Circadian+Transient` = df_ %>% filter(Circadian_transient_p.adj < 0.01) %>% .$Gene
    )
  }) %>% .[[1]] %>% 
  ggvenn::list_to_data_frame() %>% {
    df_ = .
    df_ %>% dplyr::filter(if_all(!matches("key"), ~ . == TRUE)) -> CT_genes
    df_ %>% dplyr::filter(Circadian == TRUE, Transient == FALSE, `Circadian+Transient` == FALSE) -> C_genes
    df_ %>% dplyr::filter(Circadian == FALSE, Transient == TRUE, `Circadian+Transient` == FALSE) -> T_genes
    df_ %>% dplyr::filter(Circadian == TRUE, Transient == FALSE, `Circadian+Transient` == TRUE) -> `C+CT-T_genes`
    df_ %>% dplyr::filter(Circadian == FALSE, Transient == TRUE, `Circadian+Transient` == TRUE) -> `T+CT-C_genes`
    list(
      C_genes = C_genes$key,
      `C+CT-T_genes` = `C+CT-T_genes`$key,
      CT_genes = CT_genes$key,
      `T+CT-C_genes` = `T+CT-C_genes`$key, 
      T_genes = T_genes$key
    )
  } -> venn_gene_list

venn_gene_list %>% 
#  .[1:4] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    gene_ = x
    rna.p.value %>% as.data.frame() %>% 
#      column_to_rownames("Gene") %>% 
      "rownames<-"(., .$Gene) %>% 
      .[gene_, ] %>% as_tibble() -> df_
    
    if (y == "C_genes"){
      df_ %>% dplyr::arrange(Circadian_p.adj) -> df_
    }
    if (y == "C+CT-T_genes"){
      left_join(x = df_, y = cor_df_all, by = c("Gene" = "gene")) %>% dplyr::mutate(pearson_r_abs = abs(pearson_r), .after = 9) -> df_
      df_ %>% dplyr::arrange(Circadian_p.adj, Circadian_transient_p.adj, desc(pearson_r_abs)) %>% dplyr::filter(Transient_p.adj > 0.1) -> df_
    }
    if (y == "CT_genes"){
      left_join(x = df_, y = cor_df_all, by = c("Gene" = "gene")) %>% dplyr::mutate(pearson_r_abs = abs(pearson_r), .after = 9) -> df_
      df_ %>% dplyr::arrange(Circadian_transient_p.adj, desc(pearson_r_abs)) -> df_
    }
    if (y == "T+CT-C_genes"){
      left_join(x = df_, y = cor_df_all, by = c("Gene" = "gene")) %>% dplyr::mutate(pearson_r_abs = abs(pearson_r), .after = 9) -> df_
      df_ %>% dplyr::arrange(Transient_p.adj, Circadian_transient_p.adj, desc(pearson_r_abs)) -> df_
    }
    if (y == "T_genes"){
      left_join(x = df_, y = cor_df_all, by = c("Gene" = "gene")) %>% dplyr::mutate(pearson_r_abs = abs(pearson_r), .after = 9) -> df_
      df_ %>% dplyr::arrange(Transient_p.adj, desc(pearson_r_abs)) -> df_
    }
    
#    return(df_$Gene[1:10])
    return(df_)
  }) -> venn_gene_list_1

venn_gene_list_1 %>% 
  "names<-"(.,names(.)) %>% 
  map(function(x){
    df_ = x
    genes_ = df_$Gene[1:10]
    genes_ %>% 
      map(function(x){
        gene_ = x
        output.rna.list[[gene_]] -> df_
        df_ %>% mutate(ZT = gsub("ZT(.+)", "\\1", ZT) %>% as.numeric()) -> df_
        df_ %>% 
          ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
          geom_point() +
          geom_line() + 
          theme_classic() + 
          scale_x_continuous(breaks = seq(2, 22, 4)) -> p 
        
        rna.p.value %>% filter(Gene == gene_) -> df_
        circadian_padj = df_$Circadian_p.adj %>% sprintf("%.2e", .)
        transient_padj = df_$Transient_p.adj %>% sprintf("%.2e", .)
        both = df_$Circadian_transient_p.adj %>% sprintf("%.2e", .)
        
        p + 
          ggtitle(sprintf("%s\nC p.adj=%s\nT p.adj=%s\nC&T p.adj=%s", gene_, circadian_padj, transient_padj, both)) -> p 
      }) -> p_list
    patchwork::wrap_plots(p_list, ncol = 5, guides = "collect") -> p
  }) -> p_list #Supp_Fig17

# KEGG analysis of circadian + transient genes (PN vs CV)
library(clusterProfiler)
library(org.Mm.eg.db)

readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.rna.rds") -> rna.p.value
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_rna.rds") -> output.rna.list
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref

list(p_val = rna.p.value) %>% 
  map(function(x){
    df_ = x
    list(
      `Circadian` = df_ %>% filter(Circadian_p.adj < 0.01) %>% .$Gene,
      `Transient` = df_ %>% filter(Transient_p.adj < 0.01) %>% .$Gene,
      `Circadian+Transient` = df_ %>% filter(Circadian_transient_p.adj < 0.01) %>% .$Gene
    )
  }) %>% .[[1]] %>% 
  purrr::reduce(., intersect) -> intersected_genes

c("CV", "PN") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    group_ = x
    cor_df %>% filter(zonation == group_) -> df_
    Gene_ = df_$gene
    list(Gene = Gene_) %>% 
      map(function(x){
        gene_ = x
        gene.ref[gene.ref$gene_name %in% gene_,] %>% 
          .$gene_id -> ENSEMBL_id
        
        ego_CC = enrichGO(gene = ENSEMBL_id, 
                          OrgDb = org.Mm.eg.db, 
                          ont = "CC", 
                          pAdjustMethod = "BH", 
                          minGSSize = 1, 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          readable = TRUE, 
                          keyType = "ENSEMBL")
        
        ego_BP = enrichGO(gene = ENSEMBL_id, 
                          OrgDb = org.Mm.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          minGSSize = 1, 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          readable = TRUE, 
                          keyType = "ENSEMBL")
        
        ego_MF = enrichGO(gene = ENSEMBL_id, 
                          OrgDb = org.Mm.eg.db, 
                          ont = "MF", 
                          pAdjustMethod = "BH", 
                          minGSSize = 1, 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          readable = TRUE, 
                          keyType = "ENSEMBL")
        
        select(org.Mm.eg.db, keys = ENSEMBL_id, 
               columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
               keytype = "ENSEMBL") %>% drop_na() -> entrez_id
        kk <- enrichKEGG(gene = entrez_id$ENTREZID, organism ="mmu", pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05, minGSSize = 1,
                         #readable = TRUE, 
                         use_internal_data =FALSE)
        
        list(ego_CC = ego_CC, ego_BP = ego_BP, ego_MF = ego_MF, kk = kk)
      }) -> res_go_analysis
  }) -> res_go_analysis

c("CV", "PN") %>% 
  map(function(x){
    group_ = x
    res_go_analysis[[x]]$Gene$kk@result %>% 
      .[1:8, ] %>% 
#      .[1:9, ] %>% 
      mutate(Description = gsub("(.+) -.+", "\\1", Description)) %>% 
      mutate(Description = factor(Description, levels = rev(Description))) %>% 
      ggplot(aes(x = -log(p.adjust, 10), y = Description)) + 
      geom_bar(stat = "identity", width = 0.8) + 
      scale_y_discrete(labels = scales::wrap_format(25)) + 
      theme(axis.text.y = element_text(size = 12)) + 
      theme_classic() + 
      ggtitle(group_)
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 2) # Supp Fig 19B

# Plot heatmap of the 1504 spatial + circadian genes
output.rna.list[intersected_genes] %>% 
  #  .[1:2] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    df_ = x
    gene_ = y
    0:24 %>% 
      map(function(x){
        metacell_ = sprintf("metacell_%s", x)
        df_ %>% dplyr::filter(cluster == x) -> df_
        expression = mean(df_$Mean_expression)
        pseudotime = mean(df_$pseudotime)
        data.frame(gene=gene_, expression=expression, metacell = metacell_)
      }) %>% do.call(rbind, .) -> df_
    df_ %>% pivot_wider(names_from = metacell, values_from = expression)
  }) %>% do.call(rbind, .) %>% 
  as_tibble() %>%
  column_to_rownames("gene") %>% 
  apply(., 1, function(x){scales::rescale(x, to=c(0,1))}) %>% 
  t() -> RNA_matrix

readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_gene_activity.rds") -> output.atac.list
output.atac.list[intersected_genes] %>% 
  #  .[1:2] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    df_ = x
    gene_ = y
    0:24 %>% 
      map(function(x){
        metacell_ = sprintf("metacell_%s", x)
        df_ %>% dplyr::filter(cluster == x) -> df_
        expression = mean(df_$Mean_expression)
        pseudotime = mean(df_$pseudotime)
        data.frame(gene=gene_, expression=expression, metacell = metacell_)
      }) %>% do.call(rbind, .) -> df_
    df_ %>% pivot_wider(names_from = metacell, values_from = expression)
  }) %>% do.call(rbind, .) %>% 
  as_tibble() %>%
  column_to_rownames("gene") %>% 
  apply(., 1, function(x){scales::rescale(x, to=c(0,1))}) %>% 
  t() -> ATAC_matrix

cor_df %>% mutate(cor = abs(pearson_r)) %>% 
  dplyr::arrange(ZR_p.adj, desc(cor)) %>% filter(abs(cor) > 0.7) %>% 
  head(20) %>% 
  .$gene %>% 
  list(gene = .) %>% 
  map(function(x){
    gene_ = x
    rownames(RNA_matrix) %>% ifelse(. %in% gene_, ., "") -> labels_
    pheatmap::pheatmap(RNA_matrix, 
                       cluster_rows = T, 
                       cutree_rows = 4,
                       cluster_cols = F, 
                       labels_row = labels_,
                       #                       show_rownames = F, 
                       show_colnames = F, 
                       fontsize = 5, 
                       main = "Circadian+transient: 1504 genes") -> p
  }) %>% .[[1]] -> p_heatmap

row_clusters <- cutree(p_heatmap$tree_row, k = 4)

c(4,1,2,3) %>% 
  map(function(x){
    row_clusters[row_clusters==x] %>% names() -> gene_
  }) %>% purrr::reduce(., c) %>% 
  {
    gene_ = .
    labels_ = c(
      #Single-cell spatial reconstruction reveals global division of labour in the mammalian liver
      c("Glul", "Cyp2e1", "Ass1", "Asl", "Alb", "Cyp2f2"),
      #Spatial Transcriptomics to define transcriptional patterns of zonation and structural components in the mouse liver
      c("Sds", "Cyp2f2", "Hal", "Hsd17b13", "Aldh1b1", "Glul", "Oat", "Slc1a2", "Cyp2e1", "Cyp2a5"), 
      #Space-time logic of liver gene expression at sublobular scale
      c("Glul", "Ass1", "Elovl3", "Pck1", "Acly", "Arg1")
    )
    print(labels_)
    RNA_matrix[gene_,-c(8,19,21)] -> RNA_matrix
    rownames(RNA_matrix) %>% ifelse(. %in% labels_, ., "") -> labels_
    ggplotify::as.ggplot(
      pheatmap::pheatmap(RNA_matrix, 
                         cluster_rows = F, 
                         #                       cutree_rows = 4,
                         cluster_cols = F, 
                         labels_row = labels_,
                         #                       show_rownames = F, 
                         show_colnames = F, 
                         fontsize = 5, 
                         main = "Circadian+transient: 1504 genes")       
    ) -> p
  } -> p_RNA

c(4,1,2,3) %>% 
  map(function(x){
    row_clusters[row_clusters==x] %>% names() -> gene_
  }) %>% purrr::reduce(., c) %>% 
  {
    gene_ = .
    labels_ = c(
      #Single-cell spatial reconstruction reveals global division of labour in the mammalian liver
      c("Glul", "Cyp2e1", "Ass1", "Asl", "Alb", "Cyp2f2"),
      #Spatial Transcriptomics to define transcriptional patterns of zonation and structural components in the mouse liver
      c("Sds", "Cyp2f2", "Hal", "Hsd17b13", "Aldh1b1", "Glul", "Oat", "Slc1a2", "Cyp2e1", "Cyp2a5"), 
      #Space-time logic of liver gene expression at sublobular scale
      c("Glul", "Ass1", "Elovl3", "Pck1", "Acly", "Arg1")
    )
    print(labels_)
    ATAC_matrix[gene_,-c(8,19,21)] -> ATAC_matrix
    rownames(ATAC_matrix) %>% ifelse(. %in% labels_, ., "") -> labels_
    
    ggplotify::as.ggplot(
      pheatmap::pheatmap(ATAC_matrix, 
                         cluster_rows = F, 
                         #                       cutree_rows = 4,
                         cluster_cols = F, 
                         labels_row = labels_,
                         #                       show_rownames = F, 
                         show_colnames = F, 
                         fontsize = 5, 
                         main = "Circadian+transient: 1504 genes")       
    ) -> p

  } -> p_ATAC

patchwork::wrap_plots(p_RNA, p_ATAC, ncol = 2, guides = "collect") #Supp Fig 19A

# 6. Hepatocytes spatial and circadian CREs ----

# Read data
library(tidyverse)
library(Seurat)
library(Signac)

load("~/Dropbox/singulomics/github_rda/trajectory_analysis.rda")
rm(ATres, heatdata, p_list, RNA_norm, sce);gc()

sc$pseudotime=1-sc$pseudotime

sc$group=droplevels(sc$group)

sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))

res_ = 2.5
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=res_)

pseudotime.median=aggregate(sc$pseudotime, by=list(sc$seurat_clusters), FUN=median)
pseudotime.median=cbind(pseudotime.median,rank(pseudotime.median[,2])-1)
new.cluster.id=rep(NA, ncol(sc))

for(i in 1:nrow(pseudotime.median)){
  new.cluster.id[which(sc$seurat_clusters==pseudotime.median[i,1])]=pseudotime.median[i,3]
}
sc$seurat_clusters=factor(as.numeric(new.cluster.id))
rm(pseudotime.median); rm(new.cluster.id)

sc$cluster=sc$seurat_clusters

# Link ATAC peaks to gene expression

LinkPeaks_custom = function(object, peak.assay, expression.assay, peak.slot = "counts", 
                            expression.slot = "counts", min.cells = 10, 
                            genes.use = NULL, distance = 5e+05){
  
  annot <- Annotation(object = object[[peak.assay]])
  gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annot)
  
  #  meta.features <- GetAssayData(object = object, assay = peak.assay, 
  #                                slot = "meta.features")
  peak.data <- GetAssayData(object = object, assay = peak.assay, 
                            slot = peak.slot)
  expression.data <- GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  peakcounts <- rowSums(x = peak.data > 0)
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  
  if (!is.null(x = genes.use)) {
    genes.keep <- intersect(x = names(x = genes.keep[genes.keep]), 
                            y = genes.use)
  }
  expression.data <- expression.data[genes.keep, , drop = FALSE]
  genes <- rownames(x = expression.data)
  
  gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes, ]
  
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- Signac:::DistanceToTSS(peaks = peaks, genes = gene.coords.use, 
                                                 distance = distance)
  
  return(peak_distance_matrix)
}

LinkPeaks_custom(object = sc, peak.assay = "ATAC", expression.assay = "RNA") -> peak_distance_matrix

colSums(peak_distance_matrix) %>% 
  {.[.!=0]} %>% 
  names() -> gene_select

peak_distance_matrix <- peak_distance_matrix[,gene_select]

# Generate circadian_gene_peak_expression_matrix and spatial_gene_peak_expression_matrix

##Link by zone##
args[1] %>% as.numeric() -> se
seq(1, ncol(peak_distance_matrix), length.out = 50) %>% 
round() %>%  
{
    seq_ = .
    1:49 %>% 
    map(function(i){
        start = seq_[i]
        end = seq_[i+1]-1
        if (i == 49){end = seq_[i+1]}
        data.frame(start = start, end = end)
    }) %>% do.call(rbind, .)
} -> start_end_df

start_ = start_end_df[se, "start"]
end_ = start_end_df[se, "end"]
output_dir = "./rda/linked_peak/"

sc@meta.data %>% 
  dplyr::mutate(idx = as.integer(cluster)) %>% 
  #  .$idx %>% table()
  dplyr::mutate(zonations = case_when(
    idx %in% 1:4 ~ "zone1", 
    idx %in% 5:9 ~ "zone2",
    idx %in% 10:13 ~ "zone3",
    idx %in% 14:17 ~ "zone4",
    idx %in% 18:21 ~ "zone5",
    idx %in% 22:25 ~ "zone6",
  )) %>% dplyr::select(-idx) -> sc@meta.data
factor(sc$zonations) -> sc$zonations

set.seed(123)
sc@meta.data %>% 
rownames_to_column("cellnames") %>% 
group_by(zonations) %>% 
group_map(function(x,y){
  df_ = x
  set.seed(123)
  sprintf("Rep_%s", 1:4) %>% sample(., nrow(df_), replace = T) -> Rep_
#  table(Rep_)
  df_$Reps = Rep_
  df_
}, .keep = T) %>% do.call(rbind, .) %>% as.data.frame() %>% 
column_to_rownames("cellnames") %>% 
.[rownames(sc@meta.data),] -> sc@meta.data
factor(sc$Reps) -> sc$Reps

library(furrr)
library(future)
20000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 2)

i = 0
colnames(peak_distance_matrix) %>% 
  "names<-"(.,.) %>% 
  .[start_:end_] %>% 
#  .[1:100] %>% 
#  map(function(gene_){
  future_map(function(gene_){
    i <<- i+1
    if(i%%10==0){print(i)}
    peak_distance_matrix[,gene_] -> peaks_
    peaks_ %>% {.[.==1]} %>% names() -> peaks
    levels(sc$zonations) %>% 
      map(function(zonations_){
        sc@meta.data %>% dplyr::filter(zonations == zonations_) -> df_meta
        levels(sc$Reps) %>% 
        map(function(Reps_){
          df_meta %>% dplyr::filter(Reps == Reps_) %>% rownames() -> cells
          sc$nCount_RNA[cells] -> nCount_RNA
          median(sc$nCount_RNA) -> median_nCount_RNA
          sc@assays$RNA@counts[gene_, cells] %>% {(./nCount_RNA)*median_nCount_RNA} %>% mean()
        }) %>% purrr::reduce(., c) -> RNA_expr
      }) %>% purrr::reduce(., c) -> RNA_expr

    levels(sc$zonations) %>% 
      map(function(zonations_){
        sc@meta.data %>% dplyr::filter(zonations == zonations_) -> df_meta
        levels(sc$Reps) %>% 
        map(function(Reps_){
          df_meta %>% dplyr::filter(Reps == Reps_) %>% rownames() -> cells
          sc$nCount_ATAC[cells] -> nCount_ATAC
          median(sc$nCount_ATAC) -> median_nCount_ATAC
        if (length(peaks) > 1){
            sc@assays$ATAC@counts[peaks, cells] %>% {(./nCount_ATAC)*median_nCount_ATAC} %>% rowMeans() %>% 
            as.data.frame() %>% t() %>% "rownames<-"(., sprintf("%s_%s", zonations_, Reps_)) -> ATAC_expr
        }else{
            sc@assays$ATAC@counts[peaks, cells] %>% {(./nCount_ATAC)*median_nCount_ATAC} %>% mean() %>% as.matrix() %>% 
            "colnames<-"(., peaks) %>% "rownames<-"(., sprintf("%s_%s", zonations_, Reps_)) -> ATAC_expr
        }
        return(ATAC_expr)
        }) %>% do.call(rbind, .) -> ATAC_expr
      }) %>% do.call(rbind, .) -> ATAC_expr
    
    cbind(Gene_expr = RNA_expr, ATAC_expr) -> df_
    rownames(df_) = rownames(ATAC_expr)
    return(df_)
  }) -> gene_peak_matrix_spatial
saveRDS(gene_peak_matrix_spatial, file = sprintf("%sgene_peak_matrix_spatial_%s_%s_4_reps.rds", output_dir, start_, end_))

##Linked by ZT time##
args[1] %>% as.numeric() -> se
seq(1, ncol(peak_distance_matrix), length.out = 50) %>% 
round() %>%  
{
    seq_ = .
    1:49 %>% 
    map(function(i){
        start = seq_[i]
        end = seq_[i+1]-1
        if (i == 49){end = seq_[i+1]}
        data.frame(start = start, end = end)
    }) %>% do.call(rbind, .)
} -> start_end_df

start_ = start_end_df[se, "start"]
end_ = start_end_df[se, "end"]
output_dir = "./rda/linked_peak/"

set.seed(123)
sc@meta.data %>% 
rownames_to_column("cellnames") %>% 
group_by(ZT) %>% 
group_map(function(x,y){
  df_ = x
  set.seed(123)
  sprintf("Rep_%s", 1:4) %>% sample(., nrow(df_), replace = T) -> Rep_
  df_$Reps = Rep_
  df_
}, .keep = T) %>% do.call(rbind, .) %>% as.data.frame() %>% 
column_to_rownames("cellnames") %>% 
.[rownames(sc@meta.data),] -> sc@meta.data
factor(sc$Reps) -> sc$Reps

library(furrr)
library(future)
20000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 2)

i = 0
colnames(peak_distance_matrix) %>% 
  "names<-"(.,.) %>% 
  .[start_:end_] %>% 
#  .[1:10] %>% 
  future_map(function(gene_){
    i <<- i+1
    if(i%%10==0){print(i)}
    peak_distance_matrix[,gene_] -> peaks_
    peaks_ %>% {.[.==1]} %>% names() -> peaks
    levels(sc$ZT) %>% 
      map(function(ZT_){
        sc@meta.data %>% dplyr::filter(ZT == ZT_) -> df_meta
        levels(sc$Reps) %>% 
        map(function(Reps_){
          df_meta %>% dplyr::filter(Reps == Reps_) %>% rownames() -> cells
          sc$nCount_RNA[cells] -> nCount_RNA
          median(sc$nCount_RNA) -> median_nCount_RNA
          sc@assays$RNA@counts[gene_, cells] %>% {(./nCount_RNA)*median_nCount_RNA} %>% mean()
        }) %>% purrr::reduce(., c) -> RNA_expr
      }) %>% purrr::reduce(., c) -> RNA_expr

    levels(sc$ZT) %>% 
      map(function(ZT_){
        sc@meta.data %>% dplyr::filter(ZT == ZT_) -> df_meta
        levels(sc$Reps) %>% 
        map(function(Reps_){
          df_meta %>% dplyr::filter(Reps == Reps_) %>% rownames() -> cells
          sc$nCount_ATAC[cells] -> nCount_ATAC
          median(sc$nCount_ATAC) -> median_nCount_ATAC
        if (length(peaks) > 1){
            sc@assays$ATAC@counts[peaks, cells] %>% {(./nCount_ATAC)*median_nCount_ATAC} %>% rowMeans() %>% 
            as.data.frame() %>% t() %>% "rownames<-"(., sprintf("%s_%s", ZT_, Reps_)) -> ATAC_expr
        }else{
            sc@assays$ATAC@counts[peaks, cells] %>% {(./nCount_ATAC)*median_nCount_ATAC} %>% mean() %>% as.matrix() %>% 
            "colnames<-"(., peaks) %>% "rownames<-"(., sprintf("%s_%s", ZT_, Reps_)) -> ATAC_expr
        }
        return(ATAC_expr)
        }) %>% do.call(rbind, .) -> ATAC_expr
      }) %>% do.call(rbind, .) -> ATAC_expr
    
    cbind(Gene_expr = RNA_expr, ATAC_expr) -> df_
    rownames(df_) = rownames(ATAC_expr)
    return(df_)
  }) -> gene_peak_matrix_temporal

saveRDS(gene_peak_matrix_temporal, file = sprintf("%sgene_peak_matrix_temporal_%s_%s_4_reps.rds", output_dir, start_, end_))

# load the processed data
c("spatial", "temporal") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    list.files("~/Dropbox/singulomics/github_rda/linkpeak", pattern = "4_reps\\.rds", full.names = T) -> files_
    files_ %>% {.[grepl(x, .)]} %>% gtools::mixedsort() -> files_
    files_ %>% 
      #      .[1:2] %>% 
      map(function(file_){
        readRDS(file_) -> list_
      }) %>% purrr::reduce(., c) -> list_
  }) -> process_dat_list

i = 0
names(process_dat_list$temporal) %>% 
  "names<-"(.,.) %>%
#  .[1:2] %>% 
  map(function(gene_){
    i <<- i + 1
    if(i%%100==0){print(i)}
    process_dat_list$spatial[[gene_]] -> df_spatial
    process_dat_list$temporal[[gene_]] -> df_temporal
    
    colnames(df_spatial) %>% {.[-1]} -> peaks
    peaks %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(peak_){
        cor.test(df_spatial[,"Gene_expr"], df_spatial[,peak_], method = "pearson") -> cor_test
        cor_test$p.value -> spatial_pval
        cor_test$estimate -> spatial_r
        cor.test(df_temporal[,"Gene_expr"], df_temporal[,peak_], method = "pearson") -> cor_test
        cor_test$p.value -> temporal_pval
        cor_test$estimate -> temporal_r
        data.frame(gene = gene_, peak = peak_, spatial_pval = spatial_pval, spatial_r = spatial_r, 
                   temporal_pval = temporal_pval, temporal_r = temporal_r) -> df_
      }) %>% do.call(rbind, .) -> df_
    df_ %>% dplyr::mutate(spatial_p.adj = p.adjust(spatial_pval, method = "BH"), 
                    temporal_p.adj = p.adjust(temporal_pval, method = "BH")) -> df_
  }) %>% do.call(rbind, .) %>% "rownames<-"(., NULL) -> gene_peak_cor_df

readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/results/res_2.5.rna.rds") -> rna.p.value
rna.p.value %>%
  {
    df_ = .
    Circadian_Transient_Genes = df_ %>% dplyr::filter(Circadian_transient_p.adj < 0.01) %>% .$Gene
    Circadian = df_ %>% dplyr::filter(Circadian_p.adj < 0.01, Transient_p.adj > 0.01, Circadian_transient_p.adj > 0.01) %>% .$Gene
    Transient = df_ %>% dplyr::filter(Circadian_p.adj > 0.01, Transient_p.adj < 0.01, Circadian_transient_p.adj > 0.01) %>% .$Gene
    sprintf("Circadian_Transient_Genes: %d\nCircadian: %d\nTransient: %d", length(Circadian_Transient_Genes), length(Circadian), length(Transient))
    
    rbind(
      data.frame(Gene = Circadian_Transient_Genes, gene_type = "Circadian_Transient"),
      data.frame(Gene = Circadian, gene_type = "Circadian"),
      data.frame(Gene = Transient, gene_type = "Transient")
    )
  } -> gene_type_df

left_join(x = gene_peak_cor_df, y = gene_type_df, by = c("gene" = "Gene")) %>% 
  dplyr::mutate(gene_type = ifelse(is.na(gene_type), "Others", gene_type)) %>% 
  dplyr::filter((spatial_r > 0)|(temporal_r > 0)) -> gene_peak_cor_df_2

gene_peak_cor_df_2 %>% dplyr::mutate(linked_by = case_when(
  (gene_type == "Circadian")&(temporal_r > 0.5)&(temporal_p.adj < 0.05) ~ "temporal",
  (gene_type == "Transient")&(spatial_r > 0.5)&(spatial_p.adj < 0.05) ~ "spatial",
  (gene_type == "Circadian_Transient")&(spatial_r > 0.5)&(temporal_r > 0.5)&(spatial_p.adj < 0.05)&(temporal_p.adj < 0.05) ~ "spatial_temporal",
  (gene_type == "Circadian_Transient")&(spatial_r > 0.5)&(temporal_r < 0.5)&(spatial_p.adj < 0.05) ~ "spatial",
  (gene_type == "Circadian_Transient")&(spatial_r < 0.5)&(temporal_r > 0.5)&(temporal_p.adj < 0.05) ~ "temporal",
  TRUE ~ NA
)) -> gene_peak_cor_df_2

readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
GenomicRanges::promoters(gene.ref, upstream = 0, downstream = 1) -> gene.ref.promoter
gene_peak_cor_df_2$gene %>% unique() %>% 
  "names<-"(.,.) %>% 
  #  .[c("Arntl", "Glul")] %>% 
  map(function(gene_){
    gene_peak_cor_df_2 %>% dplyr::filter(gene == gene_) -> df_
    data.frame(
      seqnames = gsub("(.+)-(.+)-(.+)", "\\1", df_$peak),
      start = as.integer(gsub("(.+)-(.+)-(.+)", "\\2", df_$peak)),
      end = as.integer(gsub("(.+)-(.+)-(.+)", "\\3", df_$peak))
    ) %>% GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) -> peak_gr
    gene_gr = gene.ref.promoter[gene.ref.promoter$gene_name == gene_,]
    GenomicRanges::distance(gene_gr, peak_gr) -> dist_
    gene_tss = gene_gr@ranges@start
    gene_strand = gene_gr %>% as.data.frame() %>% .$strand
    peak_gr %>% as.data.frame() %>% dplyr::mutate(Distance_to_TSS = dist_) -> df_1
    if (gene_strand == "+"){
      df_1 %>% dplyr::mutate(Distance_to_TSS = case_when(
        end < gene_tss ~ -1*Distance_to_TSS,
        TRUE ~ Distance_to_TSS
      )) -> df_1
    }else{
      df_1 %>% dplyr::mutate(Distance_to_TSS = case_when(
        start > gene_tss ~ -1*Distance_to_TSS,
        TRUE ~ Distance_to_TSS
      )) -> df_1
    }
    df_1 %>% dplyr::mutate(peak = sprintf("%s-%s-%s", seqnames, start, end)) %>% dplyr::select(peak, Distance_to_TSS)-> df_1
    left_join(df_, df_1, by = "peak") %>% dplyr::mutate(gene_len = gene.ref[gene.ref$gene_name == gene_]@ranges@width) -> df_
  }) %>% do.call(rbind, .) -> gene_peak_cor_df_2
gene_peak_cor_df_2$linked_by = factor(gene_peak_cor_df_2$linked_by)

gene_peak_cor_df_2 %>% drop_na() %>% 
  dplyr::filter(gene_type != "Others") %>% 
  dplyr::mutate(gene_type = case_when(
    gene_type == "Circadian_Transient" ~ "Spatial & Circadian",
    gene_type == "Transient" ~ "Spatial",
    TRUE ~ gene_type
)) %>% 
  dplyr::mutate(linked_by = case_when(
    linked_by == "spatial" ~ "Spatial",
    linked_by == "temporal" ~ "Circadian", 
    linked_by == "spatial_temporal" ~ "Spatial & Circadian", 
  )) -> gene_peak_cor_df_3
gene_peak_cor_df_3 %>% drop_na() %>% {write.csv(., "~/Downloads/supplementary_table_2.csv", row.names = F, quote = F)}

# Plot Arntl, Glul and Pck1 (Fig 4G)

## Coverage plots ##
LinkPlot_custom = function(object, region, links, extend.upstream, extend.downstream, assay = "ATAC"){
  region <- Signac:::FindRegion(object = object, region = region, sep = c("-","-"), 
                                assay = assay, extend.upstream = extend.upstream, extend.downstream = extend.downstream)
  chromosome <- seqnames(x = region)
  links.keep <- IRanges::subsetByOverlaps(x = links, ranges = region)
  
  link.df <- as.data.frame(x = links.keep)
  link.df <- link.df[link.df$start >= IRanges::start(x = region) & link.df$end <= 
                       IRanges::end(x = region), ]
  link.df$link_by = factor(link.df$link_by, levels = c("temporal", "spatial", "spatial_temporal"))

  if (nrow(x = link.df) > 0){
    link.df$group <- seq_len(length.out = nrow(x = link.df))
    df <- data.frame(x = c(link.df$start, (link.df$start + 
                                             link.df$end)/2, link.df$end), y = c(rep(x = 0, 
                                                                                     nrow(x = link.df)), rep(x = -1, nrow(x = link.df)), 
                                                                                 rep(x = 0, nrow(x = link.df))), group = rep(x = link.df$group, 
                                                                                                                             3), link_by = rep(link.df$link_by, 3))    
  }
  
  p <- ggplot(data = df) + ggforce::geom_bezier(mapping = aes_string(x = "x", 
                                                                     y = "y", group = "group", color = "link_by")) + 
    scale_color_manual(values = c("temporal" = "red", "spatial" = "blue", "spatial_temporal" = "purple"), drop = F) + 
    geom_hline(yintercept = 0, color = "grey") + 
    theme_classic() + theme(axis.ticks.y = element_blank(), 
                            axis.text.y = element_blank(), legend.position = "bottom") + ylab("Links") + xlab(label = paste0(chromosome, 
                                                                                                                             " position (bp)")) + xlim(c(IRanges::start(x = region), IRanges::end(x = region)))
  
  return(p)
}

CoveragePlot_custom = function(genes, ext.upstream_ = 35000, ext.downstream_ = 35000){
  genes %>% 
    #  c("Pck1", "Acly", "Gsta3", "Arg1", "Cyp2e1") %>% 
    "names<-"(.,.) %>% 
    map(function(gene_){
      gene_peak_cor_df_2 %>% dplyr::filter(gene == gene_) -> df_
      
      df_ %>% drop_na() -> df_
      df_ %>% dplyr::filter(abs(Distance_to_TSS) < 31000+gene.ref[gene.ref$gene_name==gene_, ]@ranges@width) -> df_
      
      if ((gene.ref[gene.ref$gene_name == gene_, ] %>% as.data.frame() %>% .$strand) == "+"){
        gene_tss = gene.ref[gene.ref$gene_name == gene_, ] %>% as.data.frame() %>% .$start
      }else{
        gene_tss = gene.ref[gene.ref$gene_name == gene_, ] %>% as.data.frame() %>% .$end
      }
      
      df_ %>% dplyr::mutate(
        seqnames = gsub("(.+)-(.+)-(.+)", "\\1", peak),
        start = as.integer(gsub("(.+)-(.+)-(.+)", "\\2", peak)),
        end = as.integer(gsub("(.+)-(.+)-(.+)", "\\3", peak))
      ) %>% 
        dplyr::mutate(middle = (start + end) / 2) %>% 
        dplyr::mutate(start_1 = case_when(
          middle < gene_tss ~ middle,
          TRUE ~ gene_tss
        )) %>% 
        dplyr::mutate(end_1 = case_when(
          middle < gene_tss ~ gene_tss,
          TRUE ~ middle
        )) %>% 
        {
          df_ = .
          data.frame(seqnames = df_$seqnames, start = df_$start_1, end = df_$end_1, link_by = df_$linked_by)
        } %>% GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) -> link_
      
#      ext.upstream = 35000
#      ext.downstream = 35000
      ext.upstream = ext.upstream_
      ext.downstream = ext.downstream_
      
      
      CoveragePlot(
        object = sc,
        region = gene_,
        features = gene_,
        expression.assay = "SCT",
        extend.upstream = ext.upstream,
        extend.downstream = ext.downstream, 
        group.by = "ZT", 
        links = FALSE, 
        peaks = FALSE,
        annotation = TRUE
      ) -> p_ZT  
      patchwork::wrap_plots(p_ZT[[1]][[1]], p_ZT[[1]][[2]], ncol = 2, widths = c(10,1)) -> p_ZT_1
      
      CoveragePlot(
        object = sc,
        region = gene_,
        features = gene_,
        expression.assay = "SCT",
        extend.upstream = ext.upstream,
        extend.downstream = ext.downstream, 
        group.by = "zonations", 
        links = FALSE
      ) -> p_zonation
      patchwork::wrap_plots(p_zonation[[1]][[1]], p_zonation[[1]][[2]], ncol = 2, widths = c(10,1)) -> p_zonation_1
      
      gene_plot <- AnnotationPlot(
        object = sc,
        region = gene_,
        extend.upstream = ext.upstream,
        extend.downstream = ext.downstream 
      )
      ExpressionPlot(object = sc, features = gene_, assay = "SCT", group.by = "celltype") -> expr_plot
      patchwork::wrap_plots(gene_plot, expr_plot, ncol = 2, widths = c(10,1)) -> gene_plot_1
      
      gene.ref[gene.ref$gene_name==gene_, ] %>% as.data.frame() %>% {sprintf("%s-%s-%s", .$seqnames, .$start, .$end)} -> gene_region
      peak_plot <- PeakPlot(
        object = sc,
        region = gene_region,
        extend.upstream = ext.upstream,
        extend.downstream = ext.downstream 
      )
      
      patchwork::wrap_plots(peak_plot, expr_plot, ncol = 2, widths = c(10,1)) -> peak_plot_1
      
      LinkPlot_custom(
        object = sc,
        region = gene_,
        links = link_,
        extend.upstream = ext.upstream,
        extend.downstream =ext.downstream 
      ) -> p_link
      patchwork::wrap_plots(p_link, expr_plot, ncol = 2, widths = c(10,1)) -> p_link_1
      
      patchwork::wrap_plots(p_ZT_1, p_zonation_1, gene_plot_1, peak_plot_1, p_link_1, ncol = 1, heights = c(10,10,2,1,2)) -> p
      
    }) -> p_list  
  return(p_list)
}

sc@meta.data %>% 
  dplyr::mutate(idx = as.integer(cluster)) %>% 
  #  .$idx %>% table()
  dplyr::mutate(zonations = case_when(
    idx %in% 1:4 ~ "zone1", 
    idx %in% 5:9 ~ "zone2",
    idx %in% 10:13 ~ "zone3",
    idx %in% 14:17 ~ "zone4",
    idx %in% 18:21 ~ "zone5",
    idx %in% 22:25 ~ "zone6",
  )) %>% dplyr::select(-idx) -> sc@meta.data

CoveragePlot_custom("Arntl")
CoveragePlot_custom("Glul")
CoveragePlot_custom("Pck1")

# For Pck1
CoveragePlot(
  object = sc,
  region = 'chr2-173138113-173138997',
  expression.assay = "SCT",
  extend.upstream = 3000,
  extend.downstream = 3000, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)
CoveragePlot(
  object = sc,
  region = 'chr2-173530926-173531838',
  expression.assay = "SCT",
  extend.upstream = 3000,
  extend.downstream = 3000, 
  group.by = "zonations", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

## Motif enrichment analysis ##
DefaultAssay(sc) <- "ATAC"
c("Arntl", "Glul", "Pck1") %>% 
  "names<-"(.,.) %>% 
  map(function(gene_){
    gene_peak_cor_df_2 %>% dplyr::filter(gene == gene_) %>% 
      {
        df_ = .
        df_ %>% drop_na() -> df_
        df_$peak -> peaks_
        Signac::FindMotifs(
          sc, 
          features = peaks_,
          background = sc@assays$ATAC@counts %>% rownames()
        ) 
      }
  }) -> enriched_motifs

enriched_motifs %>% 
  map2(.x=.,.y=names(.),.f=function(df_, gene_){
    df_ %>% dplyr::filter(pvalue < 0.05) -> df_
    if (gene_ == "Arntl"){
      motif = df_[3, "motif"]
      motif.name = df_[3, "motif.name"]
      pval = df_[3, "pvalue"]
    }else if (gene_ == "Glul"){
      motif = df_[18, "motif"]
      motif.name = df_[18, "motif.name"]
      pval = df_[18, "pvalue"]
    }else{
      motif = df_[6, "motif"]
      motif.name = df_[6, "motif.name"]
      pval = df_[6, "pvalue"]
    }
    Signac::MotifPlot(
      sc,
      motif = motif,
    ) + 
      ggtitle(sprintf("p=%.2e", pval))
  }) %>% patchwork::wrap_plots(., ncol = 3) -> p #Fig 4G
####

# Coverage plots (Clock, Cry1, Dbp, Nr1d1, Rorc) (Supp Fig 21)

# Clock
gene_peak_cor_df_2 %>% dplyr::filter(gene == "Clock") %>% drop_na()
CoveragePlot_custom(genes = "Clock", ext.upstream_ = -60000, ext.downstream_ = 0)
CoveragePlot(
  object = sc,
  region = 'chr5-76286997-76287839',
  expression.assay = "SCT",
  extend.upstream = 1000,
  extend.downstream = 1000, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

#Cry1
gene_peak_cor_df_2 %>% dplyr::filter(gene == "Cry1") %>% drop_na() %>% dplyr::arrange(Distance_to_TSS)
CoveragePlot_custom("Cry1", ext.upstream_ = -40000, ext.downstream_ = 100)

CoveragePlot(
  object = sc,
  region = 'chr10-85187236-85188127',
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 500, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

CoveragePlot(
  object = sc,
  region = 'chr10-85176335-85177231',
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 500, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

CoveragePlot(
  object = sc,
  region = 'chr10-85161454-85162299',
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 500, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

CoveragePlot(
  object = sc,
  region = 'chr10-85157420-85158122',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

CoveragePlot(
  object = sc,
  region = 'chr10-85153845-85154752',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

#Dbp
gene_peak_cor_df_2 %>% dplyr::filter(gene == "Dbp") %>% drop_na() %>% dplyr::arrange(Distance_to_TSS)
CoveragePlot_custom("Dbp", ext.upstream_ = 500, ext.downstream_ = -1000)

CoveragePlot(
  object = sc,
  region = 'chr7-45744143-45745033',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

#Nr1d1
gene_peak_cor_df_2 %>% dplyr::filter(gene == "Nr1d1") %>% drop_na() %>% dplyr::arrange(Distance_to_TSS)
CoveragePlot_custom("Nr1d1", ext.upstream_ = 500, ext.downstream_ = 500)
CoveragePlot(
  object = sc,
  region = 'chr11-98783024-98783909',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)
CoveragePlot(
  object = sc,
  region = 'chr11-98766828-98767690',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)
CoveragePlot(
  object = sc,
  region = 'chr11-98764490-98765212',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

#Rorc
gene_peak_cor_df_2 %>% dplyr::filter(gene == "Rorc") %>% drop_na() %>% dplyr::arrange(Distance_to_TSS)
CoveragePlot_custom("Rorc", ext.upstream_ = 500, ext.downstream_ = -20000)

CoveragePlot(
  object = sc,
  region = 'chr3-94398032-94398603',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)

CoveragePlot(
  object = sc,
  region = 'chr3-94398803-94399629',
  expression.assay = "SCT",
  extend.upstream = 200,
  extend.downstream = 200, 
  group.by = "ZT", 
  links = FALSE, 
  peaks = TRUE,
  annotation = TRUE
)
####

# Plot number of CREs (spatial vs spatial+circadian vs circadian)
library(ggbreak)
gene_peak_cor_df_2$linked_by %>% 
  table() %>% as.data.frame() %>% 
  "colnames<-"(., c("linked_by", "count")) %>% 
  ggplot(aes(x = linked_by, y = count, fill = linked_by)) + 
  geom_bar(stat = "identity") + 
  scale_y_break(breaks = c(1000, 9000), scales = 0.4) +  
  theme_classic() #Fig_4H ----

# Plot histone mark signals of different groups of CREs
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
full_genome <- GRanges(seqnames = seqnames(genome),
                       ranges = IRanges(start = 1, end = seqlengths(genome)))

sc@assays$ATAC@counts %>% rownames() %>% 
  {
    peaks_ = .
    data.frame(seqnames = gsub("(.+)-(.+)-(.+)", "\\1", peaks_), 
               start = gsub("(.+)-(.+)-(.+)", "\\2", peaks_) %>% as.numeric(), 
               end = gsub("(.+)-(.+)-(.+)", "\\3", peaks_) %>% as.numeric()) %>% 
      GenomicRanges::makeGRangesFromDataFrame()
  } -> ATAC_peaks_gr

regions_without_peaks <- setdiff(full_genome, ATAC_peaks_gr)

library(regioneR)
# Set parameters for random region generation
num_regions <- 10000        # Number of random regions
mean_length <-  876      # Mean length of each region
sd_length <- 0          # Standard deviation of region lengths

# Generate random regions
random_regions_non_atac_peaks <- createRandomRegions(
  nregions = num_regions,
  length.mean = mean_length,
  length.sd = sd_length,
  genome = regions_without_peaks,
  non.overlapping = TRUE
)

random_regions_genomewide <- createRandomRegions(
  nregions = num_regions,
  length.mean = mean_length,
  length.sd = sd_length,
  genome = genome,
  non.overlapping = TRUE
)

list.files("~/Dropbox/singulomics/github_rda/linkpeak/encode_mm10_histone_mark", pattern = "\\.bigwig", full.names = T) %>% 
  "names<-"(., gsub("/.+/Encode_(.+?)_.+", "\\1", .)) %>% 
  map(function(bigwig_){
    print(bigwig_)
    bw <- import.bw(bigwig_)
    
    levels(gene_peak_cor_df_2$linked_by) %>% 
      "names<-"(.,.) %>% 
      #  .[1] %>% 
      map(function(linked_by_){
        gene_peak_cor_df_2 %>% dplyr::filter(linked_by == linked_by_) -> df_
        data.frame(
          seqnames = gsub("(.+)-(.+)-(.+)", "\\1", df_$peak), 
          start = gsub("(.+)-(.+)-(.+)", "\\2", df_$peak) %>% as.numeric(), 
          end = gsub("(.+)-(.+)-(.+)", "\\3", df_$peak) %>% as.numeric()
        ) %>% GenomicRanges::makeGRangesFromDataFrame() -> peak_gr
        GenomicRanges::findOverlaps(peak_gr, bw) -> overlaps
        
        bw[subjectHits(overlaps)]@ranges@width -> bw_width
        peak_gr[queryHits(overlaps)]@ranges@width -> peak_width
        width_ <- pintersect(peak_gr[queryHits(overlaps)], bw[subjectHits(overlaps)])@ranges@width
        
        signal_values <- mcols(bw)$score[subjectHits(overlaps)]
        peaks_signal <- data.frame(
          peak_id = queryHits(overlaps),
          signal = signal_values,
          overlapping_width = width_,
          bw_width = bw_width, 
          peak_width = peak_width
        )
        
        mean_signal <- peaks_signal %>%
          group_by(peak_id) %>%
          group_by(peak_id, peak_width) %>% 
          group_map(function(x,y){
            df_ = x
            df_ %>% dplyr::mutate(norm_signal = (signal*overlapping_width)) -> df_
            normalized_signal = sum(df_$norm_signal)/df_$peak_width[1]
            data.frame(peak_id = y, mean_signal = normalized_signal)
          }, .keep = T) %>% do.call(rbind, .)
        
      }) -> mean_signal_list
    
    c("Random_peaks") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        sc@assays$ATAC@counts %>% rownames() -> peaks_
        set.seed(123)
        sample(x = peaks_, size = 10000, replace = F) -> random_peaks_
        data.frame(
          seqnames = gsub("(.+)-(.+)-(.+)", "\\1", random_peaks_), 
          start = gsub("(.+)-(.+)-(.+)", "\\2", random_peaks_) %>% as.numeric(), 
          end = gsub("(.+)-(.+)-(.+)", "\\3", random_peaks_) %>% as.numeric()
        ) %>% GenomicRanges::makeGRangesFromDataFrame() -> peak_gr
        GenomicRanges::findOverlaps(peak_gr, bw) -> overlaps
        
        bw[subjectHits(overlaps)]@ranges@width -> bw_width
        peak_gr[queryHits(overlaps)]@ranges@width -> peak_width
        width_ <- pintersect(peak_gr[queryHits(overlaps)], bw[subjectHits(overlaps)])@ranges@width
        
        signal_values <- mcols(bw)$score[subjectHits(overlaps)]
        peaks_signal <- data.frame(
          peak_id = queryHits(overlaps),
          signal = signal_values,
          overlapping_width = width_,
          bw_width = bw_width, 
          peak_width = peak_width
        )
        
        mean_signal <- peaks_signal %>%
          group_by(peak_id) %>%
          group_by(peak_id, peak_width) %>% 
          group_map(function(x,y){
            df_ = x
            df_ %>% dplyr::mutate(norm_signal = (signal*overlapping_width)) -> df_
            normalized_signal = sum(df_$norm_signal)/df_$peak_width[1]
            data.frame(peak_id = y, mean_signal = normalized_signal)
          }, .keep = T) %>% do.call(rbind, .)
        
      }) -> random_mean_signal_list
    
    c("Random_region_genomewide") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        random_regions_genomewide -> peak_gr
        GenomicRanges::findOverlaps(peak_gr, bw) -> overlaps
        
        bw[subjectHits(overlaps)]@ranges@width -> bw_width
        peak_gr[queryHits(overlaps)]@ranges@width -> peak_width
        width_ <- pintersect(peak_gr[queryHits(overlaps)], bw[subjectHits(overlaps)])@ranges@width
        
        signal_values <- mcols(bw)$score[subjectHits(overlaps)]
        peaks_signal <- data.frame(
          peak_id = queryHits(overlaps),
          signal = signal_values,
          overlapping_width = width_,
          bw_width = bw_width, 
          peak_width = peak_width
        )
        
        mean_signal <- peaks_signal %>%
          group_by(peak_id) %>%
          group_by(peak_id, peak_width) %>% 
          group_map(function(x,y){
            df_ = x
            df_ %>% dplyr::mutate(norm_signal = (signal*overlapping_width)) -> df_
            normalized_signal = sum(df_$norm_signal)/df_$peak_width[1]
            data.frame(peak_id = y, mean_signal = normalized_signal)
          }, .keep = T) %>% do.call(rbind, .)
        
      }) -> random_regions_genomewide_list
    
    c("Random_region_non_atac_peaks") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        random_regions_non_atac_peaks -> peak_gr
        GenomicRanges::findOverlaps(peak_gr, bw) -> overlaps
        
        bw[subjectHits(overlaps)]@ranges@width -> bw_width
        peak_gr[queryHits(overlaps)]@ranges@width -> peak_width
        width_ <- pintersect(peak_gr[queryHits(overlaps)], bw[subjectHits(overlaps)])@ranges@width
        
        signal_values <- mcols(bw)$score[subjectHits(overlaps)]
        peaks_signal <- data.frame(
          peak_id = queryHits(overlaps),
          signal = signal_values,
          overlapping_width = width_,
          bw_width = bw_width, 
          peak_width = peak_width
        )
        
        mean_signal <- peaks_signal %>%
          group_by(peak_id) %>%
          group_by(peak_id, peak_width) %>% 
          group_map(function(x,y){
            df_ = x
            df_ %>% dplyr::mutate(norm_signal = (signal*overlapping_width)) -> df_
            normalized_signal = sum(df_$norm_signal)/df_$peak_width[1]
            data.frame(peak_id = y, mean_signal = normalized_signal)
          }, .keep = T) %>% do.call(rbind, .)
        
      }) -> random_regions_non_atac_peaks_list
    
    
    c(mean_signal_list, random_mean_signal_list, random_regions_genomewide_list, random_regions_non_atac_peaks_list) -> mean_signal_list_1
    names(mean_signal_list_1)
    
    mean_signal_list_1 %>% 
      map2(.x=.,.y=names(.),.f=function(df_, linked_by){
        df_ %>% dplyr::mutate(linked_by = linked_by) %>% 
          dplyr::select(linked_by, mean_signal)
      }) %>% do.call(rbind, .) -> df_bw_signal
    
    gc()
    return(df_bw_signal)
    
  }) -> df_bw_signal_list
gc()

df_bw_signal_list %>% 
  map(function(df_){
    df_ %>% dplyr::mutate(linked_by = case_when(
      linked_by == "Random_peaks" ~ "Random peaks\n(ATAC)",
      linked_by == "Random_region_genomewide" ~ "Random region\n(Genomewide)",
      linked_by == "Random_region_non_atac_peaks" ~ "Random region\n(non ATAC peaks)",
      TRUE ~ linked_by
    ))
  }) -> df_bw_signal_list

df_bw_signal_list %>% 
  map(function(df_){
    levels(factor(df_$linked_by)) %>% .[c(2,3,1,4,5,6)] -> levels_
    df_ %>% dplyr::mutate(linked_by = factor(linked_by, levels = levels_)) -> df_
  }) -> df_bw_signal_list

df_bw_signal_list %>% 
  map(function(df_){
    levels(df_$linked_by)[c(1,4:6)] -> levels_
    df_ %>% dplyr::filter(linked_by %in% levels_) -> df_
    factor(df_$linked_by, levels = levels_) -> df_$linked_by
    return(df_)
  }) -> df_bw_signal_list

list(
  H3K9ac = df_bw_signal_list$H3K9ac %>% 
    ggplot(aes(x = linked_by, y = mean_signal, fill = linked_by)) + 
    geom_boxplot() + 
#    scale_y_break(breaks = c(4,20), scales = 0.5) + 
    scale_y_break(breaks = c(3,20), scales = 0.5) + 
    theme_classic() + 
    ylab("H3K9ac signal") + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + 
    xlab(NULL),
  
  H3K27ac = df_bw_signal_list$H3K27ac %>% 
    ggplot(aes(x = linked_by, y = mean_signal, fill = linked_by)) + 
    geom_boxplot() + 
#    scale_y_break(breaks = c(7,20), scales = 0.5) + 
    scale_y_break(breaks = c(6,10), scales = 0.5) + 
    theme_classic() + 
    ylab("H3K27ac signal") + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + 
    xlab(NULL),
  
  H3K4me1 = df_bw_signal_list$H3K4me1 %>% 
    ggplot(aes(x = linked_by, y = mean_signal, fill = linked_by)) + 
    geom_boxplot() + 
#    scale_y_break(breaks = c(4,20), scales = 0.5) + 
#    scale_y_break(breaks = c(3,6), scales = 0.5) + 
    scale_y_break(breaks = c(2.5,6), scales = 0.5) + 
    theme_classic() + 
    ylab("H3K4me1 signal") + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + 
    xlab(NULL),
  
  H3K4me3 = df_bw_signal_list$H3K4me3 %>% 
    ggplot(aes(x = linked_by, y = mean_signal, fill = linked_by)) + 
    geom_boxplot() + 
#    scale_y_break(breaks = c(6,20), scales = 0.5) + 
    scale_y_break(breaks = c(4.5,20), scales = 0.5) + 
    theme_classic() + 
    ylab("H3K4me3 signal") + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + 
    xlab(NULL),
  
  H3K27me3 = df_bw_signal_list$H3K27me3 %>% 
    ggplot(aes(x = linked_by, y = mean_signal, fill = linked_by)) + 
    geom_boxplot() + 
#    scale_y_break(breaks = c(4,20), scales = 0.5) + 
    scale_y_break(breaks = c(1,5), scales = 0.5) + 
#    scale_y_break(breaks = c(0.25,2), scales = 0.5) + 
    theme_classic() + 
    ylab("H3K27me3 signal")+
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) + 
    xlab(NULL)
) -> p_histone_mark

p_histone_mark$H3K27ac+p_histone_mark$H3K4me1 # Fig 4I
####