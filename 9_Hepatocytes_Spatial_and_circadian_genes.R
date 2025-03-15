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