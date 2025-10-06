library(Seurat)
library(Signac)
library(tidyverse)

readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells

colnames(sc) %>% length()
length(hepatocytes_cells)

sc <- sc[,hepatocytes_cells]
sc$celltype %>% table()

DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "multicca.umap", repel = TRUE)
DefaultAssay(sc) <- "SCT"
sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=0.2)
DimPlot(sc, group.by = "seurat_clusters", label = TRUE, reduction = "multicca.umap", repel = TRUE)
sc@meta.data$seurat_clusters %>% unique() %>% length()

sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=1)
DimPlot(sc, group.by = "seurat_clusters", label = TRUE, reduction = "multicca.umap", repel = TRUE) -> p2

# 1. Generate input data for cyclic genes detection ----
## 1.1 RNA ----
sc$nCount_RNA %>% median() -> RNA_median
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*RNA_median})) %>% 
  as.matrix() -> RNA_norm
rownames(RNA_norm) -> Gene_RNA_assay

seq(0.5, 30, 0.5) %>% "names<-"(., sprintf("res_%s", .)) %>% 
  #  .[1] %>%
  map2(.x=.,.y=names(.),.f=function(x,y){
    print(x)
    res_ = y
    sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                       resolution=x)
    sc@meta.data %>% rownames_to_column("cell_name") %>% 
      select(cell_name, seurat_clusters, ZT) -> meta_data
    
    meta_data %>% group_by(ZT, seurat_clusters) %>% 
      group_map(function(x,y){
        x$cell_name -> cell_name_
        y$ZT ->  ZT_
        y$seurat_clusters -> meta_cells_
        
        if (length(cell_name_) > 1){
          RNA_norm[,cell_name_] %>% rowMeans() -> gene_exp
        }else{
          RNA_norm[,cell_name_] -> gene_exp
        }
        
        gene_exp %>% as.data.frame() %>% rownames_to_column("Gene") %>% 
          "colnames<-"(.,c("Gene", sprintf("%s_REP%s", ZT_, as.numeric(as.character(meta_cells_))+1))) -> gene_exp
        
      }, .keep = T) %>% 
      
      #      do.call(cbind, .) -> df_
      purrr::reduce(., left_join, by = "Gene") -> df_
      write.csv(df_, file = sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_RNA.csv", res_), col.names = T, row.names = F, quote = F)  
      print(sprintf("Output: %s", sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_RNA.csv", res_)))
      rm(df_)
    gc()
  }) -> list_tmp

rm(list_tmp)
rm(RNA_norm)

## 1.2 ATAC ----
sc$nCount_ATAC %>% median() -> ATAC_median

sc@assays$ATAC_TSS@counts %>% 
  as.data.frame() %>%
  mutate(across(everything(), function(x){(x/sum(x))*ATAC_median})) %>% 
  as.matrix() -> ATAC_norm
rownames(ATAC_norm) -> Gene_ATAC_assay

seq(0.5, 30, 0.5) %>% "names<-"(., sprintf("res_%s", .)) %>% 
  #  .[1] %>%
  map2(.x=.,.y=names(.),.f=function(x,y){
    print(x)
    res_ = y
    sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                       resolution=x)
    sc@meta.data %>% rownames_to_column("cell_name") %>% 
      select(cell_name, seurat_clusters, ZT) -> meta_data
    
    meta_data %>% group_by(ZT, seurat_clusters) %>% 
      group_map(function(x,y){
        x$cell_name -> cell_name_
        y$ZT ->  ZT_
        y$seurat_clusters -> meta_cells_
        
        if (length(cell_name_) > 1){
          ATAC_norm[,cell_name_] %>% rowMeans() -> gene_exp
        }else{
          ATAC_norm[,cell_name_] -> gene_exp
        }
        
        gene_exp %>% as.data.frame() %>% rownames_to_column("Gene") %>% 
          "colnames<-"(.,c("Gene", sprintf("%s_REP%s", ZT_, as.numeric(as.character(meta_cells_))+1))) -> gene_exp
        
      }, .keep = T) %>% 
      
      #      do.call(cbind, .) -> df_
      purrr::reduce(., left_join, by = "Gene") -> df_
      write.csv(df_, file = sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_ATAC.csv", res_), col.names = T, row.names = F, quote = F)  
      print(sprintf("Output: %s", sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_ATAC.csv", res_)))
      rm(df_)
    gc()
  }) -> list_tmp

rm(ATAC_norm)
rm(list_tmp)

## 1.3 Gene activity ----
sc@assays[["ATAC"]]@fragments[[1]]@path <- "~/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
DefaultAssay(sc) <- "ATAC"
#gene.activities <- GeneActivity(sc, extend.upstream = 10000, extend.downstream = 2000)
gene.activities <- GeneActivity(sc, extend.upstream = 2000, biotypes = NULL)

gene.activities %>% colSums() %>% median() -> gene_activity_median

gene.activities %>% 
  as.data.frame() %>%
  mutate(across(everything(), function(x){(x/sum(x))*gene_activity_median})) %>% 
  as.matrix() -> gene_activity_norm
rm(gene.activities)
gc()

rownames(gene_activity_norm) -> Gene_gene_activity_assay

seq(0.5, 30, 0.5) %>% "names<-"(., sprintf("res_%s", .)) %>% 
  #  .[1] %>%
  map2(.x=.,.y=names(.),.f=function(x,y){
    print(x)
    res_ = y
    sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                       resolution=x)
    sc@meta.data %>% rownames_to_column("cell_name") %>% 
      select(cell_name, seurat_clusters, ZT) -> meta_data
    
    meta_data %>% group_by(ZT, seurat_clusters) %>% 
      group_map(function(x,y){
        x$cell_name -> cell_name_
        y$ZT ->  ZT_
        y$seurat_clusters -> meta_cells_
        
        if (length(cell_name_) > 1){
          gene_activity_norm[,cell_name_] %>% rowMeans() -> gene_exp
        }else{
          gene_activity_norm[,cell_name_] -> gene_exp
        }
        
        gene_exp %>% as.data.frame() %>% rownames_to_column("Gene") %>% 
          "colnames<-"(.,c("Gene", sprintf("%s_REP%s", ZT_, as.numeric(as.character(meta_cells_))+1))) -> gene_exp
        
      }, .keep = T) %>% 
      
      #      do.call(cbind, .) -> df_
      purrr::reduce(., left_join, by = "Gene") -> df_
    write.csv(df_, file = sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_gene_activity.csv", res_), col.names = T, row.names = F, quote = F)  
    print(sprintf("Output: %s", sprintf("~/Dropbox/singulomics/github_rda/output/Hepatocytes/Raw_data/Hep_%s_gene_activity.csv", res_)))
    rm(df_)
    gc()
  }) -> list_tmp

rm(gene_activity_norm)
rm(list_tmp)