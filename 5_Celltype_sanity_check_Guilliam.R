library(tidyverse)
library(Seurat)
library(Signac)
library(pheatmap)

setwd("~/Dropbox/singulomics")
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc=sc[,!grepl('KO',sc$group)]

# 1. Read celltype specific DEG from Guilliam's data -----
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_Hepatocyte_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> Hepatocyte_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_Fibroblast_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> Fibroblast_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_endothelial_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> endothelial_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_KC_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> KC_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_Cholangiocyte_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> Cholangiocyte_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_Tcells_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> Tcells_DEGs
read.csv(file = "./rda/Other_GES_data/Guilliams_2022_Cell_supp/mouse_Bcells_DEGs.csv",
         header = T, sep = ",", stringsAsFactors = F) -> Bcells_DEGs

DEG_list = list()
ls() %>% .[matches("DEGs", vars = .)] %>% 
  map(function(x){
    DEG_list[[x]] <<- get(x = x)
  })
names(DEG_list)

sc$celltype %>% table()

names(DEG_list) <- c("B cells", "Cholangiocytes", "Endothelial cells", "Fibroblasts", 
                     "Hepatocytes", "Kupffer cells", "T cells")
DEG_list = DEG_list[c("Hepatocytes", "Endothelial cells", "Fibroblasts", "Kupffer cells", "B cells", "T cells", "Cholangiocytes")]

####

# 2. Draw celltype specifc heatmap from identified marker genes ----
map2(DEG_list, names(DEG_list), function(x,y){
  x %>% filter(proba_not_de < 0.05) %>% .$X -> DEG
  intersect(DEG, rownames(sc@assays$SCT@data)) -> intersected_genes
  
  print(sprintf("%s: %s", y, length(DEG)))
  print(sprintf("RNA: %s", nrow(sc@assays$SCT@data)))
  print(sprintf("intersected genes: %s", length(intersected_genes)))
  
  unique(sc$celltype) %>% 
    map(function(x_){
      sc@meta.data %>% filter(celltype == x_) %>% rownames() -> cell_sele     
      sc@assays$SCT@data[intersected_genes, cell_sele] %>% rowMeans() %>% 
        as.data.frame() %>% "colnames<-"(., x_)
    }) %>% do.call(cbind, .) %>% as.matrix() -> df_
  
  df_ %>%
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() -> df_
  
#  pheatmap(df_, cluster_rows = T, cluster_cols = F, show_rownames = F, 
#           main=sprintf("%s DEG: n = %s", y, length(intersected_genes))) -> p
  pheatmap(df_, cluster_rows = F, cluster_cols = F, show_rownames = F, 
           main=sprintf("%s DEG: n = %s", y, length(intersected_genes))) -> p
  p$gtable -> p
  #  recordPlot() -> p
  return(p)
}) -> p_list

cowplot::plot_grid(plotlist = p_list, ncol = 4, nrow = 2) #Supp_Fig_5C ----

# 3. Find overlapping celltype specific DEGs (Our data vs Guilliam's data) ----
Idents(sc) %>% as.data.frame() %>% 
  rownames_to_column("cell_name") %>%
  "colnames<-"(., c("cell_name", "old")) %>% 
  left_join(x = ., 
            y = sc@meta.data %>% rownames_to_column("cell_name"), 
            by = "cell_name") %>% 
  dplyr::select(cell_name, celltype) %>% 
  {setNames(.$celltype, as.character(.$cell_name))} -> new_ident

all(names(Idents(sc)) == names(new_ident))
Idents(sc) <- new_ident

DefaultAssay(sc) <- "SCT"
sc@meta.data$celltype %>% unique() %>% 
  "names<-"(., .) %>% 
  map(function(x){
    print(x)
    FindMarkers(sc, ident.1 = x)
  }) -> sc_markers_list

list.files(path = "./rda/Other_GES_data/Guilliams_2022_Cell_supp", pattern = "DEGs", full.names = T, recursive = T) %>% 
  "names<-"(., gsub(".+/mouse_(.+?)_.+", "\\1", .)) %>% 
  #  .[1] %>% 
  map(function(x){
    read.csv(x, header = T, stringsAsFactors = F) -> df_
  }) -> celltype_DEGs
celltype_DEGs[c("Hepatocyte", "endothelial", "Fibroblast", "KC", "Bcells", "Tcells", "Cholangiocyte")] -> celltype_DEGs
names(celltype_DEGs) <- names(sc_markers_list)

ggvenn_list_ = list()
names(sc_markers_list) %>% "names<-"(., .) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype_ = y
    print(y)
    sc_markers_list[[x]] %>% 
      #      filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) -> DEG_scRNA
      dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.58) -> DEG_scRNA
#      filter(p_val_adj < 0.05) -> DEG_scRNA
    
    celltype_DEGs[[x]] %>%
      #    filter(abs(lfc_mean) > 0.5, proba_not_de < 0.05) -> DEG_Guilliams
      dplyr::filter(proba_not_de < 0.05) -> DEG_Guilliams
    DEG_Guilliams %>% dplyr::filter(X %in% rownames(sc@assays$SCT@data)) -> DEG_Guilliams
    
    sprintf("celltype: %s, scRNA: %s, Guilliams: %s", x, nrow(DEG_scRNA), nrow(DEG_Guilliams)) %>% print()
    
    list(
      sc_RNA = rownames(DEG_scRNA),
      DEG_Guilliams = DEG_Guilliams$X
    ) -> ggvenn_list
    ggvenn_list_[[celltype_]] <<- ggvenn_list
    
    ggvenn::ggvenn(ggvenn_list) + 
      ggtitle(x) + 
      theme(plot.title = element_text(hjust = 0.5))
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, ncol = 4, nrow = 2, common.legend = T, legend = "top")
####

# 4. Plot pair-wise expression of DEGs (our data vs Guillam's) ----

#Read Gulliam's seurat data
library(Seurat)
library(tidyverse)

load('~/Dropbox/singulomics/rda/reference_sc_data/LCA/lca_nuc.rda')

expression_matrix <- ReadMtx(
  mtx = "~/Downloads/Liver_atlas/rawData_mouseStSt/countTable_mouseStSt/matrix.mtx.gz", 
  features = "~/Downloads/Liver_atlas/rawData_mouseStSt/countTable_mouseStSt/features.tsv.gz",
  cells = "~/Downloads/Liver_atlas/rawData_mouseStSt/countTable_mouseStSt/barcodes.tsv.gz", 
  feature.column = 1
)

Guilliam_sc <- CreateSeuratObject(counts = expression_matrix)
rm(expression_matrix)
View(Guilliam_sc)

length(lca_nuc$annot)
#18666
lca_nuc$annot %>% names() %in% colnames(Guilliam_sc) %>% sum()
#18666

lca_nuc$annot %>% table()
c("Hepatocytes", "Endothelial cells", "Fibroblasts", "Kupffer cells", "B cells", "T cells", "Cholangiocytes") %>% 
#  .[1] %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    celltype_ = x
    lca_nuc$annot %>% {.[.==celltype_]} %>% names() -> cellnames_
    Guilliam_sc[, cellnames_] -> sc
  }) -> Guilliam_sc_list
rm(Guilliam_sc)

Guilliam_sc_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    sc = x
    print(y)
#    ncol(sc)
    sc@assays$RNA@counts %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) -> RNA_norm
    RNA_norm %>% rowMeans() -> mean_expr
  }) -> Guilliam_mean_expr
gc()

save(Guilliam_sc_list, Guilliam_mean_expr, file = "~/Dropbox/singulomics/github_rda/Guilliam_mean_expr.rda")
####

#Read our scRNA data ----
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(across(everything(), function(x){(x/sum(x))*median(sc$nCount_RNA)})) -> sc_RNA_nrom

ggvenn_list_ %>% 
#  .[1] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype_ = y
#    intersect(x$sc_RNA, x$DEG_Guilliams) -> intersected_genes
#    length(intersected_genes)
    x$DEG_Guilliams -> marker_genes_
    
    cellnames_ = sc@meta.data %>% dplyr::filter(celltype == celltype_) %>% rownames()
    sc_RNA_nrom[marker_genes_, cellnames_] %>% rowMeans() %>% 
      as.data.frame() %>% "colnames<-"(., "scRNA") %>% 
      rownames_to_column("gene") -> sc_df_
    
    Guilliam_mean_expr[[celltype_]] %>% as.data.frame() %>% 
      "colnames<-"(., "Guilliam") %>% 
      rownames_to_column("gene") -> Guilliam_df_
    Guilliam_df_ %>% 
      dplyr::filter(gene %in% marker_genes_) -> Guilliam_df_
    print(sprintf("nrow_sc: %s, nrow_Guilliam_sc: %s", nrow(sc_df_), nrow(Guilliam_df_)))
    full_join(sc_df_, Guilliam_df_, by = "gene") -> df_
  }) -> list_tmp

list_tmp %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype_ = y
    n_genes = nrow(x)
    df_ = x
#    print(head(df_))
    cor(df_$scRNA, df_$Guilliam, method = "pearson") %>% round(.,2) -> cor_
    df_ %>% 
      ggplot(aes(x = scRNA, y = Guilliam)) + 
      geom_point() + 
      theme_classic() + 
      ylab("Guilliam (Normalized RNA expr.)") + 
      xlab("scRNA (Normalized RNA expr.)") -> p
    p + ggtitle(sprintf("%s: n=%s\nr=%s", celltype_, n_genes, cor_)) -> p
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 4) #Supp_Fig_5b ----

#Plot celltype specific heatmap
rm(list=ls())

library(tidyverse)
library(Seurat)
library(Signac)
library(pheatmap)

setwd("~/Dropbox/singulomics")
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc=sc[,!grepl('KO',sc$group)]

Idents(sc) %>% as.data.frame() %>% 
  rownames_to_column("cell_name") %>%
  "colnames<-"(., c("cell_name", "old")) %>% 
  left_join(x = ., 
            y = sc@meta.data %>% rownames_to_column("cell_name"), 
            by = "cell_name") %>% 
  dplyr::select(cell_name, celltype) %>% 
  {setNames(.$celltype, as.character(.$cell_name))} -> new_ident

all(names(Idents(sc)) == names(new_ident))
Idents(sc) <- new_ident

DefaultAssay(sc) <- "SCT"
sc@meta.data$celltype %>% unique() %>% 
  "names<-"(., .) %>% 
  map(function(x){
    print(x)
    FindMarkers(sc, ident.1 = x)
  }) -> sc_markers_list

list_tmp = list()
sc_markers_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype_ = y
    df_ = x
    df_ %>% filter(p_val_adj < 0.05, avg_log2FC > 0.58) %>% rownames() -> DEG
    print(sprintf("%s DEG: %s", celltype_, length(DEG)))
    
    unique(sc$celltype) %>% 
      map(function(x_){
        sc@meta.data %>% filter(celltype == x_) %>% rownames() -> cell_sele     
        sc@assays$SCT@data[DEG, cell_sele] %>% rowMeans() %>% 
          as.data.frame() %>% "colnames<-"(., x_)
      }) %>% do.call(cbind, .) %>% as.matrix() -> df_
    
    df_ %>%
      apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
      t() -> df_
    
#    pheatmap(df_, cluster_rows = T, cluster_cols = F, show_rownames = F, 
    pheatmap(df_, cluster_rows = T, cluster_cols = F, show_rownames = F, cutree_rows = 10,
             main=sprintf("%s DEG: n = %s", y, length(DEG))) -> p
    list_tmp[[celltype_]] <<- p
    p$gtable -> p
    #  recordPlot() -> p
    return(p)
  }) -> p_list

cowplot::plot_grid(plotlist = p_list, ncol = 4, nrow = 2)

list_tmp %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype_ = y
    x$tree_row %>% stats::cutree(., 10) %>% table() %>% sort() %>% rev() -> table_
    names(table_)[1] -> clust_
    x$tree_row %>% stats::cutree(., 10) %>% {.[. == clust_]} %>% names() -> genes_
    sc_markers_list[[celltype_]][genes_, ] %>% dplyr::arrange(p_val_adj, desc(avg_log2FC)) -> df_
    genes_ = rownames(df_)
    if (celltype_ == "Kupffer cells"){
      x$tree_row %>% stats::cutree(., 10) %>% {.[. == clust_]} %>% names() %>% .[2:21] -> genes_
    }else if (celltype_ == "Endothelial cells"){
      genes_[c(1:9, 11:21)] -> genes_
    } else if (celltype_ == "Fibroblasts"){
      genes_[c(1:2,4:8,10:22)] -> genes_
    } else if (celltype_ == "T cells"){
      genes_ = genes_[1:20]
    }else{
      genes_ = genes_[1:20]
    }
    return(genes_)
  }) %>% purrr::reduce(., c) -> features_
length(features_)

sc@assays$SCT@data %>% 
#  .[1:10,1:10] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() -> sc@assays$SCT@scale.data
sc$celltype = factor(sc$celltype, levels = c("Hepatocytes", "Endothelial cells", "Fibroblasts", "Kupffer cells", "B cells", "T cells", "Cholangiocytes"))
sc$celltype

DoHeatmap(sc, features = features_, group.by = "celltype", assay = "SCT", slot = "scale.data", label = T) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                  "RdYlBu")))(100)) + 
  theme(axis.text.y = element_blank()) #Supp Fig 5A
rm(sc)
gc()