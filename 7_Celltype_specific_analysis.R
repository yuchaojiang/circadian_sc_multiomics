# 1. Downsample Hepatocytes, Endothelial cells, Fibroblasts and Kupffer cells ----
library(Seurat)
library(Signac)
library(tidyverse)
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc=sc[,!grepl('KO',sc$group)]

# Generate normalized RNA expression matrix
sc$nCount_RNA %>% median() -> RNA_median
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*RNA_median})) -> RNA_norm

# Generate motif score matrix
sc@assays$MOTIF@data %>% as.data.frame() -> MOTIF_norm

# Generate ATAC activity matrix
sc@assays[["ATAC"]]@fragments[[1]]@path <- "~/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
DefaultAssay(sc) <- "ATAC"
gene.activities <- GeneActivity(sc, extend.upstream = 2000, biotypes = NULL)

gene.activities %>% colSums() %>% median() -> gene_activity_median

gene.activities %>% 
  as.data.frame() %>%
  mutate(across(everything(), function(x){(x/sum(x))*gene_activity_median})) %>% 
  as.matrix() -> gene_activity_norm
gene_activity_norm %>% as.data.frame() -> gene_activity_norm
rm(gene.activities)
gc()

#Extract meta data
sc@meta.data -> sc_meta
rm(sc)
gc()

# Generate ZT_min
list.files(path = "~/Dropbox/singulomics/github_rda", pattern = "cellnames\\.rds", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_cellnames\\.rds", "\\1", .)) %>% 
  .[c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells")] %>% 
  map(function(x){
    readRDS(x) -> cellnames
  }) %>% purrr::reduce(., c) %>% 
  {sc_meta[., ]} %>% 
  {table(.$ZT, .$celltype)} %>% 
  apply(1, min) -> ZT_min

#Write downsampled matrix files (RNA, ATAC_activity, motif score)
list_cell_names = list()
list.files(path = "~/Dropbox/singulomics/github_rda", pattern = "cellnames\\.rds", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_cellnames\\.rds", "\\1", .)) %>% 
  .[c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells")] %>%
  map2(.x=.,.y=names(.),.f=function(x,y){
    print(x)
    readRDS(x) -> cellnames
    y -> celltype
    
    sc_meta %>% as.data.frame() %>% 
      .[cellnames, ] %>% .$ZT %>% table() -> ZT_table
    print(celltype)
    print(ZT_table)
    
    c(seq(10,200,10)) %>% 
      "names<-"(., sprintf("seed_%s", .)) %>% 
      map(function(seed_){
        set.seed(seed_)
        sc_meta %>% as.data.frame() %>% 
          .[cellnames, ] -> df_
        df_ %>% rownames_to_column("Cell") %>% 
          group_by(ZT) %>%
          group_map(function(x,y){
            df_ = x
            ZT_ = y$ZT %>% as.character()
            n = ZT_min[ZT_]
            
            set.seed(seed_) 
            sample(df_$Cell, n, replace = F) -> cellnames
            list_cell_names[[sprintf("seed_%s", seed_)]][[celltype]][[ZT_]] <<- cellnames
            
            RNA_norm[, cellnames] %>% "colnames<-"(., sprintf("%s_REP%s", ZT_, 1:ncol(.))) -> RNA_norm
            gene_activity_norm[, cellnames] %>% "colnames<-"(., sprintf("%s_REP%s", ZT_, 1:ncol(.))) -> gene_activity_norm
            MOTIF_norm[, cellnames] %>% "colnames<-"(., sprintf("%s_REP%s", ZT_, 1:ncol(.))) -> MOTIF_norm
            
            list_ = list(RNA = RNA_norm, gene_activity = gene_activity_norm, MOTIF = MOTIF_norm)
            return(list_)
            
          }, .keep = T) -> list_1
        
        c("RNA", "gene_activity", "MOTIF") %>% "names<-"(.,.) %>% 
          map(function(x){
            list_1 %>% 
              map(function(x_){
                x_[[x]] -> df_
              }) %>% purrr::reduce(., cbind) -> df_
            return(df_)
          })
      })
  }) -> list_tmp

saveRDS(list_cell_names, file = "~/Dropbox/singulomics/github_rda/output/Celltype_specific/list_cell_names.rds")

#Pseudo-bulk the downsampled data to 10 meta cells
c(3, 4, 5, 10) %>% 
  "names<-"(., sprintf("replicate_%s", .)) %>% 
  .[4] %>% 
  map(function(x){
    replicate_ = x
    list_tmp %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        celltype_ = y
        x %>% 
          map2(.x=.,.y=names(.),.f=function(x,y){
            seed_1 = y
            list_ = x
            c(1,2,3) %>% 
              .[1] %>%
              "names<-"(., sprintf("seed_%s", .)) %>% 
              map2(.x=.,.y=names(.),.f=function(x,y){
                seed_2 = y
                set.seed(x)
                list_ %>% 
                  map2(.x=.,.y=names(.),.f=function(x,y){
                    type_ = y
                    df_ = x
                    c("ZT02", "ZT06", "ZT10", "ZT14", "ZT18", "ZT22") %>% 
                      "names<-"(.,.) %>% 
                      map(function(x){
                        ZT_ = x
                        idx_ = grep(x, colnames(df_))
                        n = round(length(idx_)/replicate_)
                        c(1:replicate_) %>% 
                          map(function(x){
                            if (x != replicate_){
                              rep_ = c(rep(x, n))
                            }else{
                              rep_ = c(rep(x, length(idx_)-((x-1)*n)))
                            }
                          }) %>% purrr::reduce(., c) -> rep_
                        data.frame(idx = idx_, rep = rep_) %>% 
                          mutate(rep = sample(rep, length(rep), replace = F)) %>% 
                          group_by(rep) %>% 
                          group_map(function(x,y){
                            x$idx -> idx_
                            rep_ = y$rep
                            df_[, idx_] %>% rowMeans() %>% 
                              as.data.frame() %>% 
                              "colnames<-"(., sprintf("%s_REP%s", ZT_, rep_))
                          }, .keep = T) %>% purrr::reduce(., cbind)
                      }) %>% purrr::reduce(., cbind) -> df_1
                  })
              })
          })
      }) -> list_tmp_1
  }) -> list_tmp_1

# Output the downsampled data to 10 meta cells
list_tmp_1 %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    replicate_ = y
    x %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        celltype_ = y
        x %>% 
          map2(.x=.,.y=names(.),.f=function(x,y){
            seed_1 = y
            x %>% 
              map2(.x=.,.y=names(.),.f=function(x,y){
                seed_2 = y
                x %>% 
                  map2(.x=.,.y=names(.),.f=function(x,y){
                    type_ = y
                    df_ = x %>% rownames_to_column("Gene")
                    output_dir = sprintf("~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_data_downsampled/%s/", replicate_)
                    output_name = sprintf("%s%s_%s_%s_%s.csv", output_dir, celltype_, seed_1, gsub("seed_(\\d+)", "\\1", seed_2), type_)
                    write.csv(df_, file = output_name, row.names = F, col.names = T, quote = F)
                    print(output_name)
                  })
              })
          })
      })    
  })
# Output the downsampled data (without pseudo bulk)
list_tmp %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    celltype = y
    x %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        seed_ = y
        x %>% 
          map2(.x=.,.y=names(.),.f=function(x,y){
            df_ = x
            type_ = y
            output_dir = "~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_data"
            output_name = sprintf("%s_%s_%s.csv", celltype, seed_, type_)
            output_name = sprintf("%s/%s", output_dir, output_name)
            df_ %>% rownames_to_column("Gene") -> df_
            write.csv(df_, file = output_name, row.names = F, col.names = T, quote = F)
            print(output_name)
          })
      })
  })
####

# 2. Detect celltype specfic rhythmic genes by JTKcycle and HR ----
source("~/Dropbox/singulomics/github/Calculate_HMP.R")
c("RNA", "gene_activity", "MOTIF") %>% 
  "names<-"(.,.) %>% 
#  .[3] %>% 
  map(function(assay_){
    c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells") %>% 
      "names<-"(.,.) %>% 
      map(function(celltype_){
        list.files("~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_data_downsampled/replicate_10", full.names = T) %>% 
          {.[grepl(assay_, .)]} %>% {.[grepl(celltype_, .)]} %>% gtools::mixedsort() -> files_       
        seeds_ = c()
        files_ %>% 
          map(function(file_){
            seed_ = gsub("/.+/.+_(seed_\\d+?)_.+", "\\1", file_)
            seeds_ <<- c(seeds_, seed_)
            read.csv(file_, header = T, stringsAsFactors = F) -> df_
            Gene_ = df_$Gene
            Timepoint_ = colnames(df_)[-1] %>% gsub("ZT(\\d+)_.+", "\\1", .) %>% as.integer()
            df_ %>% dplyr::select(-Gene) %>% as_tibble() -> df_
            
            cyclic_HMP(nCores_ = 6, exp_matrix = df_, gene = Gene_, timepoints = Timepoint_, minper_ = 20) -> res_
            res_ %>% dplyr::select(matches("Gene|JTK|HR|rAMP")) -> res_
            recal_cauchy_p_and_hmp(res_) -> res_
          }) -> list_
        names(list_) = seeds_
        return(list_)
      })
  }) -> res_
gc()

#write the output from JTKcycle and HR
res_ %>% 
  map2(.x = ., .y = names(.), .f=function(list_, assay_){
    list_ %>% 
      map2(.x = ., .y = names(.), .f=function(list_, celltype_){
        list_ %>% 
        map2(.x = ., .y = names(.), .f=function(df_, seed_){
          dir_ = "/Users/chunyiptong/Dropbox/singulomics/github_rda/output/Celltype_specific/Cauchy_replicate_10"
          file_ = sprintf("%s/%s_%s_1_%s_HMP.csv", dir_, celltype_, seed_, assay_)
          write.csv(df_, file = file_, quote = F, col.names = T, row.names = F)
        })
      })
  })
#####

#### Cell partitioning using monocle3
library(monocle3)
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc@meta.data$group %>% unique()
sc=sc[,!grepl('KO',sc$group)]
DefaultAssay(sc) <- "SCT"
DimPlot(sc, group.by = "predicted.id", label = T, reduction = "multicca.umap")
sc[["UMAP"]] = sc[["multicca.umap"]]

sc_1 = as.cell_data_set(sc, assay = "SCT", default.reduction = "UMAP")
sc_1 = cluster_cells(sc_1, resolution = 1e-3)

p1 <- plot_cells(sc_1, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(sc_1, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1|p2

sc_1 <- learn_graph(sc_1, use_partition = TRUE, verbose = FALSE)

plot_cells(sc_1,
           color_cells_by = "partition",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

colnames(sc_1)[partitions(sc_1) == 2] -> endothelial_cell_ids
colData(sc_1)[endothelial_cell_ids, ] %>% as.data.frame() %>%
  filter(celltype == "Endothelial cells") %>% 
  rownames() -> endothelial_cell_ids
saveRDS(endothelial_cell_ids, 
        file = "~/Dropbox/singulomics/github_rda/Endothelial_cells_cellnames.rds")

colnames(sc_1)[partitions(sc_1) == 3] -> fibroblasts_ids
colData(sc_1)[fibroblasts_ids, ] %>% as.data.frame() %>%
  filter(celltype == "Fibroblasts") %>% 
  rownames() -> fibroblasts_ids
saveRDS(fibroblasts_ids, 
        file = "~/Dropbox/singulomics/github_rda/Fibroblasts_cellnames.rds")

colnames(sc_1)[partitions(sc_1) == 4] -> KC_ids
colData(sc_1)[KC_ids, ] %>% as.data.frame() %>%
  filter(celltype == "Kupffer cells") %>% 
  rownames() -> KC_ids
saveRDS(KC_ids, 
        file = "~/Dropbox/singulomics/github_rda/Kupffer_cells_cellnames.rds")

colnames(sc_1)[partitions(sc_1) == 5] -> Tcells_ids
colData(sc_1)[Tcells_ids, ] %>% as.data.frame() %>%
  filter(celltype == "T cells") %>% 
  rownames() -> Tcells_ids
saveRDS(Tcells_ids, 
        file = "~/Dropbox/singulomics/github_rda/T_cells_cellnames.rds")

colnames(sc_1)[partitions(sc_1) == 6] -> Bcells_ids
colData(sc_1)[Bcells_ids, ] %>% as.data.frame() %>%
  filter(celltype == "B cells") %>% 
  rownames() -> Bcells_ids
saveRDS(Bcells_ids, 
        file = "~/Dropbox/singulomics/github_rda/B_cells_cellnames.rds")

colnames(sc_1)[partitions(sc_1) == 7] -> Cholangiocytes_ids
colData(sc_1)[Cholangiocytes_ids, ] %>% as.data.frame() %>%
  filter(celltype == "Cholangiocytes") %>% 
  rownames() -> Cholangiocytes_ids
saveRDS(Cholangiocytes_ids, 
        file = "~/Dropbox/singulomics/github_rda/Cholangiocytes_cellnames.rds")

rm(sc_1)
rm(list=ls())
####

# 3. Draw dotplots of RNA expression, ATAC activity and motif score ----
ZT_min = c(ZT02 = 57, ZT06 = 33, ZT10 = 80, ZT14 = 85, ZT18 = 118, ZT22 = 121)

library(tidyverse)
library(Seurat)
library(Signac)

readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> Hep_cells
readRDS("~/Dropbox/singulomics/github_rda/Endothelial_cells_cellnames.rds") -> Endo_cells
readRDS("~/Dropbox/singulomics/github_rda/Fibroblasts_cellnames.rds") -> Fibro_cells
readRDS("~/Dropbox/singulomics/github_rda/Kupffer_cells_cellnames.rds") -> Kupffer_cells

sc[, c(Hep_cells, Endo_cells, Fibro_cells, Kupffer_cells)] -> sc

c("Hepatocytes", "Endothelial cells", "Fibroblasts", "Kupffer cells") %>% 
  map(function(x){
    sprintf("%s_ZT%02d", x, seq(02,22,4))
  }) %>% unlist() %>% gsub(" ", "_", .) -> levels_
sc@meta.data %>% mutate(celltype_ZT = sprintf("%s_%s", celltype, ZT) %>% gsub(" ", "_", .) %>% factor(., levels = levels_)) -> sc@meta.data
sc_meta = sc@meta.data
gc()

#Extract normalized RNA expression matrix
sc$nCount_RNA %>% median() -> RNA_median
sc@assays$RNA@counts %>% as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x)*RNA_median)})) %>% 
  as.matrix() -> RNA_norm

#Extract normalized ATAC acitivity matrix
sc$nCount_gene_activity %>% median() -> gene_activity_median
sc@assays$gene_activity@counts %>% as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x)*gene_activity_median)})) %>% 
  as.matrix() -> gene_activity_norm

#Extract motif score matrix
sc@assays$MOTIF@data %>% 
  as.data.frame() -> motif_norm

c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Npas2", "Dbp", "Nfil3", "Nr1d1", "Rorc", 
  "Cry1", "Ciart", "Per1", "Per2") %>% rev() -> se_genes

se_genes %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_ = x
    sc@assays$ATAC@motifs@motif.names %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame() %>%
      rownames_to_column("motif") %>% 
      "colnames<-"(., c("motif", "gene")) %>% 
      dplyr::filter(grepl(sprintf(".*%s.*", gene_), gene, ignore.case = T))
  }) %>% do.call(rbind, .) %>% 
  {
    df_ = .
    motif_ = df_$motif
    gene_ = df_$gene
    motif_norm[motif_, ] -> motif_norm
    rownames(motif_norm) <- gene_
    motif_norm
  } -> motif_norm
####

#Extract cell name for each group (ZT and celltype)
list.files(path = "~/Dropbox/singulomics/github_rda", pattern = "cellnames\\.rds", full.names = T) %>% 
  "names<-"(., gsub("/.+/(.+)_cellnames\\.rds", "\\1", .)) %>% 
  .[c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells")] %>%
  map2(.x=.,.y=names(.),.f=function(x,y){
    print(x)
    readRDS(x) -> cellnames
    y -> celltype
    
    sc_meta %>% as.data.frame() %>% 
      .[cellnames, ] %>% .$ZT %>% table() -> ZT_table
    print(celltype)
    print(ZT_table)
    
    c(seq(10,200,10)) %>% 
      "names<-"(., sprintf("seed_%s", .)) %>% 
      map(function(seed_){
        set.seed(seed_)
        sc_meta %>% as.data.frame() %>% 
          .[cellnames, ] -> df_
        df_ %>% rownames_to_column("Cell") %>% 
          group_by(ZT) %>%
          group_map(function(x,y){
            df_ = x
            ZT_ = y$ZT %>% as.character()
            n = ZT_min[ZT_]
            
            set.seed(seed_) 
            sample(df_$Cell, n, replace = F) -> cellnames
            
          }, .keep = T) %>% 
          purrr::reduce(., c) -> cellnames
      }) -> list_1
  }) -> list_tmp
rm(sc);gc()


c(seq(10,200,10)) %>% 
  "names<-"(., sprintf("seed_%s", .)) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seed_ = y
    list_tmp %>% 
      map(function(x){
        x[[seed_]]
      }) %>% purrr::reduce(., c) -> cellnames
  }) -> list_tmp

  ####

c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Npas2", "Dbp", "Nfil3", "Nr1d1", "Rorc", 
  "Cry1", "Ciart", "Per1", "Per2") %>% rev() -> se_genes

#RNA expression dot plot
list_tmp %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seed_ = y
    cells_ = x
    sc_meta[cells_, ] -> meta_df_
    se_genes %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        gene_ = x
        levels_ %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            celltype_ZT_ = x
            meta_df_ %>% filter(celltype_ZT == celltype_ZT_) %>% rownames() -> cells_
            RNA_norm[gene_, cells_] %>% mean() -> mean_expression
            non_zero_cell = sum(RNA_norm[gene_, cells_] > 0)
            non_zero_proportion = 100*(non_zero_cell/length(cells_))
            data.frame(Gene = gene_, Celltype_ZT = celltype_ZT_, Mean_expression = mean_expression, non_zero_proportion = non_zero_proportion)
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .) -> df_
    
    df_ %>% 
    "rownames<-"(., 1:nrow(.)) %>% 
      mutate(Celltype_ZT = factor(Celltype_ZT, levels = levels_)) -> df_
    
    df_ %>% 
      group_by(Gene) %>%
      group_map(function(x,y){
        df_ = x
        df_ %>% mutate(Scaled_mean_expression = scales::rescale(x = Mean_expression, to = c(0,1))) -> df_
      }, .keep = T) %>% 
      do.call(rbind, .) %>% 
      mutate(Gene = factor(Gene, levels = se_genes)) -> df_
    
      df_ %>% 
        ggplot(aes(x = Celltype_ZT, y = Gene)) + 
        geom_point(aes(size = non_zero_proportion, color = Scaled_mean_expression)) + 
        scale_size(range = c(0, 6), limits = c(min(df_$non_zero_proportion), 100)) + 
        scale_color_gradient(low = "grey", high = "blue") +
        scale_x_discrete(guide = guide_axis(angle = 90)) + 
        theme_classic() -> p_dot_scaled
      
      p_dot_scaled + ggtitle(sprintf("%s", seed_)) -> p_dot_scaled
    
  }) -> list_tmp_1
list_p_dot_RNA_scaled = list_tmp_1 ; rm(list_tmp_1)
####

#ATAC activity dot plot
list_tmp %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seed_ = y
    cells_ = x
    sc_meta[cells_, ] -> meta_df_
    se_genes %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        gene_ = x
        levels_ %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            celltype_ZT_ = x
            meta_df_ %>% filter(celltype_ZT == celltype_ZT_) %>% rownames() -> cells_
            gene_activity_norm[gene_, cells_] %>% mean() -> mean_expression
            non_zero_cell = sum(gene_activity_norm[gene_, cells_] > 0)
            non_zero_proportion = 100*(non_zero_cell/length(cells_))
            data.frame(Gene = gene_, Celltype_ZT = celltype_ZT_, Mean_expression = mean_expression, non_zero_proportion = non_zero_proportion)
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .) -> df_
    
    df_ %>% 
      "rownames<-"(., 1:nrow(.)) %>% 
      mutate(Celltype_ZT = factor(Celltype_ZT, levels = levels_)) -> df_
    
    df_ %>% 
      group_by(Gene) %>%
      group_map(function(x,y){
        df_ = x
        df_ %>% mutate(Scaled_mean_expression = scales::rescale(x = Mean_expression, to = c(0,1))) -> df_
      }, .keep = T) %>% 
      do.call(rbind, .) %>% 
      mutate(Gene = factor(Gene, levels = se_genes)) -> df_
    
    df_ %>% 
      ggplot(aes(x = Celltype_ZT, y = Gene)) + 
      geom_point(aes(size = non_zero_proportion, color = Scaled_mean_expression)) + 
      scale_size(range = c(0, 6), limits = c(min(df_$non_zero_proportion), 100)) + 
      scale_color_gradient(low = "grey", high = "red") +
      scale_x_discrete(guide = guide_axis(angle = 90)) + 
      theme_classic() -> p_dot_scaled
    
    p_dot_scaled + ggtitle(sprintf("%s", seed_)) -> p_dot_scaled
    
  }) -> list_tmp_1
list_p_dot_gene_activity_scaled = list_tmp_1 ; rm(list_tmp_1)
####

#Motif score dot plot
list_tmp %>% 
#  .[9] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    seed_ = y
    print(seed_)
    cells_ = x
    sc_meta[cells_, ] -> meta_df_
#    se_genes %>% 
    rownames(motif_norm) %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(x){
        gene_ = x
#        print(gene_)
        levels_ %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            celltype_ZT_ = x
            meta_df_ %>% filter(celltype_ZT == celltype_ZT_) %>% rownames() -> cells_
            motif_norm[gene_, cells_] %>% unlist() %>% mean() -> mean_expression
            non_zero_cell = sum(motif_norm[gene_, cells_] > 0)
            non_zero_proportion = 100*(non_zero_cell/length(cells_))
            data.frame(Gene = gene_, Celltype_ZT = celltype_ZT_, Mean_expression = mean_expression, non_zero_proportion = non_zero_proportion)
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .) -> df_
    
    df_ %>% 
      "rownames<-"(., 1:nrow(.)) %>% 
      mutate(Celltype_ZT = factor(Celltype_ZT, levels = levels_)) -> df_
    
    df_ %>% 
      group_by(Gene) %>%
      group_map(function(x,y){
        df_ = x
        df_ %>% mutate(Scaled_mean_expression = scales::rescale(x = Mean_expression, to = c(0,1))) -> df_
      }, .keep = T) %>% 
      do.call(rbind, .) %>% 
      mutate(Gene = factor(Gene, levels = rownames(motif_norm))) -> df_
    
    df_ %>% 
      ggplot(aes(x = Celltype_ZT, y = Gene)) + 
      geom_point(aes(size = non_zero_proportion, color = Scaled_mean_expression)) + 
#      scale_size(range = c(0, 6), limits = c(min(df_$non_zero_proportion), 100)) + 
      scale_size(range = c(0, 6), limits = c(0, 100)) + 
      scale_color_gradient(low = "grey", high = "purple") +
      scale_x_discrete(guide = guide_axis(angle = 90)) + 
      theme_classic() -> p_dot_scaled
    
    p_dot_scaled + ggtitle(sprintf("%s", seed_)) -> p_dot_scaled
    
  }) -> list_tmp_1
list_p_dot_motif_scaled = list_tmp_1 ; rm(list_tmp_1)
####

list(
  RNA = list_p_dot_RNA_scaled$seed_90$data, 
  ATAC_activity = list_p_dot_gene_activity_scaled$seed_90$data,
  Motif_score = list_p_dot_motif_scaled$seed_90$data
) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    assay_ = y
    df_ = x
    df_ %>% dplyr::mutate(Gene = stringr::str_to_title(Gene)) -> df_
    gene_level_ = c("Rorc", "Nr1d1", "Nfil3", "Dbp", "Npas2", "Clock", "Bhlhe41", "Bhlhe40", "Arntl")
    df_ %>% dplyr::filter(Gene %in% gene_level_) %>% 
      dplyr::mutate(Gene = factor(Gene, levels = gene_level_)) -> df_
    head(df_)
    color_ = list(RNA = "blue", ATAC_activity = "red", Motif_score = "purple")
    df_ %>% 
      ggplot(aes(x = Celltype_ZT, y = Gene)) + 
      geom_point(aes(size = non_zero_proportion, color = Scaled_mean_expression)) + 
      #      scale_size(range = c(0, 6), limits = c(min(df_$non_zero_proportion), 100)) + 
      scale_size(range = c(0, 6), limits = c(0, 100)) + 
      scale_color_gradient(low = "grey", high = color_[[assay_]]) +
      scale_x_discrete(guide = guide_axis(angle = 90)) + 
      theme_classic() -> p_dot_scaled
    p_dot_scaled + 
      xlab(NULL) + ylab(NULL) + 
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }) -> p_list
patchwork::wrap_plots(p_list, nrow = 3) #Fig_2A ----

# 4. Draw violin plot for celltype specific marker genes (RNA expression and ATAC activity) ----

library(tidyverse)
library(Seurat)
library(Signac)
library(patchwork)

# Read data
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hep_cells
readRDS("~/Dropbox/singulomics/github_rda/Endothelial_cells_cellnames.rds") -> endo_cells
readRDS("~/Dropbox/singulomics/github_rda/Fibroblasts_cellnames.rds") -> fibro_cells
readRDS("~/Dropbox/singulomics/github_rda/Kupffer_cells_cellnames.rds") -> kupffer_cells
readRDS("~/Dropbox/singulomics/github_rda/T_cells_cellnames.rds") -> t_cells
readRDS("~/Dropbox/singulomics/github_rda/B_cells_cellnames.rds") -> b_cells
readRDS("~/Dropbox/singulomics/github_rda/Cholangiocytes_cellnames.rds") -> cho_cells

sc = sc[,c(hep_cells, endo_cells, fibro_cells, kupffer_cells, t_cells, b_cells, cho_cells)]

source("~/Dropbox/singulomics/github/Stacked_Vlnplot.R")
c("Ebf1", "Slc5a1", "Stab2", "Dcn", "Egfr", "Clec4f", "Skap1") -> se_genes

StackedVlnPlot(obj = sc, features = se_genes, assay_ = "SCT") -> p_1 #Fig_1C ----
StackedVlnPlot(obj = sc, features = se_genes, assay_ = "gene_activity") -> p_2 #Fig_1C ----
####

# 5. Draw time-serise heatmap for marker genes in pseudo-bulk data (RNA expression and ) ----
library(tidyverse)
library(Seurat)
library(Signac)
#Read data
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc

#Generate normalized RNA expression matrix 
sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*median(sc@meta.data$nCount_RNA)})) -> RNA_norm
####

#Generate normalized ATAC acitivty matrix 
sc@assays[["ATAC"]]@fragments[[1]]@path <- "~/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
DefaultAssay(sc) <- "ATAC"
gene.activities <- GeneActivity(sc, extend.upstream = 2000, biotypes = NULL)

gene.activities %>% colSums() %>% median() -> gene_activity_median

gene.activities %>% 
  as.data.frame() %>%
  mutate(across(everything(), function(x){(x/sum(x))*gene_activity_median})) %>% 
  as.matrix() -> gene_activity_norm
gene_activity_norm %>% as.data.frame() -> gene_activity_norm
rm(gene.activities)
####

####
sc_meta = sc@meta.data
rm(sc)
gc()
####

#Extract rhythmic genes from Hughes liver data
read.csv("~/Dropbox/singulomics/github_rda/Hughes_Liver_circadian_genes.csv", header = T, stringsAsFactors = F) %>% 
  dplyr::filter(BH.Q < 0.01) %>% 
  dplyr::arrange(LAG) %>% 
  .$Gene.Symbol %>%
  str_split(" /// ") %>% 
  unlist() %>% 
  as.data.frame() %>% 
  "colnames<-"(., "Gene") %>% 
  dplyr::mutate(idx = 1:nrow(.)) %>% 
  dplyr::filter(Gene %in% rownames(RNA_norm)) %>% 
  dplyr::arrange(idx) %>% 
  .$Gene %>% unique() -> circadian_genes
####

Core_clock_genes = c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Npas2", "Dbp", "Nfil3", "Nr1d1", "Rorc", "Cry1", "Ciart", "Per1", "Per2")

#Extract KO cells
list_tmp = list()
sc_meta %>% dplyr::filter(grepl("KO", group)) %>% 
  {
    df_KO = .
    c("KO1", "KO2") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        group_ = x
        df_KO %>% filter(group == x) -> df_KO
        rownames(df_KO) -> cells_
        c("RNA", "ATAC") %>% 
          "names<-"(.,.) %>% 
          map(function(x){
            assay_ = x
            if (x == "RNA"){
              RNA_norm[circadian_genes, cells_] -> df_norm
            }
            if (x == "ATAC"){
                gene_activity_norm[circadian_genes, cells_] -> df_norm
            }
            list_tmp[[group_]][[assay_]] <<- df_norm
            df_norm %>% rowMeans() %>% as.data.frame() %>% "colnames<-"(., group_) -> df_norm
            df_norm
          })
      })
  } -> KO_list
####

#Extract WT samples
c("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds", 
  "~/Dropbox/singulomics/github_rda/Endothelial_cells_cellnames.rds",
  "~/Dropbox/singulomics/github_rda/Fibroblasts_cellnames.rds",
  "~/Dropbox/singulomics/github_rda/Kupffer_cells_cellnames.rds", 
  "~/Dropbox/singulomics/github_rda/Cholangiocytes_cellnames.rds",
  "~/Dropbox/singulomics/github_rda/B_cells_cellnames.rds", 
  "~/Dropbox/singulomics/github_rda/T_cells_cellnames.rds") %>% 
  map(function(x){
    readRDS(x) -> cells_
  }) %>% purrr::reduce(., c) -> WT_cells

sprintf("ZT%02d", seq(2,22,4)) %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    ZT_ = x
    sc_meta[WT_cells, ] %>% dplyr::filter(ZT == ZT_) -> df_meta
    rownames(df_meta) -> cells_
    c("RNA", "ATAC") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        assay_ = x
        if (x == "RNA"){
          RNA_norm[circadian_genes, cells_] -> df_norm
        }
        if (x == "ATAC"){
          gene_activity_norm[circadian_genes, cells_] -> df_norm
        }
        list_tmp[[ZT_]][[assay_]] <<- df_norm
        df_norm %>% rowMeans() %>% as.data.frame() %>% "colnames<-"(., ZT_) -> df_norm
        df_norm
      })
  }) -> WT_list

#Rhythmicity detection (JTKcycle and HR)
source("~/Dropbox/singulomics/github/Calculate_HMP.R")

WT_list %>% 
  map(function(list_){
    list_$RNA
  }) %>% do.call(cbind, .) %>% 
#  head() %>% 
  {
    df_ = .
    df_ %>% "colnames<-"(., sprintf("%s_REP1", colnames(.))) -> df_
    Gene_ = rownames(df_)
    df_ %>% as_tibble() -> df_
    timepoints_ = seq(2,22,4)
    cyclic_HMP(exp_matrix = df_, gene = Gene_, timepoints = timepoints_) -> res_
    res_ %>% dplyr::select(Gene, HR_phi) %>% dplyr::mutate(phase = (HR_phi/(2*pi))*24) -> res_
    res_
  } -> phase_df

circadian_genes %>% 
  as.data.frame() %>% 
  "colnames<-"(., "Gene") %>% 
  left_join(x = ., y = phase_df, by = "Gene") %>% 
  dplyr::arrange(phase) -> phase_df

phase_df$Gene -> circadian_genes
circadian_genes %>% rev() -> circadian_genes
####

#Generate heatmap matrix
c("RNA", "ATAC") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    type_ = x
    WT_list %>% 
      map(function(x){
        x[[type_]] -> df_
      }) %>% purrr::reduce(., cbind) -> df_
    apply(df_, 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
      t() %>% 
      as.data.frame() -> df_
  }) -> heatmap_mat

heatmap_mat %>% 
  map(function(df_){
    df_[phase_df$Gene, ] -> df_
  }) -> heatmap_mat

circadian_genes %>% ifelse(. %in% Core_clock_genes, ., "") -> labels_
ggplotify::as.ggplot(
  pheatmap::pheatmap(mat = heatmap_mat$RNA, cluster_rows = F, cluster_cols = F, labels_row = labels_, , color = colorRampPalette(c("#324084", "#f0e939"))(50))
) -> p_RNA

ggplotify::as.ggplot(
  pheatmap::pheatmap(mat = heatmap_mat$ATAC, cluster_rows = F, cluster_cols = F, labels_row = labels_, , color = colorRampPalette(c("#324084", "#f0e939"))(50))
) -> p_ATAC

patchwork::wrap_plots(p_RNA, p_ATAC, ncol = 2) + patchwork::plot_annotation(title = "1997 circadian genes (Hughes et al. 2009)") #Fig_1D ----
####

# 7. Generate celltype specifc rhythmic genes table (RNA expression and ATAC activity) ----
c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  map(function(celltype_){
    c("RNA", "gene_activity") %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(assay_){
        list.files("~/Dropbox/singulomics/github_rda/output/Celltype_specific/Cauchy_replicate_10", full.names = T) %>% 
          {.[grepl(celltype_, .)]} %>% {.[grepl(assay_, .)]} %>% gtools::mixedsort() -> files_
        files_ %>% 
#          .[1:2] %>% 
          map(function(file_){
            seed_ = gsub("/.+/.+_(seed_\\d+?)_.+", "\\1", file_)
            read.csv(file_, header = T, stringsAsFactors = F) -> df_
            df_ %>% dplyr::select(matches("Gene|JTK|HR|MetaCycle_meta2d_Base|MetaCycle_meta2d_AMP|MetaCycle_meta2d_rAMP")) -> df_
            df_[is.na(df_)] <- 1
            readRDS("~/Dropbox/singulomics/github_rda/output/Celltype_specific/intersect_genes_high_exp.rds") -> highly_expressed_genes
            df_ %>% column_to_rownames("Gene") %>% .[highly_expressed_genes,] -> df_
            colnames(df_) = sprintf("%s_%s", seed_, colnames(df_))
            df_ %>% rownames_to_column("Gene") -> df_
          }) %>% purrr::reduce(., left_join, by = "Gene") -> df_
        source("~/Dropbox/singulomics/github/Calculate_HMP.R")
        recal_cauchy_p_and_hmp(df_) -> df_
        df_ %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% arrange(cauchy_p) -> df_ 
        df_ %>% .[,1:3] -> df_pval
        
        radian_to_phase = function(radian){
          phase = (radian/(2*pi))*24
          return(phase)
        }
        df_ %>% dplyr::select(matches("HR_phi")) -> df_phi
        df_phi %>% dplyr::mutate(across(everything(), function(x){radian_to_phase(radian = x)})) -> df_phi
        phase = rowMeans(df_phi)
        df_pval %>% dplyr::mutate(Phase = phase) -> df_pval
        
        df_ %>% dplyr::select(matches("meta2d_AMP")) -> df_AMP
        Mean_Amp = apply(df_AMP, 1, function(x){mean(x)})
        SD_Amp = apply(df_AMP, 1, function(x){sd(x)})
        df_pval %>% dplyr::mutate(`Amplitude(Mean)` = Mean_Amp, `Amplitude(SD)` = SD_Amp) -> df_pval
        
        df_ %>% dplyr::select(matches("meta2d_Base")) -> df_Base
        Mean_Base = apply(df_Base, 1, function(x){mean(x)})
        SD_Base = apply(df_Base, 1, function(x){sd(x)})
        df_pval %>% dplyr::mutate(`Base(Mean)` = Mean_Base, `Base(SD)` = SD_Base) -> df_pval
        
        df_ %>% dplyr::select(matches("meta2d_rAMP")) -> df_rAMP
        Mean_rAMP = apply(df_rAMP, 1, function(x){mean(x)})
        SD_rAMP = apply(df_rAMP, 1, function(x){sd(x)})
        df_pval %>% dplyr::mutate(`rAMP(Mean)` = Mean_rAMP, `rAMP(SD)` = SD_rAMP) -> df_pval
        
        list.files("~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_data_downsampled/replicate_10", full.names = T) %>% 
          {.[grepl(assay_, .)]} %>% 
          {.[grepl(celltype_, .)]} %>% 
          #          .[1:2] %>% 
          map(function(file_){
            read.csv(file_, header = T, stringsAsFactors = F) -> df_2
            df_2 %>% pivot_longer(-Gene, names_to = "ZT", values_to = "value") -> df_2
            df_2$ZT = gsub("(ZT\\d+)_.+", "\\1", df_2$ZT)
            df_2 %>% group_by(Gene, ZT) %>% summarise(Mean = mean(value, na.rm = T)) %>% ungroup() -> df_2
          }) %>% do.call(rbind, .) %>% group_by(Gene, ZT) %>% summarise(Mean_expr = mean(Mean, na.rm = T)) %>% ungroup() -> df_2
        df_2 %>% pivot_wider(names_from = ZT, values_from = Mean_expr) -> df_2
        left_join(df_pval, df_2, by = "Gene") -> df_pval
      })
  }) -> df_list_supp_table_1 #Supp_table_1

# 7. Upset plots and venndiagram of cell celltype specfic genes (RNA expression and ATAC activity) ----
library(ggupset)
library(ggbreak)

c("RNA", "gene_activity") %>% 
  "names<-"(.,.) %>% 
  map(function(assay_){
    df_list_supp_table_1 %>% 
      map(function(list_){
        list_[[assay_]] %>% .$Gene
      })
  }) -> ggvenn_list

#RNA
ggvenn_list$RNA %>% 
  ggvenn::list_to_data_frame() %>% 
#  .[1:100, ] %>%  
  group_by(key) %>% 
  group_map(function(x,y){
    celltype_ = list(names(x)[x == T])
    tibble(Gene = y$key, Celltype = celltype_)
  }) %>% do.call(rbind, .) %>% 
  ggplot(aes(x=Celltype)) + 
  geom_bar() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(n_intersections = 30) + 
  theme_classic() +
#  scale_y_break(c(200,2000), scales = 0.5) +
  ggtitle("RNA") -> p1

ggvenn_list$RNA %>% 
  ggvenn::list_to_data_frame() %>% 
  #  .[1:100, ] %>%  
  group_by(key) %>% 
  group_map(function(x,y){
    celltype_ = list(names(x)[x == T])
    tibble(Gene = y$key, Celltype = celltype_)
  }) %>% do.call(rbind, .) %>% 
  ggplot(aes(x=Celltype)) + 
  geom_bar() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(n_intersections = 30) + 
  theme_classic() +
  scale_y_break(c(200,2000), scales = 0.5) +
  ggtitle("RNA") -> p1_1

ggvenn_list$gene_activity %>% 
  ggvenn::list_to_data_frame() %>% 
  #  .[1:100, ] %>%  
  group_by(key) %>% 
  group_map(function(x,y){
    celltype_ = list(names(x)[x == T])
    tibble(Gene = y$key, Celltype = celltype_)
  }) %>% do.call(rbind, .) %>% 
  ggplot(aes(x=Celltype)) + 
  geom_bar() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(n_intersections = 30) + 
  theme_classic() +
  ggtitle("Gene activity") -> p2

ggvenn_list$gene_activity %>% 
  ggvenn::list_to_data_frame() %>% 
  #  .[1:100, ] %>%  
  group_by(key) %>% 
  group_map(function(x,y){
    celltype_ = list(names(x)[x == T])
    tibble(Gene = y$key, Celltype = celltype_)
  }) %>% do.call(rbind, .) %>% 
  ggplot(aes(x=Celltype)) + 
  geom_bar() + 
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  scale_x_upset(n_intersections = 30) + 
  theme_classic() +
  scale_y_break(c(150,400), scales = 0.5) +
  ggtitle("Gene activity") -> p2_2
  
ggpubr::ggarrange(p1, p2, p1_1, p2_2, ncol = 2)
library(patchwork)
wrap_plots(p1, p1_1, ncol = 1, nrow = 2) #Fig_2B
wrap_plots(p2, p2_2, ncol = 1, nrow = 2) #Fig_2B

list(RNA = ggvenn_list$RNA$Hepatocytes, ATAC_activity = ggvenn_list$gene_activity$Hepatocytes) %>% 
  ggvenn::ggvenn() -> p_3 #Fig_2D

####

# 8. Celltype specific phase enrichment analysis ----

c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells") %>% 
  "names<-"(.,.) %>% 
  map(function(celltype_){
    df_list_supp_table_1[[celltype_]]$RNA %>% dplyr::select(Gene, Phase) -> df_
    output_file = sprintf("~/Dropbox/singulomics/github_rda/output/Celltype_specific/PSEA/input/%s_celltype_specific.txt", celltype_)
    write.table(df_, file = output_file, quote = F, row.names = F, col.names = F, sep = "\t")
    print(output_file)
  })

c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    celltype_ = x
    file_ = sprintf("~/Dropbox/singulomics/github_rda/output/Celltype_specific/PSEA/output/Celltype_specific/%s/results.txt", celltype_)
    read.table(file = file_, header = T, sep = "\t", stringsAsFactors = F) -> df_
    
    df_ %>% 
      filter(Kuiper.q.value..vs..uniform. < 0.05) %>% dplyr::arrange(Kuiper.q.value..vs..uniform.) -> df_
    df_ %>% dplyr::filter(grepl("GOBP_.+", Set.ID)) -> df_
    df_ %>% mutate(Celltype = celltype_, .before = 1) -> df_

  }) -> PSEA_res_list

PSEA_res_list %>% 
  map2(.x=.,.y=names(.),.f=function(df_, celltype){
    if (celltype == "Hepatocytes"){df_ %>% dplyr::filter(grepl("GOBP_PEPTIDE_METABOLIC_PROCESS|GOBP_REGULATION_OF_DEFENSE_RESPONSE|GOBP_LIPID_METABOLIC_PROCESS|GOBP_SMALL_MOLECULE_METABOLIC_PROCESS", Set.ID)) -> df_}
    if (celltype == "Endothelial_cells"){df_ %>% dplyr::filter(grepl("GOBP_CELL_CELL_ADHESION|GOBP_POSITIVE_REGULATION_OF_LOCOMOTION|CELL_CELL_SIGNALING|GOBP_ANGIOGENESIS", Set.ID)) -> df_}
    if (celltype == "Fibroblasts"){df_ %>% dplyr::filter(grepl("GOBP_CYTOSKELETON_ORGANIZATION|GOBP_SUPRAMOLECULAR_FIBER_ORGANIZATION|GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS|GOBP_REGULATION_OF_CELL_GROWTH", Set.ID)) -> df_}
    if (celltype == "Kupffer_cells"){df_ %>% dplyr::filter(grepl("GOBP_REGULATION_OF_DEFENSE_RESPONSE|GOBP_ENDOCYTOSIS|GOBP_LEUKOCYTE_DIFFERENTIATION|GOBP_ACTIVATION_OF_IMMUNE_RESPONSE", Set.ID)) -> df_}
    return(df_)
  }) -> PSEA_res_list_selected

generate_random_string <- function(length) {
  # Define the characters to sample from
  characters <- c(letters, LETTERS, 0:9)  # Lowercase, uppercase, and digits
  
  # Generate the random string
  random_string <- paste(sample(characters, length, replace = TRUE), collapse = "")
  
  return(random_string)
}

c(3,3,3,3,7) %>% 
  map(function(x){
    i = x
    1:i %>% 
      map(function(x){
        generate_random_string(20)
      }) %>% unlist() -> random_string
  }) -> random_string

i = 0
PSEA_res_list_selected %>% 
  map2(.x=.,.y=names(.),.f=function(df_, celltype){
    df_ = df_[,c(1,2,7,9)]
    colnames(df_) = c("Celltype", "Set.ID", "qvalue", "phase")
    
    i <<- i + 1
    data.frame(Celltype = "unknown", Set.ID = random_string[[i]], qvalue = 0, phase = 1) -> df_1
    rbind(df_, df_1)
  }) %>% do.call(rbind, .) %>% 
  {
    df_ = .
    data.frame(Celltype = "unknown", Set.ID = random_string[[5]], qvalue = 0, phase = 1) -> df_1
    rbind(df_, df_1)
  } %>% 
  dplyr::mutate(phase = round(phase)) %>% 
  mutate(Set.ID = factor(Set.ID, levels = rev(unique(Set.ID)))) %>% 
  ggplot(aes(x = phase, y = Set.ID, fill = qvalue)) + 
  geom_tile() + 
  coord_polar(theta = "x") + 
  scale_x_continuous(limits = c(0, 24), breaks = c(2,6,10,14,18,22)) #Fig_2C ----

####

# 9. Hepatocytes rhythmic ATAC peaks TF motif ----

library(tidyverse)
library(Seurat)
library(Signac)

# Load in the integrated and QC'ed data
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells

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

sc$nCount_ATAC %>% median() -> median_
sc@meta.data$cluster %>% unique() %>% 
  as.character() %>% 
  as.numeric() %>% 
  sort() %>% 
  "names<-"(.,sprintf("cluster_%s", .)) %>% 
  map(function(cluster_){
    rep_ = cluster_+1
    sc@meta.data %>% dplyr::filter(cluster == cluster_) -> df_meta
    seq(2,22,4) %>% 
      "names<-"(., sprintf("ZT%s", .)) %>% 
      map(function(ZT_){
        df_meta %>% dplyr::filter(ZT == ZT_) %>% rownames() -> cells_
        #        head(cells_)
        sc@assays$ATAC@counts[,cells_] -> df_
        df_ %>% as.data.frame() %>% 
          dplyr::mutate(across(everything(), function(x){(x/sum(x))*median_})) -> df_
        
        rowMeans(df_) %>% as.data.frame() %>% "colnames<-"(., sprintf("ZT%s_REP%s", ZT_, rep_)) -> df_
      }) %>% do.call(cbind, .) -> df_
    #    head(df_)
  }) %>% do.call(cbind, .) -> ATAC_df
gsub(".+\\.(.+)", "\\1", colnames(ATAC_df)) -> colnames(ATAC_df)

source("~/Dropbox/singulomics/github/Calculate_HMP.R")
ATAC_df %>% 
  {
    df_ = .
    gene_ = rownames(df_)
    timepoints_ = gsub("ZT(\\d+).+", "\\1", colnames(df_)) %>% as.integer()
    mat_ = df_ %>% as_tibble()
    cyclic_HMP(exp_matrix = mat_, gene = gene_, timepoints = timepoints_) -> res_
    res_
  } -> ATAC_df_res

ATAC_df_res %>% dplyr::select(matches("Gene|JTK|HR")) %>% recal_cauchy_p_and_hmp() -> ATAC_df_res
ATAC_df_res %>% dplyr::mutate(phase = (HR_phi/(2*pi))*24) -> ATAC_df_res
ATAC_df_res %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% dim()
save(ATAC_df, ATAC_df_res, file = "~/Downloads/roseplot_data.rda")

ATAC_df_res %>% dplyr::filter(cauchy_BH.Q < 0.05) %>% 
  dplyr::mutate(group = case_when(
    (phase > 0)&(phase <= 3) ~ "ZT0-3", 
    (phase > 3)&(phase <= 6) ~ "ZT3-6",
    (phase > 6)&(phase <= 9) ~ "ZT6-9",
    (phase > 9)&(phase <= 12) ~ "ZT9-12",
    (phase > 12)&(phase <= 15) ~ "ZT12-15",
    (phase > 15)&(phase <= 18) ~ "ZT15-18",
    (phase > 18)&(phase <= 21) ~ "ZT18-21",
    (phase > 21)&(phase <= 24) ~ "ZT21-24"
  )) -> ATAC_df_res_sig

DefaultAssay(sc) = "ATAC"
ATAC_df_res_sig$group %>% 
  unique() %>% 
  gtools::mixedsort() %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  map(function(group_){
    ATAC_df_res_sig %>% dplyr::filter(group == group_) -> df_
    #    nrow(df_)
    enriched.motifs = FindMotifs(object = sc, features = df_$Gene, background = ATAC_df_res$Gene) -> df_
  }) -> enriched.TF.motifs

enriched.TF.motifs %>% 
  map(function(group_){
    group_ %>% dplyr::mutate(motif.name = toupper(motif.name)) -> enriched.TF.motifs_
    genes_ = sc@assays$RNA@counts %>% rownames() %>% toupper()
    enriched.TF.motifs_ %>% dplyr::filter(motif.name %in% genes_) -> enriched.TF.motifs_
  }) -> enriched.TF.motifs.mm10

  ATAC_df_res_sig$group %>% unique() %>% 
  gtools::mixedsort() %>% 
  map(function(group_){
    ATAC_df_res_sig %>% dplyr::filter(group == group_) -> df_
    data.frame(group = group_, freq = nrow(df_))
  }) %>% do.call(rbind, .) %>% 
  {
    df_ = .
    df_$group -> levels_
    df_$group %>% factor(levels = levels_) -> df_$group
    df_
  } %>% 
  ggplot(aes(x = group, y = freq, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(start = 0) #Fig_2F

####