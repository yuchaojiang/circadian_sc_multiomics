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