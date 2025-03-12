library(tidyverse)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(monocle3)

# 1. Cell partitioning using monocle3 ----
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc=sc[,!grepl('KO',sc$group)] #exclude KO cells
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
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

colnames(sc_1)[partitions(sc_1) == 1] -> hepatocyte_ids
colData(sc_1)[hepatocyte_ids, ] %>% as.data.frame() %>%
  filter(celltype == "Hepatocytes") %>% 
  rownames() -> hepatocyte_ids
sc_1[,hepatocyte_ids] -> sc_hepatocyte

sc_hepatocyte = cluster_cells(sc_hepatocyte, resolution = 0.35e-4)
plot_cells(sc_hepatocyte, color_cells_by = "cluster", show_trajectory_graph = FALSE)
plot_cells(sc_hepatocyte, color_cells_by = "partition", show_trajectory_graph = FALSE)

sc_hepatocyte <- learn_graph(sc_hepatocyte, use_partition = TRUE, verbose = FALSE, close_loop = FALSE)

plot_cells(sc_hepatocyte,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           show_trajectory_graph = F) +
  theme(legend.position = "top")

sc_hepatocyte <- order_cells(sc_hepatocyte, 
                             root_cells = colnames(sc_hepatocyte[,clusters(sc_hepatocyte) == 1]))

plot_cells(sc_hepatocyte,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "red") + 
  #  scale_color_gradientn(colors = c("red", "blue")) +
  theme(legend.position = "top") -> p1

plot_cells(sc_hepatocyte,
           color_cells_by = "ZT",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           show_trajectory_graph = F) +
  theme(legend.position = "top") -> p2

ggpubr::ggarrange(p1, p2, ncol = 2, legend = "bottom")

sc_hep_graph_test_results <- graph_test(sc_hepatocyte,
                                        neighbor_graph = "principal_graph",
                                        cores = 6)

sc_hep_graph_test_results %>% 
  filter(q_value < 0.05) %>% 
  dplyr::arrange(desc(morans_I), q_value) %>% 
  rownames() -> deg_ids

rowData(sc_hepatocyte)$gene_name <- rownames(sc_hepatocyte)
rowData(sc_hepatocyte)$gene_short_name <- rowData(sc_hepatocyte)$gene_name

plot_cells(sc_hepatocyte,
           #           genes=se_genes,
           genes=head(deg_ids, 12),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

map(deg_ids[1:12], function(x) FeaturePlot(sc, features = x, reduction = "multicca.umap")) -> p_list
ggpubr::ggarrange(plotlist = p_list, ncol = 4, nrow = 3)

transient_genes = c("Cyp2f2", "Cyp2e1", "Glul")

sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
  as.matrix() -> RNA_norm

transient_genes %>% "names<-"(.,.) %>% 
  #  .[1] %>%
  map(function(x){
    RNA_norm[x, ] %>% as.data.frame() %>% rownames_to_column("cell_name") %>% "colnames<-"(.,c("cell_name", "exp")) -> gene_exp
    pseudotime(sc_hepatocyte) %>% as.data.frame() %>% rownames_to_column("cell_name") %>% "colnames<-"(.,c("cell_name", "pseudotime")) -> pseudotime
    left_join(pseudotime, gene_exp, by = "cell_name") %>% 
      arrange(pseudotime) %>% 
      group_by(pseudotime) %>%
      group_map(function(x,y){
        mean(x$exp) -> mean_
        y$pseudotime -> pseudotime
        data.frame(pseudotime = pseudotime, mean = mean_)
      }, .keep = T) %>% 
      do.call(rbind, .) -> df_
    df_ %>% 
      ggplot(aes(x = pseudotime, y = mean)) +
      geom_smooth() + 
      ylab("TP10K") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      ggtitle(x) +
      theme(plot.title = element_text(hjust = 0.5))
  }) -> p_list
ggpubr::ggarrange(plotlist = p_list, ncol = 3, nrow = 1)


pseudotime(sc_hepatocyte) %>% as.data.frame() %>% "colnames<-"(., "pseudotime") %>% 
  rownames_to_column("cell_name") %>% 
  {left_join(x = ., y = sc@meta.data %>% rownames_to_column("cell_name"), by = "cell_name")} %>% 
  dplyr::select(pseudotime, ZT) %>% 
  arrange(pseudotime) %>%
  distinct() %>% 
  ggplot(aes(x = pseudotime, group = ZT, color = ZT)) + 
  geom_density() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

pseudotime(sc_hepatocyte) %>% as.data.frame() %>% "colnames<-"(., "pseudotime") %>% 
  rownames_to_column("cell_name") %>% 
  {left_join(x = ., y = sc@meta.data %>% rownames_to_column("cell_name"), by = "cell_name")} %>% 
  dplyr::select(pseudotime, ZT) %>% 
  arrange(pseudotime) %>%
  distinct() %>% 
  ggplot(aes(x = ZT, y = pseudotime, fill = ZT)) +
  geom_violin() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggarrange(p1,p2,p3, ncol = 3, nrow = 1, legend = "bottom")

sc_hep_graph_test_results %>% 
  filter(q_value < 0.05) %>% 
  dplyr::arrange(desc(morans_I), q_value) %>% 
  rownames() %>% .[1:35] -> deg_ids_tmp

all(deg_ids_tmp %in% rownames(sc))

pseudotime(sc_hepatocyte) %>% floor() %>% unique() %>% 
  "names<-"(., sprintf("T%s", .)) %>%
  sort() %>% 
  map2(.x = ., .y = names(.), function(x,y){
    pseudotime(sc_hepatocyte) %>% floor() %>% 
      {.[. == x]} %>% 
      names() -> cells_
    RNA_norm[deg_ids_tmp, cells_] %>% 
      #      .[1:5,] %>% 
      #      .[c("Cyp2f2", "Cyp2e1", "Glul"), ] %>%
      rowMeans() %>% 
      as.data.frame() %>% 
      "colnames<-"(., y)
  }) %>% do.call(cbind, .) %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() -> df_tmp

pheatmap::pheatmap(df_tmp, cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = T, cutree_rows = 2)

# 2. Use SlingShot to infer pseudo time in Hepatocytes ----
colnames(sc_hepatocyte) %>% 
  write_rds(., "~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds")

library(slingshot)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(tidyverse)

readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc

colnames(sc) %>% length()
length(hepatocytes_cells)

sc <- sc[,hepatocytes_cells]

DimPlot(sc, group.by = "group", label = TRUE, reduction = "multicca.umap", repel = TRUE)
DefaultAssay(sc)='SCT'
sce=as.SingleCellExperiment(sc, assay = "SCT")

reducedDims(sce)

sce <- slingshot(sce,reducedDim = "MULTICCA.UMAP",
                 allow.breaks = FALSE)
summary(sce$slingPseudotime_1)

pseudotime=sce$slingPseudotime_1
pseudotime=pseudotime/max(pseudotime)
#pseudotime = 1 - pseudotime
pseudotime=round(pseudotime,4)
sc=AddMetaData(sc,metadata = pseudotime, col.name = 'pseudotime')

library(tradeSeq)
library(BiocParallel)
library(batchtools)

# fit negative binomial GAM
#register(BatchtoolsParam(workers=4))
#BPPARAM <- BiocParallel::bpparam()
bpparam()
multicoreParam <- MulticoreParam(workers = 4)
multicoreParam
sce <- fitGAM(sce, parallel = T, BPPARAM = multicoreParam)
#load("~/Downloads/fitGAM_results.rda")

# test for dynamic expression
ATres <- associationTest(sce)
ATres %>% 
  dplyr::arrange(pvalue, desc(waldStat)) %>% 
  head(35) %>% rownames() -> topgenes
ATres %>% 
  dplyr::arrange(pvalue, desc(waldStat)) %>% 
  head(1000) %>% rownames() -> topgenes_1000
FeaturePlot(sc, features = c("Cryl1"), reduction = "multicca.umap")

pst.ord <- order(sce$slingPseudotime_1, na.last = NA, decreasing = T)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatdata_genomewide <- assays(sce)$counts[, pst.ord]
#heatclus <- sce$GMM[pst.ord]

class(heatdata)
pheatmap::pheatmap(log1p(heatdata), cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F)
pheatmap::pheatmap(log1p(heatdata_genomewide), cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F)

rm(sce, pseudotime)
hist(sc$pseudotime)

p1=FeaturePlot(sc, features = 'pseudotime', reduction = "multicca.umap")
p2=DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "multicca.umap", repel = TRUE)
p3=VlnPlot(sc, features='pseudotime', group.by='ZT', pt.size = 0)
ggpubr::ggarrange(p1,p2,p3, ncol = 3, nrow = 1, legend = "top")

transient_genes = c("Cyp2f2", "Cyp2e1", "Glul")

sc@assays$RNA@counts %>% 
  as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*10^6})) %>% 
  as.matrix() -> RNA_norm

transient_genes %>% "names<-"(.,.) %>% 
  #  .[1] %>%
  map(function(x){
    RNA_norm[x, ] %>% as.data.frame() %>% rownames_to_column("cell_name") %>% "colnames<-"(.,c("cell_name", "exp")) -> gene_exp
    gene_exp %>% mutate(pseudotime = sc@meta.data$pseudotime) -> gene_exp
    gene_exp %>% arrange(pseudotime) %>% 
      group_by(pseudotime) %>% 
      summarise(mean = mean(exp)) -> df_
    
    df_ %>% 
      ggplot(aes(x = pseudotime, y = mean)) +
      geom_smooth() + 
      ylab("TP10K") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      ggtitle(x) +
      theme(plot.title = element_text(hjust = 0.5))
  }) -> p_list
ggpubr::ggarrange(plotlist = p_list, ncol = 3, nrow = 1)

save.image(file = "~/Dropbox/singulomics/github_rda/trajectory_analysis.rda")