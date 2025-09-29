setwd("~/Dropbox/singulomics")

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
library(scDblFinder)
library(mbkmeans)
library(ggpubr)

load('rda/2_QC.rda')

# 1. Labels transfers (7 cell types)----
############################################
# Liver cell atlas: nuclei nucseq
############################################
load('./rda/reference_sc_data/LCA/lca_nuc.rda')
lca_nuc[["percent.mt"]] <- PercentageFeatureSet(lca_nuc, pattern = "^mt-")
VlnPlot(lca_nuc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0)
lca_nuc <- subset(lca_nuc, subset = nFeature_RNA < 5000 &
                    nCount_RNA<20000 & percent.mt < 1)
VlnPlot(lca_nuc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,  pt.size = 0)

DefaultAssay(lca_nuc) <- "RNA"
lca_nuc <- SCTransform(lca_nuc, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DimPlot(lca_nuc, reduction = "umap.rna", group.by = 'annot', label = TRUE) -> p1
DimPlot(lca_nuc, reduction = "umap.rna", group.by = 'sample', label = TRUE) -> p2
ggpubr::ggarrange(p1, p2, ncol = 1, nrow = 2, legend = "top")


###Plot lca_nuc_subset integrated umap ----
lca_nuc_subset = lca_nuc[,lca_nuc$annot %in% c('Endothelial cells','Fibroblasts','Hepatocytes','Kupffer cells', "B cells", "Cholangiocytes", 
                                                "T cells")]
lca_nuc_subset@meta.data %>% dplyr::filter(sample != "ABU21") %>% rownames() -> se_
lca_nuc_subset = lca_nuc_subset[,se_]

lca_nuc_subset_splited = SplitObject(lca_nuc_subset, split.by = "sample")

for (i in 1:length(lca_nuc_subset_splited)) {
  lca_nuc_subset_splited[[i]] <- SCTransform(lca_nuc_subset_splited[[i]], 
                                             vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"), assay="RNA", verbose = TRUE, return.only.var.genes = FALSE)
}

#Find common SCT assay genes among all groups
v.common <- rownames(lca_nuc_subset_splited[[1]]$SCT@data)
for (  i in 1:length(lca_nuc_subset_splited)){
  if (i == 1){
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(v.common))) 
  }else{
    v.common <- intersect(v.common, rownames(lca_nuc_subset_splited[[i]]$SCT@data))
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(rownames(lca_nuc_subset_splited[[i]]$SCT@data))))
  }
} 
length( v.common );
cat( "length(v.common)=", length(v.common), "\n" )

v.inte.features <- vector()
if ( 1 ){                                     # find variable features/genes common to all
  v.var.common <- lca_nuc_subset_splited[[1]]@assays$SCT@var.features
  print(sprintf("condition %s number of SCT var.features: %s", 1, length(v.var.common)))
  for (i in 2:length(lca_nuc_subset_splited)){
    v.var.common <- c(v.var.common, lca_nuc_subset_splited[[i]]@assays$SCT@var.features)
    print(sprintf("condition %s number of SCT var.features: %s", i, length(lca_nuc_subset_splited[[i]]@assays$SCT@var.features)))
  }
  length(unique(v.var.common))
  v.var.common <- names(which(table(v.var.common) == length(lca_nuc_subset_splited)))
  
  cat("length(v.var.common)=", length(v.var.common), "\n")
  length(v.var.common)
  v.inte.features <- v.var.common
} else {
  v.inte.features <- SelectIntegrationFeatures(object.list = lca_nuc_subset_splited, nfeatures = 3000)
}

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
#k.weight.use = 50
k.weight.use = 80
lca_nuc_subset_splited.anchors <- FindIntegrationAnchors(object.list = lca_nuc_subset_splited, dims = 1:nDims, 
                                                         assay=rep( "SCT", length(lca_nuc_subset_splited)), 
                                                         anchor.features=v.inte.features, k.anchor=k.anchor.use, 
                                                         k.filter=k.filter.use, k.score=k.score.use, verbose=TRUE )
lca_nuc_integrated = IntegrateData(anchorset = lca_nuc_subset_splited.anchors, dims = 1:nDims, features.to.integrate=v.common, 
                                   k.weight = k.weight.use, verbose = T, new.assay.name = "RNA_integrated")

View(lca_nuc_integrated)


VariableFeatures(lca_nuc_integrated, assay="RNA_integrated") <- v.inte.features
lca_nuc_integrated@assays$RNA_integrated@key <- "rnaintegrated_"

DefaultAssay(lca_nuc_integrated) <- "RNA_integrated"
lca_nuc_integrated <- ScaleData(lca_nuc_integrated, verbose=TRUE, assay="RNA_integrated", features=rownames(lca_nuc_integrated$RNA_integrated@data))

lca_nuc_integrated <- RunPCA(lca_nuc_integrated, assay="RNA_integrated", npcs = nDims, verbose = FALSE, 
                             reduction.name = "rna_integrated_pca", 
                             reduction.key = "rnaintegratedpc_")
lca_nuc_integrated <- RunUMAP(lca_nuc_integrated, reduction = "rna_integrated_pca", dims = 1:nDims, 
                              reduction.name = "rna_integrated_umap", 
                              reduction.key = "rnaintegratedumap_", 
                              assay = "RNA_integrated")
lca_nuc_integrated <- FindNeighbors(lca_nuc_integrated, assay="RNA_integrated", reduction = "rna_integrated_pca", dims = 1:nDims, force.recalc=T, 
                                    graph.name = c("RNA_integrated_nn", "RNA_integrated_snn"))
lca_nuc_integrated <- FindClusters(lca_nuc_integrated, assay="RNA_integrated", resolution = 1.0, graph.name = "RNA_integrated_snn")

lca_nuc_integrated@meta.data %>% mutate(rna_integrated_cluster = seurat_clusters) -> lca_nuc_integrated@meta.data

DimPlot(lca_nuc_integrated, group.by = "rna_integrated_cluster", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(lca_nuc_integrated, group.by = "annot", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE) -> p1 # Supp_Fig_4A ----
DimPlot(lca_nuc_integrated, group.by = "sample", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE) -> p2 # Supp_Fig_4A ----
ggpubr::ggarrange(p1, p2, ncol = 2, legend = "top")
save(lca_nuc_integrated, file = "~/Dropbox/singulomics/github_rda/lca_nuc_integrated.rda")
####


#Transfer cell type labels
DefaultAssay(sc) <- "SCT"
DefaultAssay(lca_nuc) <- "SCT"

sc.anchors <- FindTransferAnchors(reference = lca_nuc, query = sc, normalization.method = "SCT", reference.reduction = "pca", dims = 1:30)
predictions <- TransferData(anchorset = sc.anchors, refdata = lca_nuc$annot, dims = 1:30)
sc <- AddMetaData(sc, metadata = predictions)

hist(sc$prediction.score.max,20)
sc$prediction.score.max %>% 
  as.data.frame() %>% 
  "colnames<-"(., "Prediction_score") %>% 
  ggplot(aes(x = Prediction_score)) + 
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") -> p_hist_1

sc=sc[,sc$prediction.score.max > 0.5]
sc$prediction.score.max %>% 
  as.data.frame() %>% 
  "colnames<-"(., "Prediction_score") %>% 
  ggplot(aes(x = Prediction_score)) + 
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") + 
  theme_classic() -> p_hist_2

c("B cells", "Cholangiocytes", "Endothelial cells", "Fibroblasts", "Hepatocytes", "Kupffer cells", "T cells") %>% 
  "names<-"(.,.) %>% 
  purrr::map(function(x){
    celltype_ = x
    sc@meta.data %>% dplyr::filter(predicted.id == celltype_) %>% .$prediction.score.max %>% 
      as.data.frame() %>% "colnames<-"(., "Prediction_score") %>% 
      "colnames<-"(., "Prediction_score") %>% 
      ggplot(aes(x = Prediction_score)) + 
      geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") + 
      theme_classic() -> p
    p + ggtitle(celltype_) -> p
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 7)

DimPlot(sc, reduction = "umap.rna", group.by='predicted.id', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA")
DimPlot(sc, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA") 

rm(lca_nuc, predictions, sc.anchors)

table(sc$predicted.id, sc$group)

############################################
# QC
############################################
# # Same predictions from both references
# sc=sc[,sc$predicted.id==sc$predicted.id.lca]

# Only keep major cell types
sc=sc[,sc$predicted.id %in% c('Endothelial cells','Fibroblasts','Hepatocytes','Kupffer cells', "B cells", "Cholangiocytes", 
                              "T cells")]

dim(sc@assays$RNA); length(gene.ref)
dim(sc@assays$ATAC); length(peak.ref)

table(sc$predicted.id, sc$group)

############################################
# Renormalize
############################################

# 2. RNA analysis (Normailization) ----
DefaultAssay(sc) <- "RNA"
sc <- SCTransform(sc, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# 3. ATAC analysis (Normalization) ----
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sc) <- "ATAC"
sc <- RunTFIDF(sc)
sc <- FindTopFeatures(sc, min.cutoff = 'q0')
sc <- RunSVD(sc)
sc <- RunUMAP(sc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sc <- FindMultiModalNeighbors(sc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sc <- RunUMAP(sc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sc <- FindClusters(sc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

sc$celltype=sc$predicted.id

# Some initial visualization
p1 <- DimPlot(sc, reduction = "umap.rna", group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA") 
p2 <- DimPlot(sc, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA") 
p3 <- DimPlot(sc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc ATAC")
p4 <- DimPlot(sc, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc ATAC")
p5 <- DimPlot(sc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p6 <- DimPlot(sc, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

ggarrange(plotlist = list(p1, p3), nrow = 2, ncol = 1, legend = "bottom", common.legend = T)
rm(p1,p2,p3,p4,p5,p6)

# 3 .scRNA seq integration -----
#load('rda/3_transferred.rda')
sc@meta.data$ZT %>% table()

sc_splited = SplitObject(sc, split.by = "group")
sc@assays$RNA@data %>% rownames() %>% length()

for (i in 1:length(sc_splited)) {
  sc_splited[[i]] <- SCTransform(sc_splited[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC"), assay="RNA", verbose = TRUE, return.only.var.genes = FALSE)
}

#Find common SCT assay genes among all groups
v.common <- rownames(sc_splited[[1]]$SCT@data)
for (  i in 1:length(sc_splited)){
  if (i == 1){
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(v.common))) 
  }else{
    v.common <- intersect(v.common, rownames(sc_splited[[i]]$SCT@data))
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(rownames(sc_splited[[i]]$SCT@data))))
  }
} 
length( v.common );
cat( "length(v.common)=", length(v.common), "\n" )

v.inte.features <- vector()
if ( 1 ){                                     # find variable features/genes common to all
  v.var.common <- sc_splited[[1]]@assays$SCT@var.features
  print(sprintf("condition %s number of SCT var.features: %s", 1, length(v.var.common)))
  for (i in 2:length(sc_splited)){
    v.var.common <- c(v.var.common, sc_splited[[i]]@assays$SCT@var.features)
    print(sprintf("condition %s number of SCT var.features: %s", i, length(sc_splited[[i]]@assays$SCT@var.features)))
  }
  length(unique(v.var.common))
  v.var.common <- names(which(table(v.var.common) == length(sc_splited)))
  
  cat("length(v.var.common)=", length(v.var.common), "\n")
  length(v.var.common)
  v.inte.features <- v.var.common
} else {
  v.inte.features <- SelectIntegrationFeatures(object.list = sc_splited, nfeatures = 3000)
}

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
k.weight.use = 100
sc_splited.anchors <- FindIntegrationAnchors(object.list = sc_splited, dims = 1:nDims, 
                                             assay=rep( "SCT", length(sc_splited)), 
                                             anchor.features=v.inte.features, k.anchor=k.anchor.use, 
                                             k.filter=k.filter.use, k.score=k.score.use, verbose=TRUE )

#Run the integration in HPC
IntegrateData_slurm = function(sc_splited.anchors, nDims, v.common, k.weight.use){
  sc_integrated <- IntegrateData(anchorset = sc_splited.anchors, dims = 1:nDims, features.to.integrate=v.common, 
                                 k.weight = k.weight.use, verbose = T, new.assay.name = "RNA_integrated")
  return(sc_integrated)
}

rm(sc_splited)

sjob = slurm_call(IntegrateData_slurm, jobname = "anchor_integration_RNA_T_cell_etc", 
                  list(sc_splited.anchors = sc_splited.anchors, 
                       nDims = nDims, 
                       v.common = v.common, 
                       k.weight.use = k.weight.use), submit = FALSE, libPaths = .libPaths(), 
                  slurm_options = list(time = "10:00:00", ntasks = "4", mem = "150G"))

# Get RNA sc_integrated from HPC's results
readRDS("~/Dropbox/singulomics/_rslurm_anchor_integration_RNA_T_cell_etc/results_0.RDS") -> sc_integrated

VariableFeatures(sc_integrated, assay="RNA_integrated") <- v.inte.features
sc_integrated@assays$RNA_integrated@key <- "rnaintegrated_"

DefaultAssay(sc_integrated ) <- "RNA_integrated"
sc_integrated <- ScaleData(sc_integrated, verbose=TRUE, assay="RNA_integrated", features=rownames(sc_integrated$RNA_integrated@data))

sc_integrated <- RunPCA(sc_integrated, assay="RNA_integrated", npcs = nDims, verbose = FALSE, 
                        reduction.name = "rna_integrated_pca", 
                        reduction.key = "rnaintegratedpc_")
sc_integrated <- RunUMAP(sc_integrated, reduction = "rna_integrated_pca", dims = 1:nDims, 
                         reduction.name = "rna_integrated_umap", 
                         reduction.key = "rnaintegratedumap_", 
                         assay = "RNA_integrated")
sc_integrated <- FindNeighbors(sc_integrated, assay="RNA_integrated", reduction = "rna_integrated_pca", dims = 1:nDims, force.recalc=T, 
                               graph.name = c("RNA_integrated_nn", "RNA_integrated_snn"))
sc_integrated <- FindClusters(sc_integrated, assay="RNA_integrated", resolution = 1.0, graph.name = "RNA_integrated_snn")

sc_integrated@meta.data %>% mutate(rna_integrated_cluster = seurat_clusters) -> sc_integrated@meta.data

DimPlot(sc_integrated, group.by = "rna_integrated_cluster", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "celltype", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "group", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)

table(sc_integrated$celltype, sc_integrated$rna_integrated_cluster)
#                     0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
#B cells              0    0    0    0    0    1    0    0    1    0    1    0    0    0    0   22  139    0
#Cholangiocytes       0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0    0
#Endothelial cells    0    0    0 2208    0    0    0    0    0    0    0    0    0    0  154   65    0    0
#Fibroblasts          0    0    0    0    0    0    0    0    0    0 1200    1    0    0    0   31    0    0
#Hepatocytes       3916 3824 2567    1 2158 1807 1791 1747 1491 1411   52  917   52  748  222   11   46    7
#Kupffer cells        0    0    0    0    0    0    0    0    0    1    2    0  808    0    2  161    1    1
#T cells              0    0    0    0    0    0    0    0    0    0    0    0    1    0    0   26    0  172
#
#.                    18   19   20   21   22   23
#B cells              0    0    0    0    0    0
#Cholangiocytes       0  110    1    0    0    0
#Endothelial cells  110    0    0    0   55   35
#Fibroblasts          0    0    0   62    0    0
#Hepatocytes          6    2   88    1    1    2
#Kupffer cells        0    0    0    0    0    0
#T cells              0    0    0    0    0    0

table(sc_integrated$group, sc_integrated$rna_integrated_cluster)

## 3.1 Recluster integrated scRNA assay ----
cluster_sele = c(0:13, 16:19) # exclude minor clusters (low cell numbers) and clusters mixed with different cell types
sc_integrated@meta.data %>% rownames_to_column(var = "cellnames") %>% 
  filter(seurat_clusters %in% cluster_sele) %>% 
  group_by(seurat_clusters) %>% 
  group_map(function(x,y){
    cbind(x,y) -> df_
    table(df_$celltype) %>% sort(decreasing = T) -> celltype_
    df_ %>% filter(celltype == names(celltype_)[1]) -> df_
    df_
  }) %>% do.call(rbind, .) %>% 
  .$cellnames -> cell_sele

sc_integrated <- subset(sc_integrated, cells = cell_sele) # select cells from selected clusters

sc_integrated@meta.data %>% 
  .$seurat_cluster %>% table()

VariableFeatures(sc_integrated, assay="RNA_integrated") <- v.inte.features

#Recluster RNA_integrated with cell filtering
DefaultAssay(sc_integrated ) <- "RNA_integrated"
sc_integrated <- ScaleData(sc_integrated, verbose=TRUE, assay="RNA_integrated", features=rownames(sc_integrated$RNA_integrated@data))

sc_integrated <- RunPCA(sc_integrated, assay="RNA_integrated", npcs = nDims, verbose = FALSE, 
                        reduction.name = "rna_integrated_pca", 
                        reduction.key = "rnaintegratedpc_")
sc_integrated <- RunUMAP(sc_integrated, reduction = "rna_integrated_pca", dims = 1:nDims, 
                         reduction.name = "rna_integrated_umap", 
                         reduction.key = "rnaintegratedumap_", 
                         assay = "RNA_integrated")
sc_integrated <- FindNeighbors(sc_integrated, assay="RNA_integrated", reduction = "rna_integrated_pca", dims = 1:nDims, force.recalc=T, 
                               graph.name = c("RNA_integrated_nn", "RNA_integrated_snn"))
sc_integrated <- FindClusters(sc_integrated, assay="RNA_integrated", resolution = 1.0, graph.name = "RNA_integrated_snn")

sc_integrated@meta.data %>% mutate(rna_integrated_cluster = seurat_clusters) -> sc_integrated@meta.data

DimPlot(sc_integrated, group.by = "rna_integrated_cluster", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "celltype", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "group", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "ZT", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)

saveRDS(object = sc_integrated, file = "~/Downloads/sc_integrated_rna.rds")
rm(sc_integrated)
###

# 4. scATAC seq Integration ----
readRDS(file = "~/Downloads/sc.rds") -> sc

DefaultAssay(sc) <- "ATAC"
sc_splited = SplitObject(sc, split.by = "group")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = sc_splited,
  anchor.features = rownames(sc),
  reduction = "rlsi",
  dims = 2:50
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = sc[["lsi"]],
  new.reduction.name = "atac_integrated_lsi",
  dims.to.integrate = 1:50 
)
#key = "atacintegratedlsi_"

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "atac_integrated_lsi", 
                      reduction.name = "atac_integrated_umap", 
                      reduction.key = "atacintegratedumap_", 
                      dims = 2:50)
integrated <- FindNeighbors(integrated, assay="ATAC", reduction = "atac_integrated_lsi", dims = 1:50, force.recalc=T, 
                            graph.name = c("atac_integrated_nn", "atac_integrated_snn"))
integrated <- FindClusters(integrated, assay="ATAC", resolution = 1.0, graph.name = "atac_integrated_snn")
integrated@meta.data %>% mutate(atac_integrated_cluster = seurat_clusters) -> integrated@meta.data

DimPlot(integrated, group.by = "atac_integrated_cluster", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(integrated, group.by = "celltype", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(integrated, group.by = "group", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(integrated, group.by = "ZT", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)

saveRDS(object = integrated, "~/Downloads/sc_integrated_atac.rds")

integrated_atac <- integrated
rm(integrated)
table(integrated_atac$celltype, integrated_atac$atac_integrated_cluster)
table(integrated_atac$group, integrated_atac$atac_integrated_cluster)

## 4.1 Recluster integrated scATAC assay ----
cluster_sele = c(0:12, 14, 15, 17)
integrated_atac@meta.data %>% rownames_to_column(var = "cellnames") %>%
  filter(seurat_clusters %in% cluster_sele) %>%
  group_by(seurat_clusters) %>%
  group_map(function(x,y){
    cbind(x,y) -> df_
    df_$celltype %>% table() %>% sort(decreasing = T) -> celltype_
    df_ %>% filter(celltype == names(celltype_)[1]) -> df_
    df_
  }) %>%
  do.call(rbind, .) %>% .$cellnames -> cell_sele
#25890 cells
integrated_atac <- subset(integrated_atac, cells = cell_sele)

integrated_atac <- RunUMAP(integrated_atac, reduction = "atac_integrated_lsi", dims = 2:50, reduction.name = "atac_integrated_umap", 
                           reduction.key = "atacintegratedumap_", assay = "ATAC")
integrated_atac <- FindNeighbors(integrated_atac, assay="ATAC", reduction = "atac_integrated_lsi", dims = 1:50, force.recalc=T, 
                                 graph.name = c("ATAC_integraed_nn", "ATAC_integrated_snn"))
integrated_atac <- FindClusters(integrated_atac, assay="ATAC", resolution = 1.0, graph.name = "ATAC_integrated_snn")

integrated_atac@meta.data %>% mutate(atac_integrated_cluster = seurat_clusters) -> integrated_atac@meta.data

DimPlot(integrated_atac, group.by = "group", label = T, reduction = "atac_integrated_umap")
DimPlot(integrated_atac, group.by = "ZT", label = T, reduction = "atac_integrated_umap")
DimPlot(integrated_atac, group.by = "atac_integrated_cluster", label = T, reduction = "atac_integrated_umap")
DimPlot(integrated_atac, group.by = "celltype", label = T, reduction = "atac_integrated_umap")

saveRDS(object = integrated_atac, "~/Downloads/sc_integrated_atac.rds")
readRDS("~/Downloads/sc_integrated_rna.rds") -> integrated_rna

source("~/Dropbox/singulomics/0_misc.R")

table(integrated_atac$celltype, integrated_atac$atac_integrated_cluster)
integrated_atac=keep_major_celltype(integrated_atac) # keep major cell types in each cluster
table(integrated_atac$celltype, integrated_atac$atac_integrated_cluster)

table(integrated_rna$celltype, integrated_rna$rna_integrated_cluster)
integrated_rna=keep_major_celltype(integrated_rna) # keep major cell types in each cluster
table(integrated_rna$celltype, integrated_rna$rna_integrated_cluster)

DimPlot(integrated_rna, group.by = "celltype", label = T, reduction = "rna_integrated_umap") + ggtitle("Integrated RNA")-> p1
DimPlot(integrated_atac, group.by = "celltype", label = T, reduction = "atac_integrated_umap") + ggtitle("Integrated ATAC") -> p2
ggarrange(plotlist = list(p1, p2), nrow = 2, ncol = 1, legend = "bottom", common.legend = T)

# 5. Intersected cells in integrated RNA and integrated ATAC ----
# Match the cells in the right orders
cells.to.keep=intersect(colnames(integrated_rna),colnames(integrated_atac))

sc=sc[,match(cells.to.keep, colnames(sc))] #subset sc

integrated_rna=integrated_rna[,match(cells.to.keep, colnames(integrated_rna))]
integrated_atac=integrated_atac[,match(cells.to.keep, colnames(integrated_atac))]

dim(sc); dim(integrated_atac); dim(integrated_rna)
all(colnames(sc)==colnames(integrated_atac))
all(colnames(integrated_atac)==colnames(integrated_rna))

# 6. Assign the rna_integrated_pca and rna_integrated_umap to sc ----
integrated_rna@reductions$rna_integrated_pca@assay.used <- "RNA"
sc[["rna_integrated_pca"]] <- SeuratObject::CreateDimReducObject(embeddings = integrated_rna@reductions$rna_integrated_pca@cell.embeddings, 
                                                                 loadings = integrated_rna@reductions$rna_integrated_pca@feature.loadings,
                                                                 projected = integrated_rna@reductions$rna_integrated_pca@feature.loadings.projected,
                                                                 assay = integrated_rna@reductions$rna_integrated_pca@assay.used,
                                                                 stdev = integrated_rna@reductions$rna_integrated_pca@stdev,
                                                                 key = integrated_rna@reductions$rna_integrated_pca@key, 
                                                                 jackstraw = integrated_rna@reductions$rna_integrated_pca@jackstraw, 
                                                                 misc = integrated_rna@reductions$rna_integrated_pca@misc)
#check cell matching
all(colnames(sc) == rownames(sc@reductions$rna_integrated_pca@cell.embeddings))
all(colnames(sc) == rownames(sc@reductions$rna_integrated_pca))

integrated_rna@reductions$rna_integrated_umap@assay.used <- "RNA"
sc[["rna_integrated_umap"]] <- SeuratObject::CreateDimReducObject(embeddings = integrated_rna@reductions$rna_integrated_umap@cell.embeddings, 
                                                                  loadings = integrated_rna@reductions$rna_integrated_umap@feature.loadings,
                                                                  projected = integrated_rna@reductions$rna_integrated_umap@feature.loadings.projected,
                                                                  assay = integrated_rna@reductions$rna_integrated_umap@assay.used,
                                                                  stdev = integrated_rna@reductions$rna_integrated_umap@stdev,
                                                                  key = integrated_rna@reductions$rna_integrated_umap@key, 
                                                                  jackstraw = integrated_rna@reductions$rna_integrated_umap@jackstraw, 
                                                                  misc = integrated_rna@reductions$rna_integrated_umap@misc, 
                                                                  global = TRUE)

#check cell matching
all(colnames(sc) == rownames(sc@reductions$rna_integrated_umap@cell.embeddings))
all(colnames(sc) == rownames(sc@reductions$rna_integrated_umap))

#Assign rna_integrated_cluster to sc@metadata
all(rownames(sc@meta.data) == rownames(integrated_rna@meta.data))
sc@meta.data %>% mutate(rna_integrated_cluster = integrated_rna@meta.data$rna_integrated_cluster) -> sc@meta.data

DimPlot(sc, group.by = "rna_integrated_cluster", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "group", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE)
rm(integrated_rna)

# 7. Assign the atac_integrated_lsi and atac_integrated_umap to sc ----
integrated_atac@reductions$atac_integrated_lsi@assay.used <- "ATAC"
sc[["atac_integrated_lsi"]] <- SeuratObject::CreateDimReducObject(embeddings = integrated_atac@reductions$atac_integrated_lsi@cell.embeddings, 
                                                                  loadings = integrated_atac@reductions$atac_integrated_lsi@feature.loadings,
                                                                  projected = integrated_atac@reductions$atac_integrated_lsi@feature.loadings.projected,
                                                                  assay = integrated_atac@reductions$atac_integrated_lsi@assay.used,
                                                                  stdev = integrated_atac@reductions$atac_integrated_lsi@stdev,
                                                                  key = integrated_atac@reductions$atac_integrated_lsi@key, 
                                                                  jackstraw = integrated_atac@reductions$atac_integrated_lsi@jackstraw, 
                                                                  misc = integrated_atac@reductions$atac_integrated_lsi@misc)
#check cell matching
all(colnames(sc) == rownames(sc@reductions$atac_integrated_lsi@cell.embeddings))
all(colnames(sc) == rownames(sc@reductions$atac_integrated_lsi))

integrated_atac@reductions$atac_integrated_umap@assay.used <- "ATAC"
sc[["atac_integrated_umap"]] <- SeuratObject::CreateDimReducObject(embeddings = integrated_atac@reductions$atac_integrated_umap@cell.embeddings, 
                                                                   loadings = integrated_atac@reductions$atac_integrated_umap@feature.loadings,
                                                                   projected = integrated_atac@reductions$atac_integrated_umap@feature.loadings.projected,
                                                                   assay = integrated_atac@reductions$atac_integrated_umap@assay.used,
                                                                   stdev = integrated_atac@reductions$atac_integrated_umap@stdev,
                                                                   key = integrated_atac@reductions$atac_integrated_umap@key, 
                                                                   jackstraw = integrated_atac@reductions$atac_integrated_umap@jackstraw, 
                                                                   misc = integrated_atac@reductions$atac_integrated_umap@misc, 
                                                                   global = TRUE)

#check cell matching
all(colnames(sc) == rownames(sc@reductions$atac_integrated_umap@cell.embeddings))
all(colnames(sc) == rownames(sc@reductions$atac_integrated_umap))

#Assign atac_integrated_cluster to sc@metadata
all(rownames(sc@meta.data) == rownames(integrated_atac@meta.data))
sc@meta.data %>% mutate(atac_integrated_cluster = integrated_atac@meta.data$atac_integrated_cluster) -> sc@meta.data

DimPlot(sc, group.by = "atac_integrated_cluster", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "group", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)
DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)

#rm(integrated_rna)
rm(integrated_atac)

# 8. ATAC and RNA (Joint clustering by Concensus concensus PCA, MultiCCA and WNN) ----

library(Destin2)
DefaultAssay(sc) <- "ATAC"
reductions <- list("rna_integrated_pca", "atac_integrated_lsi")
reduction_dims <- list(1:50, 2:50)

sc <- GetConsensusPCA(sc, reduction.list=reductions,
                      dims.list = reduction_dims,
                      reduction.name = 'consensus.pca',
                      reduction.key = "consensuspca_", 
                      assay = "ATAC", n.consensus.pc = 50)

sc <- GetMultiCCA(sc, reduction.list=reductions,
                  dims.list = reduction_dims,
                  reduction.name = 'multicca.pca',
                  reduction.key = "multiccapca_",
                  assay = "ATAC", n.cca = 50)

sc <- FindMultiModalNeighbors(sc, reduction.list=reductions,
                              dims.list = reduction_dims, verbose=TRUE, 
                              knn.graph.name = "integrated_wknn", 
                              snn.graph.name = "integrated_wsnn", 
                              weighted.nn.name = "integrated_weighted.nn")

sc <- RunUMAP(sc, dims = 1:50, reduction='consensus.pca', 
              reduction.name='consensus.umap', reduction.key = 'consensusumap_', verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:50, reduction='multicca.pca', 
              reduction.name='multicca.umap', reduction.key = 'multiccaumap_', verbose = FALSE)
sc <- RunUMAP(sc, nn.name = "integrated_weighted.nn", reduction.name = "integrated.wnn.umap", reduction.key = "integratedwnnUMAP_")

DefaultAssay(sc)='ATAC'
sc <- FindNeighbors(object = sc, reduction = 'consensus.pca', dims = 1:50, verbose = FALSE, graph.name=c('consensuspca_nn','consensuspca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='consensuspca_snn')
sc@meta.data %>% mutate(consensus_cluster = seurat_clusters) -> sc@meta.data

sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn')
sc@meta.data %>% mutate(cc_cluster = seurat_clusters) -> sc@meta.data

sc <- FindClusters(sc, graph.name = "integrated_wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
sc@meta.data %>% mutate(wnn_cluster = seurat_clusters) -> sc@meta.data

DimPlot(sc, group.by = "predicted.id", label = T, reduction = "consensus.umap")
DimPlot(sc, group.by = "ZT", label = T, reduction = "consensus.umap")
DimPlot(sc, group.by = "group", label = T, reduction = "consensus.umap")
DimPlot(sc, group.by = "consensus_cluster", label = T, reduction = "consensus.umap")

DimPlot(sc, group.by = "predicted.id", label = T, reduction = "multicca.umap")
DimPlot(sc, group.by = "ZT", label = T, reduction = "multicca.umap")
DimPlot(sc, group.by = "group", label = T, reduction = "multicca.umap")
DimPlot(sc, group.by = "cc_cluster", label = T, reduction = "multicca.umap")

DimPlot(sc, group.by = "predicted.id", label = T, reduction = "integrated.wnn.umap")
DimPlot(sc, group.by = "ZT", label = T, reduction = "integrated.wnn.umap")
DimPlot(sc, group.by = "group", label = T, reduction = "integrated.wnn.umap")
DimPlot(sc, group.by = "wnn_cluster", label = T, reduction = "integrated.wnn.umap")
####

DimPlot(sc, group.by = "predicted.id", label = T, reduction = "consensus.umap") + ggtitle("Concensus PCA") -> p1
DimPlot(sc, group.by = "predicted.id", label = T, reduction = "multicca.umap") + ggtitle("MultiCCA") -> p2
ggarrange(plotlist = list(p1, p2), nrow = 2, ncol = 1, legend = "bottom", common.legend = T)

# 9. Generate Gene Activity ----
DefaultAssay(sc) <- "ATAC"
sc@assays[["ATAC"]]@fragments[[1]]@path <- "~/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
#gene.activities <- GeneActivity(sc, extend.upstream = 2000, biotypes = NULL)
gene.activities <- GeneActivity(sc, extend.upstream = 2000)

sc[['gene_activity']] <- CreateAssayObject(counts = gene.activities)
sc <- NormalizeData(
  object = sc,
  assay = 'gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(sc$nCount_gene_activity)
)

DefaultAssay(sc) <- 'gene_activity'

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc, features = rownames(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc), 
             reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:50, reduction='activity.pca', 
              reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

DimPlot(sc, group.by = "predicted.id", label = T, reduction = "activity.umap")
DimPlot(sc, group.by = "group", label = T, reduction = "activity.umap")

## 9.1 Gene activity integration ----
sc_splited = SplitObject(sc, split.by = "group")

for (i in 1:length(sc_splited)) {
  sc_splited[[i]] <- SCTransform(sc_splited[[i]], vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", 
                                                                      "nCount_gene_activity", "nFeature_gene_activity"), assay="RNA", 
                                 verbose = TRUE, return.only.var.genes = FALSE, 
                                 new.assay.name = "gene_activity_SCT")
}

v.common <- rownames(sc_splited[[1]]$gene_activity_SCT@data)
for (  i in 1:length(sc_splited)){
  if (i == 1){
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(v.common))) 
  }else{
    v.common <- intersect(v.common, rownames(sc_splited[[i]]$gene_activity_SCT@data))
    print(sprintf("condtion %s number of gene in SCT matrix: %s", i, length(rownames(sc_splited[[i]]$gene_activity_SCT@data))))
  }
} 
length( v.common );
cat( "length(v.common)=", length(v.common), "\n" )

v.inte.features <- vector()
if ( 1 ){                                     # find variable features/genes common to all
  v.var.common <- sc_splited[[1]]@assays$gene_activity_SCT@var.features
  print(sprintf("condition %s number of SCT var.features: %s", 1, length(v.var.common)))
  for (i in 2:length(sc_splited)){
    v.var.common <- c(v.var.common, sc_splited[[i]]@assays$gene_activity_SCT@var.features)
    print(sprintf("condition %s number of SCT var.features: %s", i, length(sc_splited[[i]]@assays$gene_activity_SCT@var.features)))
  }
  length(unique(v.var.common))
  v.var.common <- names(which(table(v.var.common) == length(sc_splited)))
  
  cat("length(v.var.common)=", length(v.var.common), "\n")
  length(v.var.common)
  v.inte.features <- v.var.common
} else {
  v.inte.features <- SelectIntegrationFeatures(object.list = sc_splited, nfeatures = 3000)
}

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
k.weight.use = 100

sc_splited.anchors <- FindIntegrationAnchors(object.list = sc_splited, dims = 1:nDims, 
                                             assay=rep( "gene_activity_SCT", length(sc_splited)), 
                                             anchor.features=v.inte.features, k.anchor=k.anchor.use, 
                                             k.filter=k.filter.use, k.score=k.score.use, verbose=TRUE )

sc_integrated <- IntegrateData(anchorset = sc_splited.anchors, dims = 1:nDims, features.to.integrate=v.common, 
                               k.weight = k.weight.use, verbose = T, new.assay.name = "gene_activity_integrated")

VariableFeatures(sc_integrated, assay="gene_activity_integrated") <- v.inte.features

DefaultAssay(sc_integrated ) <- "gene_activity_integrated"
sc_integrated <- ScaleData(sc_integrated, verbose=TRUE, assay="gene_activity_integrated", features=rownames(sc_integrated$gene_activity_integrated@data))

sc_integrated <- RunPCA(sc_integrated, assay="gene_activity_integrated", npcs = nDims, verbose = FALSE, 
                        reduction.name = "gene_activity_integrated_pca", reduction.key = "geneactivityintegratedpc_")
sc_integrated <- RunUMAP(sc_integrated, reduction = "gene_activity_integrated_pca", dims = 1:nDims, 
                         reduction.name = "gene_activity_integrated_umap", 
                         reduction.key = "geneactivityintegratedumap_")
sc_integrated <- FindNeighbors(sc_integrated, assay="gene_activity_integrated", reduction = "gene_activity_integrated_pca", 
                               dims = 1:nDims, force.recalc=T, graph.name = c("integrated_gene_activity_nn", "integrated_gene_activity_snn"))
sc_integrated <- FindClusters(sc_integrated, assay="gene_activity_integrated", resolution = 1.0, graph.name = "integrated_gene_activity_snn")

sc_integrated@meta.data %>% mutate(gene_activity_cluster = seurat_clusters) -> sc_integrated@meta.data

DimPlot(sc_integrated, group.by = "gene_activity_cluster", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "celltype", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "group", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE)

table(sc_integrated$group, sc_integrated$gene_activity_cluster)
table(sc_integrated$celltype, sc_integrated$gene_activity_cluster)

table(sc_integrated$group, sc_integrated$gene_activity_cluster) %>% 
  as.data.frame() %>% 
  ggplot(., aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")

table(sc_integrated$celltype, sc_integrated$gene_activity_cluster) %>% 
  as.data.frame() %>% 
  ggplot(., aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="fill", stat="identity")

ncol(sc) == ncol(sc_integrated)
all(colnames(sc) == colnames(sc_integrated))

sc_integrated@reductions$gene_activity_integrated_pca@assay.used <- "gene_activity"
sc_gene_activity <- sc_integrated
rm(sc_integrated)
#all(colnames(sc) == colnames(sc_gene_activity))
sc[["gene_activity_integrated_pca"]] <- SeuratObject::CreateDimReducObject(embeddings = sc_gene_activity@reductions$gene_activity_integrated_pca@cell.embeddings, 
                                                                           loadings = sc_gene_activity@reductions$gene_activity_integrated_pca@feature.loadings,
                                                                           projected = sc_gene_activity@reductions$gene_activity_integrated_pca@feature.loadings.projected,
                                                                           assay = sc_gene_activity@reductions$gene_activity_integrated_pca@assay.used,
                                                                           stdev = sc_gene_activity@reductions$gene_activity_integrated_pca@stdev,
                                                                           key = sc_gene_activity@reductions$gene_activity_integrated_pca@key, 
                                                                           jackstraw = sc_gene_activity@reductions$gene_activity_integrated_pca@jackstraw, 
                                                                           misc = sc_gene_activity@reductions$gene_activity_integrated_pca@misc)
all(colnames(sc) == rownames(sc@reductions$gene_activity_integrated_pca))
all(colnames(sc) == rownames(sc@reductions$gene_activity_integrated_pca@cell.embeddings))

sc_gene_activity@reductions$gene_activity_integrated_umap@assay.used <- "gene_activity"
sc[["gene_activity_integrated_umap"]] <- SeuratObject::CreateDimReducObject(embeddings = sc_gene_activity@reductions$gene_activity_integrated_umap@cell.embeddings, 
                                                                            loadings = sc_gene_activity@reductions$gene_activity_integrated_umap@feature.loadings,
                                                                            projected = sc_gene_activity@reductions$gene_activity_integrated_umap@feature.loadings.projected,
                                                                            assay = sc_gene_activity@reductions$gene_activity_integrated_umap@assay.used,
                                                                            stdev = sc_gene_activity@reductions$gene_activity_integrated_umap@stdev,
                                                                            key = sc_gene_activity@reductions$gene_activity_integrated_umap@key, 
                                                                            jackstraw = sc_gene_activity@reductions$gene_activity_integrated_umap@jackstraw, 
                                                                            misc = sc_gene_activity@reductions$gene_activity_integrated_umap@misc, 
                                                                            global = TRUE)
all(colnames(sc) == rownames(sc@reductions$gene_activity_integrated_umap))
all(colnames(sc) == rownames(sc@reductions$gene_activity_integrated_umap@cell.embeddings))

DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE)

# 10. ChromVar ----
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)

genome=BSgenome.Mmusculus.UCSC.mm10
ensdb=EnsDb.Mmusculus.v79
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
DefaultAssay(sc)='ATAC'
sc <- AddMotifs(
  object = sc,
  genome = genome,
  pfm = pfm
)

# Computing motif activities
sc <- RunChromVAR(
  object = sc,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(sc) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
#Choice of 50 PCs here is arbitrary
pca_motif<-prcomp(t(sc@assays[["MOTIF"]]@data))$x[,1:50]
sc[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
sc <- RunUMAP(sc, reduction = "motif.pca", dims = 1:nDims, reduction.name = "motif_umap")

DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "motif_umap", repel = TRUE)
DimPlot(sc, group.by = "group", label = TRUE, reduction = "motif_umap", repel = TRUE)

# 10.1 ChoromVar Integration ----
sc_splited = SplitObject(sc, split.by = "group")
sc_splited[[1]]

nDims = 50
k.anchor.use = 5
k.filter.use = 199
k.score.use = 30
k.weight.use = 100

sc_splited.anchors <- FindIntegrationAnchors(object.list = sc_splited, dims = 1:nDims, 
                                             assay=rep( "MOTIF", length(sc_splited)), 
                                             k.anchor=k.anchor.use, 
                                             k.filter=k.filter.use, k.score=k.score.use, verbose=TRUE )

save(sc_splited.anchors, nDims, v.common, k.weight.use, file = "~/Downloads/sc_splited.anchors.chromVar.rda")
rm(sc_splited.anchors, sc_splited)

#run in server
sc_integrated <- IntegrateData(anchorset = sc_splited.anchors, dims = 1:nDims, 
                               k.weight = k.weight.use, verbose = T, new.assay.name = "MOTIF_integrated")

pca_motif<-prcomp(t(sc_integrated@assays[["MOTIF_integrated"]]@data))$x[,1:50]
sc_integrated[["motif_integrated.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifintegratedpca_", assay = "MOTIF_integrated")
sc_integrated <- RunUMAP(sc_integrated, reduction = "motif_integrated.pca", dims = 1:nDims, reduction.name = "motif_integrated.umap",
                         reduction.key = "motifintegratedumap_")
sc_integrated <- FindNeighbors(sc_integrated, assay="MOTIF_integrated", reduction = "motif_integrated.pca", 
                               dims = 1:nDims, force.recalc=T, graph.name = c("integratedmotif_nn", "integratedmotif_snn"))
sc_integrated <- FindClusters(sc_integrated, assay="MOTIF_integrated", resolution = 1.0, graph.name = "integratedmotif_snn")

DimPlot(sc_integrated, group.by = "celltype", label = TRUE, reduction = "motif_integrated.umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "group", label = TRUE, reduction = "motif_integrated.umap", repel = TRUE)
DimPlot(sc_integrated, group.by = "motif_integrated_cluster", label = TRUE, reduction = "motif_integrated.umap", repel = TRUE)

sc_motif = sc_integrated
rm(sc_integrated)
sc_motif@reductions$motif_integrated.pca@assay.used <- "MOTIF"
sc[["motif_integrated_pca"]] <- SeuratObject::CreateDimReducObject(embeddings = sc_motif@reductions$motif_integrated.pca@cell.embeddings, 
                                                                   loadings = sc_motif@reductions$motif_integrated.pca@feature.loadings,
                                                                   projected = sc_motif@reductions$motif_integrated.pca@feature.loadings.projected,
                                                                   assay = sc_motif@reductions$motif_integrated.pca@assay.used,
                                                                   stdev = sc_motif@reductions$motif_integrated.pca@stdev,
                                                                   key = sc_motif@reductions$motif_integrated.pca@key, 
                                                                   jackstraw = sc_motif@reductions$motif_integrated.pca@jackstraw, 
                                                                   misc = sc_motif@reductions$motif_integrated.pca@misc)

#all(colnames(sc) == rownames(sc@reductions$motif_integrated_pca))
#all(colnames(sc) == rownames(sc@reductions$motif_integrated_pca@cell.embeddings))

sc_motif@reductions$motif_integrated.umap@assay.used <- "MOTIF"
sc[["motif_integrated_umap"]] <- SeuratObject::CreateDimReducObject(embeddings = sc_motif@reductions$motif_integrated.umap@cell.embeddings, 
                                                                    loadings = sc_motif@reductions$motif_integrated.umap@feature.loadings,
                                                                    projected = sc_motif@reductions$motif_integrated.umap@feature.loadings.projected,
                                                                    assay = sc_motif@reductions$motif_integrated.umap@assay.used,
                                                                    stdev = sc_motif@reductions$motif_integrated.umap@stdev,
                                                                    key = sc_motif@reductions$motif_integrated.umap@key, 
                                                                    jackstraw = sc_motif@reductions$motif_integrated.umap@jackstraw, 
                                                                    misc = sc_motif@reductions$motif_integrated.umap@misc, 
                                                                    global = TRUE)

options(ggrepel.max.overlaps = Inf)
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "motif_integrated_umap", repel = TRUE) + 
  ggtitle("chromVar Motif Score") -> p_motif
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE) + 
  ggtitle("Gene Activity") -> p_gene_activity
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE) + 
  ggtitle("scRNA-seq") -> p_rna
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE) + 
  ggtitle("scATAC-seq") -> p_atac
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "multicca.umap", repel = TRUE) + 
  ggtitle("MultiCCA") -> p_cca
ggpubr::ggarrange(p_rna, p_atac, p_cca, p_gene_activity, p_motif, ncol = 3, nrow = 2, common.legend = T, legend = "right")

options(ggrepel.max.overlaps = Inf)
DimPlot(sc, group.by = "group", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE)

#save.image(file = "~/Dropbox/singulomics/rda/5_6_umap_keep_other_cells.rda")
saveRDS(sc, "~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds")


#Add ATAC_TSS assay ----
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc

#Generate ATAC_TSS matrix ----
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
fpath = "/Users/chunyiptong/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
cells_ = colnames(sc)
fragments_ <- CreateFragmentObject(fpath, cells = cells_)
upstream = 2000 
downstream = 2000
gene.tss.ref = GenomicRanges::promoters(gene.ref, upstream = upstream, downstream = downstream) 

atac_counts_df <- Signac::FeatureMatrix(fragments = fragments_, 
                                        features = gene.tss.ref, 
                                        cells = cells_)

atac_counts_df %>% rownames() -> seq_
atac_counts_df %>% as.data.frame() %>% mutate(seq_ = seq_) -> atac_counts_df

gene.tss.ref %>% as.data.frame() %>% 
  dplyr::select(seqnames, start, end, gene_name) %>% 
  mutate(seq_ = sprintf("%s-%s-%s", seqnames, start, end)) %>% 
  dplyr::select(seq_, gene_name) %>% 
  {right_join(x = .,
              y = atac_counts_df,
              by = "seq_")} -> atac_counts_df

atac_counts_df %>% distinct() -> atac_counts_df
dim(atac_counts_df)
rm(seq_)

atac_counts_df$gene_name -> rownames_
rownames(atac_counts_df) <- rownames_
atac_counts_df %>% dplyr::select(-seq_, -gene_name) -> atac_counts_df 
rm(rownames_)
ATAC_counts = as.matrix(atac_counts_df)
rm(atac_counts_df)

sc[['ATAC_TSS']] <- CreateAssayObject(counts = ATAC_counts)

##Plot number of cells in different groups of each celltype ----
sc@meta.data$celltype %>% unique() %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    sc@meta.data %>% dplyr::filter(celltype == x) %>% 
      dplyr::select(celltype, group) %>% 
      group_by(group) %>% 
      group_map(function(x_,y_){
        nrow(x_) -> n_
        data_frame(celltype = x, group = y_$group, n_cell = n_)
      }, .keep = T) %>% 
      do.call(rbind, .) %>% 
      ggplot(aes(x = group, y = n_cell, fill = group)) + 
      geom_bar(stat = "identity") + 
      geom_text(aes(label=n_cell), position=position_dodge(width=0.9), vjust=-0.25) +
      #      theme_minimal() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      ggtitle(x) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.title.x = element_blank()) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }) -> p_list #Supp_Fig_4B ----

ggpubr::ggarrange(plotlist = p_list, ncol = 2, nrow = 4, common.legend = T, legend = "top")

c("B cells", "Cholangiocytes", "Endothelial cells", "Fibroblasts", "Hepatocytes", "Kupffer cells", "T cells") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    celltype_ = x
    sc@meta.data %>% dplyr::filter(celltype == celltype_) %>% .$prediction.score.max %>% 
      as.data.frame() %>% "colnames<-"(., "Prediction_score") %>% 
      "colnames<-"(., "Prediction_score") %>% 
      ggplot(aes(x = Prediction_score)) + 
#      geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") + 
      geom_histogram(fill = "skyblue", color = "black") + 
      theme_classic() -> p
    p + ggtitle(celltype_) -> p
  }) -> p_list #Supp_Fig_4C ----
patchwork::wrap_plots(p_list, ncol = 7)

saveRDS(sc, "~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds")

### cell reclustering without KO cells ---- 
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc$group %>% table()
sc@meta.data %>% dplyr::filter(!grepl("KO", group)) %>% rownames() -> cell_sele
sc = sc[,cell_sele]

#RNA
DefaultAssay(sc) <- "RNA"
sc <- SCTransform(sc, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#ATAC
DefaultAssay(sc) <- "ATAC"
sc <- RunTFIDF(sc)
sc <- FindTopFeatures(sc, min.cutoff = 'q0')
sc <- RunSVD(sc)
sc <- RunUMAP(sc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Gene activity
DefaultAssay(sc) <- 'gene_activity'

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
sc <- ScaleData(sc, features = rownames(sc))
sc <- RunPCA(sc, features = VariableFeatures(sc), 
             reduction.name = 'activity.pca',  reduction.key = 'activitypca_', verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:50, reduction='activity.pca', 
              reduction.name='activity.umap', reduction.key = 'activityumap_', verbose = FALSE)

#MOTIF score
DefaultAssay(sc)='ATAC'
# Computing motif activities
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)

genome=BSgenome.Mmusculus.UCSC.mm10
ensdb=EnsDb.Mmusculus.v79
sc <- RunChromVAR(
  object = sc,
  genome = genome,
  new.assay.name = 'MOTIF'
)

DefaultAssay(sc) <- 'MOTIF'
#Using PCA dimensional reduction on motif modality and put it into Seurat Object
#Choice of 50 PCs here is arbitrary
nDims = 50
pca_motif<-prcomp(t(sc@assays[["MOTIF"]]@data))$x[,1:50]
sc[["motif.pca"]] <- CreateDimReducObject(embeddings = pca_motif, key = "motifpca_", assay = "MOTIF")
sc <- RunUMAP(sc, reduction = "motif.pca", dims = 1:nDims, reduction.name = "motif_umap")

#Supp_Fig_3 and Fig_1B ----
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "umap.rna", repel = TRUE) + 
  ggtitle("scRNA-seq") -> p1
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "umap.atac", repel = TRUE) + 
  ggtitle("scATAC-seq") -> p2
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "activity.umap", repel = TRUE) + 
  ggtitle("Gene Activity") -> p3
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "motif_umap", repel = TRUE) + 
  ggtitle("MOTIF Score") -> p4
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, legend = "none") -> p #18.03X4.35

DimPlot(sc, group.by = "group", label = TRUE, reduction = "umap.rna", repel = TRUE) + 
  ggtitle("scRNA-seq") -> p1
DimPlot(sc, group.by = "group", label = TRUE, reduction = "umap.atac", repel = TRUE) + 
  ggtitle("scATAC-seq") -> p2
DimPlot(sc, group.by = "group", label = TRUE, reduction = "activity.umap", repel = TRUE) + 
  ggtitle("Gene Activity") -> p3
DimPlot(sc, group.by = "group", label = TRUE, reduction = "motif_umap", repel = TRUE) + 
  ggtitle("MOTIF Score") -> p4
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, legend = "none") -> p

DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE) + 
  ggtitle("scRNA-seq") -> p1
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE) + 
  ggtitle("scATAC-seq") -> p2
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE) + 
  ggtitle("Gene Activity") -> p3
DimPlot(sc, group.by = "celltype", label = TRUE, reduction = "motif_integrated_umap", repel = TRUE) + 
  ggtitle("MOTIF Score") -> p4
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, legend = "none") -> p

DimPlot(sc, group.by = "group", label = TRUE, reduction = "rna_integrated_umap", repel = TRUE) + 
  ggtitle("scRNA-seq") -> p1
DimPlot(sc, group.by = "group", label = TRUE, reduction = "atac_integrated_umap", repel = TRUE) + 
  ggtitle("scATAC-seq") -> p2
DimPlot(sc, group.by = "group", label = TRUE, reduction = "gene_activity_integrated_umap", repel = TRUE) + 
  ggtitle("Gene Activity") -> p3
DimPlot(sc, group.by = "group", label = TRUE, reduction = "motif_integrated_umap", repel = TRUE) + 
  ggtitle("MOTIF Score") -> p4
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, legend = "none") -> p
