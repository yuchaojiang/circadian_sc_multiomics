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

load('rda/1_Read_Data.rda')

#######################################
## read, feature, and mitochondria
#######################################

VlnPlot(sc, features = c("nCount_ATAC", "nFeature_ATAC","TSS.enrichment","nucleosome_signal",
                         "nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 4,
        log = TRUE, pt.size = 0) + NoLegend()

VlnPlot(sc, group.by='group',features = c("nCount_ATAC", "nFeature_ATAC","TSS.enrichment","nucleosome_signal",
                         "nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 4,
        log = TRUE, pt.size = 0) + NoLegend()

##Draw Violin plots (nCount_ATAC, nFeature_ATAC, TSS.enrichment, nucleosome_signal, nCount_RNA, nFeature_RNA, percent.mt) before QC ----
c("nCount_ATAC", "nFeature_ATAC","TSS.enrichment","nucleosome_signal",
  "nCount_RNA", "nFeature_RNA","percent.mt") %>% 
#  .[1] %>% 
  map(function(x_){
    sc[[x_]] -> df_
    df_ %>% "colnames<-"("value") %>% dplyr::mutate(group = x_) -> df_
    head(df_)
    ggplot(df_, aes(x = group, y = value, fill = group)) + 
    geom_violin() + 
    xlab(NULL) + ylab(NULL) + 
    theme_classic() + theme(legend.position = "none") -> p
#    if (x_ != "percent.mt"){
      p + scale_y_continuous(trans = "log10") -> p
#    }
    return(p)
  }) -> p_list
patchwork::wrap_plots(p_list, ncol = 4) -> p_pre_QC #Fig_S2A ----

sc <- subset(
  x = sc,
  subset = nCount_ATAC < 1e5 &
    nCount_ATAC > 1e3 &
    nFeature_ATAC > 500 &
    TSS.enrichment >3 &
    TSS.enrichment <20 &
    nucleosome_signal < 2 &
    nCount_RNA < 1e5 &
    nCount_RNA > 1e3 &
    nFeature_RNA > 500 &
    percent.mt < 5
)

VlnPlot(sc, features = c("nCount_ATAC", "nFeature_ATAC","TSS.enrichment","nucleosome_signal",
                         "nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 4,
        log = TRUE, pt.size = 0) + NoLegend()

##Draw Violin plots (nCount_ATAC, nFeature_ATAC, TSS.enrichment, nucleosome_signal, nCount_RNA, nFeature_RNA, percent.mt) after QC ----
c("nCount_ATAC", "nFeature_ATAC","TSS.enrichment","nucleosome_signal",
  "nCount_RNA", "nFeature_RNA","percent.mt") %>% 
  #  .[1] %>% 
  map(function(x_){
    sc[[x_]] -> df_
    df_ %>% "colnames<-"("value") %>% dplyr::mutate(group = x_) -> df_
    head(df_)
    ggplot(df_, aes(x = group, y = value, fill = group)) + 
      geom_violin() + 
      xlab(NULL) + ylab(NULL) + 
      theme_classic() + theme(legend.position = "none") -> p
    #    if (x_ != "percent.mt"){
    p + scale_y_continuous(trans = "log10") -> p
    #    }
    return(p)
  }) -> p_list

patchwork::wrap_plots(p_list, ncol = 4) -> p_after_QC #Fig_S2B ----

median(sc$nCount_ATAC)
median(sc$nCount_RNA)
median(sc$nFeature_ATAC)
median(sc$nFeature_RNA)

VlnPlot(sc, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC"), ncol = 4,
        log = TRUE, pt.size = 0) + NoLegend()


dim(sc@assays$RNA); length(gene.ref)
dim(sc@assays$ATAC); length(peak.ref)

table(sc$group)

# #######################################
# ## read, feature, and mitochondria
# #######################################
# 
# par(mfrow=c(1,2))
# 
# plot(sc$nCount_RNA,sc$nFeature_RNA, pch=16, cex=0.6, col=2, main='RNA QC')
# spl <- smooth.spline(sc$nCount_RNA, sc$nFeature_RNA)
# spl.pred = predict(spl, sc$nCount_RNA)$y
# points(sc$nCount_RNA, spl.pred, cex=0.6, pch=16, col='blue')
# 
# spl.resid=sc$nFeature_RNA-spl.pred
# spl.filter=spl.resid < median(spl.resid)-6*mad(spl.resid)  | spl.resid > median(spl.resid)+6*mad(spl.resid)
# sc=sc[,!spl.filter]
# points(sc$nCount_RNA,sc$nFeature_RNA, pch=16, cex=0.6)
# 
# 
# plot(sc$nCount_ATAC,sc$nFeature_ATAC, pch=16, cex=0.6, col=2, main='ATAC QC')
# spl <- smooth.spline(sc$nCount_ATAC, sc$nFeature_ATAC)
# spl.pred = predict(spl, sc$nCount_ATAC)$y
# points(sc$nCount_ATAC, spl.pred, cex=0.6, pch=16, col='blue')
# 
# spl.resid=sc$nFeature_ATAC-spl.pred
# spl.filter=spl.resid < median(spl.resid)-6*mad(spl.resid)  | spl.resid > median(spl.resid)+6*mad(spl.resid)
# sc=sc[,!spl.filter]
# points(sc$nCount_ATAC,sc$nFeature_ATAC, pch=16, cex=0.6)
# 
# rm(spl, spl.resid, spl.pred, spl.filter)

#######################################
## Remove doublets
#######################################
set.seed(1234)
# https://github.com/plger/scDblFinder
sce.rna <- SingleCellExperiment(list(counts=sc@assays$RNA@counts))
sce.rna <- scDblFinder(sce.rna, samples = sc$group, dbr = 0.005)
table(sce.rna$scDblFinder.class, sc$group)

# https://plger.github.io/scDblFinder/articles/scATAC.html
sce.atac <- SingleCellExperiment(list(counts=sc@assays$ATAC@counts))
sce.atac <- scDblFinder(sce.atac, samples = sc$group, aggregateFeatures=TRUE, nfeatures=100, dbr = 0.005) # Using the top 1000 features
table(sce.atac$scDblFinder.class, sc$group)

sc=sc[,!(sce.rna$scDblFinder.class=='doublet' & sce.atac$scDblFinder.class =='doublet')]
rm(sce.atac, sce.rna)

#######################################
# At least ~0.1% of cells expressing the gene or have peak accessibility
#######################################
rna_counts=sc@assays$RNA@counts
atac_counts=sc@assays$ATAC@counts

cell.filter=(ncol(rna_counts)*0.001)
rna_counts_temp=rna_counts>0
gene.ref=gene.ref[rowSums(rna_counts_temp)>cell.filter ]
rna_counts=rna_counts[rowSums(rna_counts_temp)>cell.filter,]
rm(rna_counts_temp)

atac_counts_temp=atac_counts>0
peak.ref=peak.ref[rowSums(atac_counts_temp)>cell.filter]
atac_counts=atac_counts[rowSums(atac_counts_temp)>cell.filter,]
rm(atac_counts_temp, cell.filter)


#######################################
# Create new object
#######################################
sc.old=sc
sc <- CreateSeuratObject(counts = rna_counts)
rm(rna_counts)
sc[["percent.mt"]] = sc.old$percent.mt
sc[['ZT']]=sc.old$ZT
sc[['group']]=sc.old$group

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = "aggregate_analysis/atac_fragments.tsv.gz",
  annotation = annotations
)
sc[["ATAC"]] <- chrom_assay
rm(chrom_assay); rm(atac_counts); rm(annotations); rm(sc.old)


# RNA analysis
DefaultAssay(sc) <- "RNA"
sc <- SCTransform(sc, verbose = FALSE, return.only.var.genes = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sc) <- "ATAC"
sc <- RunTFIDF(sc)
sc <- FindTopFeatures(sc, min.cutoff = 'q0')
sc <- RunSVD(sc)
sc <- RunUMAP(sc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sc <- FindMultiModalNeighbors(sc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sc <- RunUMAP(sc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sc <- FindClusters(sc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

dim(sc@assays$RNA); length(gene.ref)
dim(sc@assays$ATAC); length(peak.ref)
table(sc$group)

# Some initial visualization
p1 <- DimPlot(sc, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA") + NoLegend()
p2 <- DimPlot(sc, reduction = "umap.rna", group.by='group', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc RNA") 
p3 <- DimPlot(sc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc ATAC")+ NoLegend()
p4 <- DimPlot(sc, reduction = "umap.atac", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X sc ATAC")
p5 <- DimPlot(sc, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")+ NoLegend()
p6 <- DimPlot(sc, reduction = "wnn.umap", group.by = "group", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1+p2
p3+p4

rm(p1,p2,p3,p4,p5,p6)

save.image(file='rda/2_QC.rda')