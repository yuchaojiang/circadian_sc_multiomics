library(tidyverse)
library(Seurat)
library(Signac)
library(TRIPOD)


# 1. load data ----
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
readRDS("~/Dropbox/singulomics/github_rda/Hepatocytes_cellnames.rds") -> hepatocytes_cells

colnames(sc) %>% length()
length(hepatocytes_cells)

sc <- sc[,hepatocytes_cells]
sc$celltype %>% table()

DimPlot(sc, group.by = "ZT", label = TRUE, reduction = "multicca.umap", repel = TRUE)
DefaultAssay(sc) <- "SCT"

# Choose resolution: do not need too many metacells, especially because we need to characterize the distributions of cells within the metacells
resolution=1
sc <- FindNeighbors(object = sc, reduction = 'multicca.pca', dims = 1:50, verbose = FALSE, graph.name=c('multicca_nn','multicca_snn'))
sc <- FindClusters(object = sc, algorithm = 1, verbose = FALSE, graph.name='multicca_snn',
                   resolution=resolution)
DimPlot(sc, group.by = "seurat_clusters", label = TRUE, reduction = "multicca.umap", repel = TRUE)

sc$group=droplevels(sc$group)
table(sc$group)

DimPlot(sc, group.by = "celltype", label = T, reduction = "multicca.umap") |
  DimPlot(sc, group.by = "seurat_clusters", label = T, reduction = "multicca.umap") 

# Need to construct metacells by ZT: each ZT has different metacells
table(paste0(sc$seurat_clusters, '_',sc$ZT)) %>% {.[names(.) %>% gtools::mixedsort()]}

# 0_ZT02  0_ZT06  0_ZT10  0_ZT14  0_ZT18  0_ZT22  1_ZT02  1_ZT06  1_ZT10  1_ZT14  1_ZT18  1_ZT22  2_ZT02  2_ZT06  2_ZT10  2_ZT14  2_ZT18  2_ZT22 
# 253     346     308     571     371     452     278     244     387     542     362     436     244     451     333     444     244     423 
# 3_ZT02  3_ZT06  3_ZT10  3_ZT14  3_ZT18  3_ZT22  4_ZT02  4_ZT06  4_ZT10  4_ZT14  4_ZT18  4_ZT22  5_ZT02  5_ZT06  5_ZT10  5_ZT14  5_ZT18  5_ZT22 
# 377     327     240     569     302     324     232     258     341     437     197     263     308     165     248     109     210     388 
# 6_ZT02  6_ZT06  6_ZT10  6_ZT14  6_ZT18  6_ZT22  7_ZT02  7_ZT06  7_ZT10  7_ZT14  7_ZT18  7_ZT22  8_ZT02  8_ZT06  8_ZT10  8_ZT14  8_ZT18  8_ZT22 
# 239      90     243      75     190     437     136     228     109     185     243     255     164     171     133     234     156     179 
# 9_ZT02  9_ZT06  9_ZT10  9_ZT14  9_ZT18  9_ZT22 10_ZT02 10_ZT06 10_ZT10 10_ZT14 10_ZT18 10_ZT22 
# 132     136      43      77      87      77      37      37      22      33      25      22 

length(table(paste0(sc$seurat_clusters, '_',sc$ZT)))
# 66

i = 0
levels(sc$seurat_clusters) %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    cluster_ = x
    sc@meta.data %>% dplyr::filter(seurat_clusters == cluster_) -> meta
    levels(sc$ZT) %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        ZT_ = x
        meta %>% dplyr::filter(ZT == ZT_) -> meta
        meta  %>% rownames_to_column("cellname") -> meta 
        meta[,c("cellname", "seurat_clusters")] -> meta
        meta$new_cluster = as.character(i)
        i <<- i+1
        return(meta)
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>% as_tibble() %>% column_to_rownames("cellname") -> meta_
meta_[rownames(sc@meta.data), ] -> meta_
sc$seurat_clusters <- meta_$new_cluster
factor(sc$seurat_clusters) %>% levels() %>% gtools::mixedsort() -> levels_
factor(sc$seurat_clusters, levels = levels_)  -> sc$seurat_clusters

DimPlot(sc, group.by = "seurat_clusters", label = T, reduction = "multicca.umap") 

# 2. Run TRIPOD
getObjectsForModelFit_custom = function(object, chr, biotype){
  DefaultAssay(object) <- "ATAC"
  transcripts.gr <- Signac:::CollapseToLongestTranscript(ranges = Annotation(object))
  
  if (biotype == "all"){
    transcripts.gr <- transcripts.gr 
  }else{
    transcripts.gr <- transcripts.gr[transcripts.gr$gene_biotype == 
                                       biotype]    
  }
  
  transcripts.gr <- transcripts.gr[as.character(seqnames(transcripts.gr)) %in% 
                                     chr]
  transcripts.gr <- sort(transcripts.gr)
  peaks.gr <- object@assays$ATAC@ranges
  motifxTF <- unlist(object@assays$ATAC@motifs@motif.names)
  motifxTF <- cbind(names(motifxTF), toupper(motifxTF))
  colnames(motifxTF) <- c("motif", "TF")
  peakxmotif <- object@assays$ATAC@motifs@data
  sel <- motifxTF[, 2] %in% toupper(rownames(object@assays$RNA))
  peakxmotif <- peakxmotif[, sel]
  motifxTF <- motifxTF[sel, ]
  motifxTF[, 2] <- rownames(object@assays$RNA)[match(motifxTF[, 
                                                              2], toupper(rownames(object@assays$RNA)))]
  genes <- intersect(transcripts.gr$gene_name, rownames(object@assays$SCT))
  DefaultAssay(object) <- "RNA"
  transcripts.gr <- transcripts.gr[match(genes, transcripts.gr$gene_name)]
  peakxmotif <- peakxmotif[, motifxTF[, 2] %in% genes]
  motifxTF <- motifxTF[motifxTF[, 2] %in% genes, ]
  result <- list(transcripts.gr = transcripts.gr, peaks.gr = peaks.gr, 
                 motifxTF = motifxTF, peakxmotif = peakxmotif)
  return(result)
}

tripod.sc <- getObjectsForModelFit_custom(object = sc, chr = paste0("chr", 1:19), biotype = "all")

transcripts.gr <- tripod.sc$transcripts.gr
peaks.gr <- tripod.sc$peaks.gr
motifxTF <- tripod.sc$motifxTF
peakxmotif <- tripod.sc$peakxmotif

filterSeuratObject_custom = function(object, tripod.object){
  genes <- tripod.object$transcripts.gr$gene_name
  motifxTF <- tripod.object$motifxTF
  object@assays$RNA <- subset(object@assays$RNA, features = match(genes, 
                              rownames(object@assays$RNA)))
  object@assays$SCT <- subset(object@assays$SCT, features = match(genes, 
                              rownames(object@assays$SCT)))
  object@assays$MOTIF <- subset(object@assays$MOTIF, 
                                   features = match(motifxTF[, 1], 
                                  rownames(object@assays$MOTIF)))
  return(object)
}
sc <- filterSeuratObject_custom(object = sc, tripod.object = tripod.sc)

# Recluster

metacell.sc <- getMetacellMatrices(object = sc,
                                    cluster.name = "seurat_clusters",
                                    min.num = 0
)

metacell.rna <- metacell.sc$rna
metacell.peak <- metacell.sc$peak
rm(metacell.sc)

sc$ZT=factor(sc$ZT)
color.sc <- getColors(object = sc, reduction = "multicca.umap", celltype.col.name = "ZT", cluster.col.name = "seurat_clusters")

metacell.celltype <- color.sc$metacell$celltype
metacell.celltype.col <- color.sc$metacell$color

DefaultAssay(sc) <- "SCT"
hvg.sc <- VariableFeatures(sc)

sc.RNA.genes = sc@assays$SCT@counts %>% rownames()

save(transcripts.gr, file = '~/Dropbox/singulomics/github_rda/TRIPOD/transcripts.gr.rda')
save(peaks.gr, file='~/Dropbox/singulomics/github_rda/TRIPOD/peaks.gr.rda')
save(motifxTF, file='~/Dropbox/singulomics/github_rda/TRIPOD/motifxTF.rda')
save(peakxmotif, file='~/Dropbox/singulomics/github_rda/TRIPOD/peakxmotif.rda')
save(metacell.rna, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
save(metacell.peak, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')
save(metacell.celltype, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.celltype.rda')
save(metacell.celltype.col, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.celltype.col.rda')
save(hvg.sc, file='~/Dropbox/singulomics/github_rda/TRIPOD/hvg.sc.rda')
save(sc.RNA.genes, file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.RNA.genes.rda')
save(sc, file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')
rm(list=ls())
####

# 3. Process TRIPOD outputs ----
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(patchwork)
library(BiocParallel)
library(dendextend)
library(TRIPOD)
library(tidyverse)

load(file = '~/Dropbox/singulomics/github_rda/TRIPOD/transcripts.gr.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/peaks.gr.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/motifxTF.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/peakxmotif.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.celltype.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.celltype.col.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/hvg.sc.rda')
load(file="~/Dropbox/singulomics/github_rda/TRIPOD/sc.RNA.genes.rda")
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')

##Replace TF in metacell.rna with motif score##
##Make metacell.motif_score_scaled##
sc@assays$MOTIF@data -> motif_score_mat
row.names(motif_score_mat) %>% 
#  .[1] %>% 
  purrr::map(function(motif_){
    motifxTF[motif_, 2]
  }) %>% purrr::reduce(., c) -> TF_gene
rownames(motif_score_mat) <- TF_gene

sc@meta.data -> sc_meta
levels(sc_meta$seurat_clusters) %>% 
#  .[1:2] %>% 
  purrr::map(function(cluster_){
    sc_meta %>% dplyr::filter(seurat_clusters == cluster_) %>% rownames() -> cells_
    motif_score_mat[, cells_] %>% rowMeans() %>% as.matrix() %>% "colnames<-"(., sprintf("metacell_%s", cluster_))
  }) %>% do.call(cbind, .) %>% t() %>% as.matrix() -> metacell.motif_score
metacell.motif_score %>% apply(., 2, function(x){scales::rescale(x = x, to = c(0,1))}) -> metacell.motif_score_scaled
save(metacell.motif_score, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.motif_score.rda')
save(metacell.motif_score_scaled, file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.motif_score_scaled.rda')

getXYMatrices_custom = function (gene.name, ext.upstream, ext.downstream = NULL, transcripts.gr, 
                                 peaks.gr, metacell.rna, metacell.peak, peakxmotif, motifxTF, 
                                 metacell.celltype = NULL, metacell.celltype.col = NULL) 
{
  if (is.null(ext.downstream)) 
    ext.downstream <- ext.upstream
  transcripts.ext.gr <- promoters(transcripts.gr, upstream = ext.upstream, 
                                  downstream = ext.downstream + 1)
  Yg <- metacell.rna[, gene.name]
  transcripts.gr.g <- transcripts.gr[transcripts.gr$gene_name == 
                                       gene.name]
  transcripts.ext.gr.g <- transcripts.ext.gr[transcripts.ext.gr$gene_name == 
                                               gene.name]
  peaks.gr.g <- subsetByOverlaps(peaks.gr, transcripts.ext.gr.g)
  Xt <- metacell.peak[, overlapsAny(peaks.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name == 
                                                                   gene.name]), drop = FALSE]
  peakxmotif.g <- peakxmotif[overlapsAny(peaks.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name == 
                                                                        gene.name]), , drop = FALSE]
  peakxmotif.g <- peakxmotif.g[, apply(peakxmotif.g, 2, sum) > 
                                 0, drop = FALSE]
  TF.g <- motifxTF[match(colnames(peakxmotif.g), motifxTF[, 
                                                          1]), 2]
  Yj <- metacell.motif_score_scaled[, TF.g, drop = FALSE]
  Yj_ <- metacell.rna[, TF.g, drop = FALSE]
  
  TF.filter1 <- apply(Yj, 2, sum) != 0
  TF.filter2 <- TF.g != gene.name
  TF.filter3 <- apply(Yj_, 2, sum) > 0 
  TF.filter <- TF.filter1 & TF.filter2 & TF.filter3
  Yj <- Yj[, TF.filter, drop = FALSE]
  TF.g <- TF.g[TF.filter, drop = FALSE]
  peakxmotif.g <- peakxmotif.g[, TF.filter, drop = FALSE]
  peakxmotif.g <- as.matrix(peakxmotif.g)
  nonzero.peakxmotif.g <- which(peakxmotif.g != 0, arr.ind = T)
  rm(TF.filter)
  peakxmotif.g <- as(peakxmotif.g, "dgCMatrix")
  results <- list(gene.name = gene.name, ext.upstream = ext.upstream, 
                  ext.downstream = ext.downstream, Yg = Yg, Xt = Xt, Yj = Yj, 
                  peakxmotif.g = peakxmotif.g, nonzero.peakxmotif.g = nonzero.peakxmotif.g, 
                  TF.g = TF.g, peaks.gr.g = peaks.gr.g, metacell.celltype = metacell.celltype, 
                  metacell.celltype.col = metacell.celltype.col)
  return(results)
}

genes_ = sc@assays$SCT@counts %>% rownames()
ext.upstream <- ext.downstream <- 2e5
xymats.list <- bplapply(
  genes_,
  getXYMatrices_custom,
  ext.upstream = ext.upstream,
  transcripts.gr = transcripts.gr,
  peaks.gr = peaks.gr,
  metacell.rna = metacell.rna,
  metacell.peak = metacell.peak,
  peakxmotif = peakxmotif,
  motifxTF = motifxTF,
  metacell.celltype = metacell.celltype,
  metacell.celltype.col = metacell.celltype.col
)
names(xymats.list) <- genes_

xymats.list %>% 
  map(function(x){
    x$peaks.gr.g %>% as.data.frame() %>% nrow() -> nrow_
    if (nrow_ > 0){
      x$gene
    }else{
      NULL
    }
  }) %>% purrr::reduce(., c) -> gene_
xymats.list <- xymats.list[gene_]

save(xymats.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.list_new_scaled.rda')
gc()

##Fit TRIPOD Xt##
safe_fitModel <- function(x) {
  tryCatch(
    fitModel(xymats = x, model.name = "TRIPOD", match.by = "Xt"),
    error = function(e) NULL
  )
}
xymats.tripod.Xt.list <- bplapply(xymats.list, safe_fitModel)
save(xymats.tripod.Xt.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Xt.list_new_scaled.rda')
gc()
####

##Fit TRIPOD Yj##
safe_fitModel <- function(x) {
  tryCatch(
    fitModel(xymats = x, model.name = "TRIPOD", match.by = "Yj"),
    error = function(e) NULL
  )
}
xymats.tripod.Yj.list <- bplapply(xymats.list, safe_fitModel)
save(xymats.tripod.Yj.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Yj.list_new.rda')
gc()
####

load("~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Yj.list_new_scaled.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Xt.list_new_scaled.rda")

fdr.thresh <- 0.01
list(
  TRIPOD.Yj = xymats.tripod.Yj.list_scaled,
  TRIPOD.Xt = xymats.tripod.Xt.list_scaled
) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    model_ = y
    print(model_)
    list_ = x
    c("positive", "negative") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        sign_ = x
        if (grepl("TRIPOD", model_)){
          c(1, 2) %>% 
            "names<-"(., sprintf("level_%s", .)) %>% 
            map(function(x){
              level_ = x
              getTrios(
                xymats.list = list_,
                fdr.thresh = fdr.thresh,
                sign = sign_,
                model.name = "TRIPOD",
                level = level_
              ) -> res_
              return(res_)
            })
        }else{
          getTrios(
            xymats.list = list_,
            fdr.thresh = fdr.thresh,
            sign = sign_,
            model.name = model_
          ) -> res_
          return(res_)
        }
      })
  }) -> list_res
list_trios_res_scaled = list_res ; rm(list_res)
save(list_trios_res_scaled, file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res_new_scaled.rda')
####

##Add addition information to the TRIPOD dataframe##
library(Seurat)
library(Signac)
library(tidyverse)

load(file='~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res_new.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.motif_score.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')

readRDS(file="~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
gene.ref %>% as.data.frame() -> gene.ref

sc@meta.data[,c("ZT", "seurat_clusters")] %>% distinct() %>% dplyr::arrange(seurat_clusters) %>% 
  mutate(metacell = sprintf("metacell_%s", seurat_clusters)) %>% 
  "rownames<-"(., 1:nrow(.)) -> metacell_ZT

source("~/Dropbox/singulomics/github/Calculate_HMP.R")
list(rna = metacell.rna, peak = metacell.peak, motif.score = metacell.motif_score) %>% 
  #list(peak = metacell.peak) %>% 
  map(function(x){
    df_ = x 
    df_meta = metacell_ZT %>% dplyr::arrange(ZT) %>% mutate(ZT = sprintf("%s_REP%s", ZT, 1:11))
    t(x) %>% as.data.frame() %>% .[, df_meta$metacell] -> df_
    colnames(df_) = df_meta$ZT
    #    head(df_) -> df_
    colnames(df_) %>% gsub("ZT(\\d+)_.+", "\\1", .) %>% as.integer() -> timepoint
    rownames(df_) -> genes
    df_ %>% as_tibble() -> df_
    cyclic_HMP(exp_matrix = df_, gene = genes, timepoints = timepoint) -> res_
    res_ %>% dplyr::select(matches("Gene|JTK|HR")) -> res_
    recal_cauchy_p_and_hmp(res_) -> res_
    return(res_)
  }) -> list_pval
save(list_pval, file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_circadian_pval_new.rda')

list_trios_res_scaled %>% 
  map2(.x=.,.y=names(.), .f=function(x,y){
    model_ = y
    x %>% 
      map2(.x=.,.y=names(.), .f=function(x,y){
        sign_ = y
        if (class(x) == "data.frame"){
          res_df = x
          res_df %>% mutate(model = model_, sign = sign_, level = NA) -> res_df
        }else{
          x %>% 
            map2(.x=.,.y=names(.), .f=function(x,y){
              level_ = y
              res_df = x
              res_df %>% mutate(model = model_, sign = sign_, level = level_) -> res_df
            }) %>% do.call(rbind, .) -> res_df
        }
      }) %>% do.call(rbind, .) -> res_df
  }) %>% do.call(rbind, .) -> trios_res_scaled_df

list(trios_res_df = trios_res_scaled_df) %>% 
  map(function(x){
    df_ = x
    gene.pos.df = gene.ref %>% mutate(gene.pos = sprintf("%s-%s-%s", seqnames, start, end)) %>% 
      mutate(gene_strand = strand) %>% 
      dplyr::select(gene_name, gene.pos, gene_strand)
    tf.pos.df = gene.ref %>% mutate(tf.pos = sprintf("%s-%s-%s", seqnames, start, end)) %>% 
      dplyr::select(gene_name, tf.pos)
    left_join(df_, gene.pos.df, by = c("gene" = "gene_name")) -> df_
    left_join(df_, tf.pos.df, by = c("TF" = "gene_name")) -> df_
    gene_tss = df_$gene.pos %>% gsub(".+?-(\\d+?)-.+", "\\1", .) %>% as.numeric()
    peak_start = df_$peak %>% gsub(".+?-(\\d+?)-.+", "\\1", .) %>% as.numeric()
    peak_end = df_$peak %>% gsub(".+?-(\\d+?)-(\\d+)", "\\2", .) %>% as.numeric() 
    peak_middle = (peak_start + peak_end)/2
    dist = gene_tss - peak_middle
    df_ %>% mutate(peak_to_tss_dist = dist) -> df_
  }) %>% .[[1]] -> trios_res_scaled_df

trios_res_scaled_df$Model_1 = "TRIPOD"

##Add CRE annotations (peak biotypes)##
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
gene.ref[gene.ref$gene_biotype == "lincRNA"] -> gene.ref.lincRNA
trios_res_scaled_df %>% 
  {
    df_ = .
    GenomicRanges::GRanges(
      seqnames = df_$peak %>% gsub("(.+)-.+-.+", "\\1", .),
      IRanges::IRanges(
        start = df_$peak %>% gsub("(.+)-(.+)-.+", "\\2", .) %>% as.integer(),
        end = df_$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
      )
    ) %>% unique() -> trios_gr
    ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, gene.ref.lincRNA)
  } -> ol
ol$overlappingPeaks[[1]][,c(2,3,4,14,15)] %>% 
  mutate(peak = sprintf("%s-%s-%s", seqnames, start, end)) %>% 
  dplyr::select(peak, gene_name, gene_biotype) %>% 
  "colnames<-"(., gsub("gene", "peak", colnames(.))) %>% 
  "rownames<-"(., 1:nrow(.)) -> lincRNA_map

read.table(file = "~/Dropbox/singulomics/github_rda/TRIPOD/eRNAbase/eRNA_merge_mouse_nonred.bed", header = F, sep = "\t", stringsAsFactors = F) -> df_
GenomicRanges::GRanges(
  seqnames = df_$V1,
  IRanges::IRanges(
    start = df_$V2,
    end = df_$V3
  )
) %>% unique() -> eRNA_gr

trios_res_scaled_df %>% 
  {
    df_ = .
    GenomicRanges::GRanges(
      seqnames = df_$peak %>% gsub("(.+)-.+-.+", "\\1", .),
      IRanges::IRanges(
        start = df_$peak %>% gsub("(.+)-(.+)-.+", "\\2", .) %>% as.integer(),
        end = df_$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
      )
    ) %>% unique() -> trios_gr
    ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, eRNA_gr)
  } -> ol

ol$overlappingPeaks[[1]][,c(2,3,4)] %>% 
  mutate(peak = sprintf("%s-%s-%s", seqnames, start, end)) %>% 
  mutate(peak_name = NA, peak_biotype = "eRNA") %>% 
  dplyr::select(peak, peak_name, peak_biotype) %>% 
  "rownames<-"(., 1:nrow(.)) -> eRNA_map

rbind(lincRNA_map, eRNA_map) -> map_df
map_df$peak %>% unique() -> unique_peaks

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 6)

unique_peaks %>% 
  "names<-"(.,.) %>% 
  future_map(function(x){
    peak_ = x
    map_df %>% dplyr::filter(peak == peak_) -> df_
    data.frame(peak = peak_, 
               peak_name = df_$peak_name %>% {.[!is.na(.)]} %>% unique() %>% paste(., collapse = ";"), 
               peak_biotype = df_$peak_biotype %>% unique() %>% sort() %>% paste(., collapse = ";"))
  }) -> map_df
map_df %>% do.call(rbind, .) -> map_df

full_join(x = trios_res_scaled_df, y = map_df, by = "peak") -> trios_res_scaled_df_1

trios_res_scaled_df_1 %>% 
  mutate(peak_biotype = case_when(
    is.na(peak_biotype) ~ "No annotation", 
    TRUE ~ peak_biotype
  )) -> trios_res_scaled_df_1

##Add circadian pval##
load(file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_circadian_pval_new.rda')
list_pval$rna %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("gene", "gene_cauchy_BH.Q", "gene_phase")) -> gene_map
list_pval$motif.score %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("TF", "TF_cauchy_BH.Q", "TF_phase")) -> TF_map
list_pval$peak %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("peak", "peak_cauchy_BH.Q", "peak_phase")) -> peak_map

dim(trios_res_scaled_df_1)
gene_map %>% 
  left_join(x = trios_res_scaled_df_1, y = ., by = "gene") %>% 
  left_join(x = ., y = TF_map, by = "TF") %>% 
  left_join(x = ., y = peak_map, by = "peak") -> trios_res_scaled_df_1

rm(gene_map, TF_map, peak_map)
save(trios_res_scaled_df_1, file = "~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1_new_scaled.RData")
####

#Add Chromatin state annotations##
c("15-State", "18-State") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    list.files(sprintf("~/Dropbox/singulomics/github_rda/TRIPOD/Mouse_liver_chromHMM/%s", x), pattern = "\\.bed$", full.names = T) %>% 
      map(function(x){
        print(x)
        group_ = gsub("/.+/ChromHMM_mouse_liver_(.+)_EN.+", "\\1", x)
        read.table(x, header = F, sep = "\t", stringsAsFactors = F) -> df_
        df_ = df_[,c("V1", "V2", "V3", "V4")]
        colnames(df_) = c("chr", "start", "end", group_)
        gc()
        return(df_)
      }) %>% purrr::reduce(., full_join, by = c("chr", "start", "end")) -> ChromHMM_df
  }) -> ChromHMM_df_list

trios_res_scaled_df_1 %>% 
  {
    trios_df = .
    chr_ = trios_df$peak %>% gsub("(.+)-(.+)-(.+)", "\\1", .)
    start_ = trios_df$peak %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer()
    end_ = trios_df$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
    data.frame(chr = chr_, start = start_, end = end_) %>% 
      GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% unique() -> trios_df
    trios_df
  } -> trios_gr

ChromHMM_df_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    print(group_)
    df_ = x %>% dplyr::select(chr, start, end, postnatal_0) %>% dplyr::filter(!is.na(postnatal_0))
    GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> ChromHMM_gr
    if (group_ == "15-State"){
      ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, ChromHMM_gr, minoverlap = 200) -> ol
    }else{
      ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, ChromHMM_gr) -> ol 
    }
    return(ol)
  }) -> ol_list

ol_list %>% 
  map(function(x){
    ol = x
    ol$overlappingPeaks[[1]] -> ol
    
    colnames(ol) -> col_
    col_[2:6] = col_[2:6] %>% sprintf("%s_1", .)
    col_[8:12] = col_[8:12] %>% sprintf("%s_2", .)
    colnames(ol) = col_
    
    overlapped_bp_ = function(start_1, end_1, start_2, end_2){
      sort(c(start_1, end_1, start_2, end_2)) %>% {.[3]-.[2]} -> overlapped_bp
      return(overlapped_bp)
    }
    
    ol %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(overlapped_bp = overlapped_bp_(start_1, end_1, start_2, end_2)) -> ol
  }) -> ol_list_1
ol_list = ol_list_1; rm(ol_list_1)

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 6)

ol_list %>% 
  #  .[2] %>% 
  map(function(x){
    #    ol = x[1:200, ]
    ol = x
    ol %>% mutate(peak = sprintf("%s-%s-%s", seqnames_1, start_1, end_1), .before = 1) -> ol
    ol$peak %>% unique() %>% 
      #      map(function(x){
      future_map(function(x){
        peak_ = x
        ol %>% filter(peak == peak_) -> df_
        max_overlap = df_$overlapped_bp %>% max()
        df_ %>% dplyr::filter(overlapped_bp == max_overlap) %>% .$postnatal_0 %>% unique() %>% paste(., collapse = ";") -> max_overlap_state
        data.frame(peak = peak_, chromatin_state = max_overlap_state)
      }) %>% do.call(rbind, .)
  }) -> mapping_df_list

save(ChromHMM_df_list, ol_list, mapping_df_list, trios_gr, 
     file = "~/Dropbox/singulomics/github_rda/TRIPOD/ChromHMM_dat_new_scaled.RData")

left_join(x = trios_res_scaled_df_1, y = mapping_df_list$`15-State`, by = "peak") %>% 
  "colnames<-"(., gsub("^chromatin_state$", "chromatin_state_15", colnames(.))) %>% 
  left_join(x = ., y = mapping_df_list$`18-State`, by = "peak") %>% 
  "colnames<-"(., gsub("^chromatin_state$", "chromatin_state_18", colnames(.))) -> trios_res_scaled_df_1

trios_res_scaled_df_1 %>% 
  dplyr::mutate(chromatin_state_15 = case_when(
    is.na(chromatin_state_15) ~ "no_annotation",
    TRUE ~ chromatin_state_15
  )) %>% 
  dplyr::mutate(chromatin_state_18 = case_when(
    is.na(chromatin_state_18) ~ "no_annotation",
    TRUE ~ chromatin_state_18
  )) -> trios_res_scaled_df_1

save(trios_res_scaled_df_1, file = "~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1_new_scaled.RData")

setdiff(
  colnames(metacell.peak), 
  trios_res_scaled_df_1$peak %>% unique()
) %>% 
  {
    peak = .
    chr_ = peak %>% gsub("(.+)-(.+)-(.+)", "\\1", .)
    start_ = peak %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer()
    end_ = peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
    data.frame(chr = chr_, start = start_, end = end_) %>% 
      GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = T) %>% unique() -> df_
    df_
  } -> rest_peak_gr

ChromHMM_df_list %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    print(group_)
    df_ = x %>% dplyr::select(chr, start, end, postnatal_0) %>% dplyr::filter(!is.na(postnatal_0))
    GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> ChromHMM_gr
    if (group_ == "15-State"){
      ChIPpeakAnno::findOverlapsOfPeaks(rest_peak_gr, ChromHMM_gr, minoverlap = 200) -> ol
    }else{
      ChIPpeakAnno::findOverlapsOfPeaks(rest_peak_gr, ChromHMM_gr) -> ol 
    }
    return(ol)
  }) -> ol_list_rest

ol_list_rest %>% 
  map(function(x){
    ol = x
    ol$overlappingPeaks[[1]] -> ol
    
    colnames(ol) -> col_
    col_[2:6] = col_[2:6] %>% sprintf("%s_1", .)
    col_[8:12] = col_[8:12] %>% sprintf("%s_2", .)
    colnames(ol) = col_
    
    overlapped_bp_ = function(start_1, end_1, start_2, end_2){
      sort(c(start_1, end_1, start_2, end_2)) %>% {.[3]-.[2]} -> overlapped_bp
      return(overlapped_bp)
    }
    
    ol %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(overlapped_bp = overlapped_bp_(start_1, end_1, start_2, end_2)) -> ol
  }) -> ol_list_rest

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 6)

ol_list_rest %>% 
  #  .[2] %>% 
  map(function(x){
    #    ol = x[1:200, ]
    ol = x
    ol %>% mutate(peak = sprintf("%s-%s-%s", seqnames_1, start_1, end_1), .before = 1) -> ol
    ol$peak %>% unique() %>% 
      #      map(function(x){
      future_map(function(x){
        peak_ = x
        ol %>% filter(peak == peak_) -> df_
        max_overlap = df_$overlapped_bp %>% max()
        df_ %>% dplyr::filter(overlapped_bp == max_overlap) %>% .$postnatal_0 %>% unique() %>% paste(., collapse = ";") -> max_overlap_state
        data.frame(peak = peak_, chromatin_state = max_overlap_state)
      }) %>% do.call(rbind, .)
  }) -> mapping_df_list_rest

names(mapping_df_list) %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    rbind(
      mapping_df_list[[x]], 
      mapping_df_list_rest[[x]]
    ) -> mapping_df_all
  }) -> mapping_df_list_all

save(ChromHMM_df_list, ol_list, mapping_df_list, trios_gr, mapping_df_list_all,
     file = "~/Dropbox/singulomics/github_rda/TRIPOD/ChromHMM_dat_new_scaled.RData")
####

##Plot Chromatin state validataion (Hypergeometric test)##
load("~/Dropbox/singulomics/github_rda/TRIPOD/ChromHMM_dat_new_scaled.RData")

c("chromatin_state_15", "chromatin_state_18") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    trios_res_df_1 %>% 
      dplyr::filter(Model_1 == "TRIPOD") -> df_
    df_[[x]] %>% str_split(";") %>% unlist() -> chromatin_state
    table(chromatin_state) %>% 
      {(./sum(.))*100} %>% 
      as.data.frame() -> df_
    df_$chromatin_state %>% {.[!grepl("no_annotation", .)]} %>% as.character() %>% c(., "no_annotation") -> level_
    df_ %>% mutate(chromatin_state = factor(chromatin_state, level = level_)) -> df_
    df_ %>% 
      ggplot(aes(x = chromatin_state, y = Freq, fill = chromatin_state)) + 
      geom_bar(stat = "identity") + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ylab("Percentage") + 
      xlab(NULL)
  }) -> p_list

mapping_df_list_all %>%
  map(function(x){
    df_ = x
    colnames(metacell.peak) %>% 
      as.data.frame() %>% 
      "colnames<-"(., "peak") %>% 
      #      head(500) %>% 
      left_join(x = ., y = df_, by = "peak") -> df_
    df_ %>% mutate(chromatin_state = case_when(
      is.na(chromatin_state) ~ "no_annotation",
      TRUE ~ chromatin_state
    )) -> df_
    
    df_$chromatin_state %>% 
      str_split(";") %>%
      unlist() %>% 
      table() %>% 
      {
        (./sum(.))*100
      } -> chromatin_state
    chromatin_state %>% 
      as.data.frame() %>% 
      "colnames<-"(., c("chromatin_state", "Freq")) -> df_
    df_$chromatin_state %>% {.[!grepl("no_annotation", .)]} %>% as.character() %>% c(., "no_annotation") -> level_
    df_ %>% dplyr::mutate(chromatin_state = factor(chromatin_state, level = level_)) -> df_
    df_ %>% 
      ggplot(aes(x = chromatin_state, y = Freq, fill = chromatin_state)) + 
      geom_bar(stat = "identity") + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ylab("Percentage") + 
      xlab(NULL)
  }) -> p_list_1

trios_res_scaled_df_1 %>% dplyr::filter(Model_1 == "TRIPOD") %>% .$peak %>% unique() %>% length() -> N_TRIPOD_unique_peaks
mapping_df_list_all %>%
  map(function(x){
    set.seed(123)
    df_ = x
    colnames(metacell.peak) %>% 
      sample(., N_TRIPOD_unique_peaks, replace = F) %>%
      as.data.frame() %>% 
      "colnames<-"(., "peak") %>% 
      #      head(500) %>% 
      left_join(x = ., y = df_, by = "peak") -> df_
    df_ %>% mutate(chromatin_state = case_when(
      is.na(chromatin_state) ~ "no_annotation",
      TRUE ~ chromatin_state
    )) -> df_
    
    df_$chromatin_state %>% 
      str_split(";") %>%
      unlist() %>% 
      table() %>% 
      {
        (./sum(.))*100
      } -> chromatin_state
    chromatin_state %>% 
      as.data.frame() %>% 
      "colnames<-"(., c("chromatin_state", "Freq")) -> df_
    df_$chromatin_state %>% {.[!grepl("no_annotation", .)]} %>% as.character() %>% c(., "no_annotation") -> level_
    df_ %>% dplyr::mutate(chromatin_state = factor(chromatin_state, level = level_)) -> df_
    df_ %>% 
      ggplot(aes(x = chromatin_state, y = Freq, fill = chromatin_state)) + 
      geom_bar(stat = "identity") + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      ylab("Percentage") + 
      xlab(NULL)
  }) -> p_list_2

ggplot() + 
  geom_bar(data = p_list$chromatin_state_18$data, aes(x = chromatin_state, y = Freq, fill = chromatin_state), stat = "identity") + 
  geom_bar(data = p_list_1$`18-State`$data, aes(x = chromatin_state, y = Freq), stat = "identity", width = 0.5) + 
  geom_bar(data = p_list_2$`18-State`$data, aes(x = chromatin_state, y = Freq), stat = "identity", width = 0.2, fill = "grey") + 
  theme_classic() + 
  ylab("Percentage") + 
  xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p18_state

list(
  TRIPOD = p18_state$layers[[1]]$data,
  All_peaks = p18_state$layers[[2]]$data,
  Random_peaks = p18_state$layers[[3]]$data
) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    df_ = x
    colnames(df_) = c("choromatine_state", "Percentage")
    df_[["group"]] <- y
    #    head(df_)
    df_
  }) %>% do.call(rbind, .) %>% 
  dplyr::arrange(choromatine_state) %>% 
  "rownames<-"(., 1:nrow(.)) -> p18_state_df

chromatin = p18_state_df

chromatin=cbind(chromatin, newstate=NA)

chromatin$newstate[grep('Tss',chromatin$choromatine_state)]='TSS'
chromatin$newstate[grep('Tx',chromatin$choromatine_state)]='Transcription'
chromatin$newstate[grep('Enh',chromatin$choromatine_state)]='Enhancer'
chromatin$newstate[grep('Quies',chromatin$choromatine_state)]='Quiescent'
chromatin$newstate[grep('Repr',chromatin$choromatine_state)]='Repressive'
chromatin$newstate[grep('Het',chromatin$choromatine_state)]='Repressive'

chromatin=chromatin[-grep('no_annotation',chromatin$choromatine_state),]

# Need to rescale the prop to sum up to 1.
all_peaks=chromatin[chromatin$group=='All_peaks',]
all_peaks$Percentage=all_peaks$Percentage/sum(all_peaks$Percentage)
all_peaks=aggregate(Percentage~newstate+group, all_peaks, sum)

tripod=chromatin[chromatin$group=='TRIPOD',]
tripod$Percentage=tripod$Percentage/sum(tripod$Percentage)
tripod=aggregate(Percentage~newstate+group, tripod, sum)

to_plot=rbind(all_peaks, tripod)
to_plot$newstate=factor(to_plot$newstate, levels= c('Transcription', 'TSS',
                                                    'Enhancer', 'Repressive',
                                                    'Quiescent'))

library(ggplot2)
p=ggplot(data=to_plot, aes(x=newstate, y=Percentage, fill=group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() + xlab("chromHMM states")+ scale_fill_manual(values=c('#999999','#E69F00')) #Fig_5F

ggsave(plot = p, filename = "~/Dropbox/singulomics/github_rda/TRIPOD/chromatin_state_validation_new_scaled.pdf") #Fig_5F

## Hi-C validation ##
trios_res_df_2 = trios_res_scaled_df_1
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/GSE155161", pattern = "washU", full.names = T) %>% 
  gtools::mixedsort() %>% 
  as.data.frame() %>% 
  "colnames<-"(., c("file")) %>% 
  mutate(ZT = gsub("/.+/.+(ZT\\d+)All.+", "\\1", file)) %>% 
  {
    df_ = .
    unique(df_$ZT) %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(x){
        ZT_ = x
        df_ %>% dplyr::filter(ZT == ZT_) %>% .$file -> file_
        print(file_)
        read.table(file_, sep = "\t", header = F, stringsAsFactors = F) -> df_
        score_ = df_$V3
        
        peak_df_1 = data.frame(V1 = gsub("(.+),.+,.+", "\\1", df_$V1), 
                               V2 = gsub("(.+),(.+),.+", "\\2", df_$V1) %>% as.integer(),
                               V3 = gsub("(.+),(.+),(.+)", "\\3", df_$V1) %>% as.integer()
                               )
        peak_df_1 %>% mutate(V4 = ZT_, V5 = 1:nrow(.)) -> peak_df_1
        write.table(peak_df_1, "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed", quote = F, row.names = F, col.names = F, sep = "\t")
        input_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed" 
        liftOver_ = "~/Documents/UCSC_tools/liftOver"
        chain_ = "~/Dropbox/singulomics/github_rda/TRIPOD/mm9ToMm10.over.chain"
        output_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm10.bed"
        unlifted_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_unlifted.bed"
        command_ = sprintf("%s %s %s %s %s", liftOver_, input_, chain_, output_, unlifted_)
        system(command_)
        read.table(output_, header = F, sep = "\t", stringsAsFactors = F) -> peak_df_1
        unlink(c(input_, output_, unlifted_))
        
        peak_df_2 = data.frame(V1 = gsub("(.+),.+,.+", "\\1", df_$V2), 
                               V2 = gsub("(.+),(.+),.+", "\\2", df_$V2) %>% as.integer(),
                               V3 = gsub("(.+),(.+),(.+)", "\\3", df_$V2) %>% as.integer()
        )
        peak_df_2 %>% mutate(V4 = ZT_, V5 = 1:nrow(.)) -> peak_df_2
        write.table(peak_df_2, "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed", quote = F, row.names = F, col.names = F, sep = "\t")
        input_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed" 
        liftOver_ = "~/Documents/UCSC_tools/liftOver"
        chain_ = "~/Dropbox/singulomics/github_rda/TRIPOD/mm9ToMm10.over.chain"
        output_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm10.bed"
        unlifted_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_unlifted.bed"
        command_ = sprintf("%s %s %s %s %s", liftOver_, input_, chain_, output_, unlifted_)
        system(command_)
        read.table(output_, header = F, sep = "\t", stringsAsFactors = F) -> peak_df_2
        unlink(c(input_, output_, unlifted_))
        
#        rbind(peak_df_1[,1:3], peak_df_2[,1:3]) %>% distinct() -> peak_df_map
#        list(peak_pair_1 = peak_df_1, peak_pair_2 = peak_df_2)
        rbind(peak_df_1, peak_df_2) -> peak_df
        peak_df$V6 = score_[peak_df$V5]
        return(peak_df)
      }) %>% do.call(rbind, .)
  } -> list_GSE155161_Promoter_Capture_HI_C_df
list_GSE155161_Promoter_Capture_HI_C_df[, c(1,2,3)] %>% 
  distinct() %>% nrow() %>% {./2} -> total_peak_pair

list(df_ = list_GSE155161_Promoter_Capture_HI_C_df) %>% 
  map(function(x){
    df_ = x
    GenomicRanges::GRanges(seqnames = df_$V1, IRanges::IRanges(start = df_$V2, end = df_$V3)) -> gr_
    gr_ = unique(gr_)
    ChIPpeakAnno::findOverlapsOfPeaks(gr_, gene.tss.ref.unique, maxgap = 0) -> ol
    ol$overlappingPeaks[[1]] -> df_
    df_[,c(2,3,4,13,14)] -> df_
  }) %>% .[[1]] -> mapping_df

left_join(x = list_GSE155161_Promoter_Capture_HI_C_df, y = mapping_df, by = c("V1"="seqnames", "V2"="start", "V3" = "end")) -> list_GSE155161_Promoter_Capture_HI_C_df_1
dim(list_GSE155161_Promoter_Capture_HI_C_df_1)
list_GSE155161_Promoter_Capture_HI_C_df_1 %>% 
  mutate(peak = paste(V1, V2, V3, sep = ":")) %>%
  "colnames<-"(., c("chr", "start", "end", "ZT", "idx", "score", "gene_name", "gene_biotype", "peak")) %>% 
#  dplyr::filter(ZT == "ZT0", idx %in% 1:100) %>% 
  group_by(ZT, idx, score) %>% 
  group_map(function(x,y){
#    print(y)
    df_ = x
    if (all(is.na(df_$gene_name))|all(!is.na(df_$gene_name))){
#      print(df_$gene_name)
      data.frame(gene = character(), biotype = character(), CRE = character(), gene_peak = character(), ZT = character(), idx = numeric(), score = numeric()) -> df_out
    }else{
      gene_ = df_$gene_name %>% {.[!is.na(.)]}
#      gene_ = df_$gene_name %>% {.[!is.na(.)]} %>% paste(., collapse = ";")
#      print(sprintf("gene: %s", gene_))
      biotype_ = df_$gene_biotype %>% {.[!is.na(.)]}
#      biotype_ = df_$gene_biotype %>% {.[!is.na(.)]} %>% paste(., collapse = ";") 
#      print(sprintf("biotype: %s", biotype_)) 
      CRE_ = df_ %>% dplyr::filter(is.na(gene_name)) %>% .$peak %>% unique() %>% paste(., collapse = ",")
      gene_peak_ = df_ %>% dplyr::filter(!is.na(gene_name)) %>% .$peak %>% unique() %>% paste(., collapse = ",")
#      print(sprintf("CRE: %s", CRE_)) 
      data.frame(gene = gene_, biotype = biotype_, CRE = CRE_, gene_peak = gene_peak_, ZT = y$ZT, idx = y$idx, score = y$score) -> df_out
    }
    return(df_out)
  }, .keep = T) %>% do.call(rbind, .) -> GSE155161_annotated_gene_promoter_capture
GSE155161_annotated_CHI_C = GSE155161_annotated_gene_promoter_capture 
rm(GSE155161_annotated_gene_promoter_capture)
save(GSE155161_annotated_CHI_C, file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE155161_CHI-C_annotated.rda")

load(file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE155161_CHI-C_annotated.rda")
trios_res_df_2 = trios_res_scaled_df_1
intersect(
  trios_res_df_2 %>% filter(Model_1 == "TRIPOD") %>% .$gene %>% unique(),
  GSE155161_annotated_CHI_C$gene %>% unique()
) -> intersected_genes

options(bedtools.path = "/Users/chunyiptong/anaconda3/envs/bedtools_env/bin")
i = 0
c(
  trios_res_df_2 %>% filter(Model_1 == "TRIPOD") %>% .$gene %>% unique(),
  GSE155161_annotated_CHI_C$gene %>% unique()
) %>% unique() %>% 
  "names<-"(., .) %>% 
  #  .[1] %>% 
  map(function(x){
    gene_ = x
    i <<- i+1
    print(i)
    if (gene_ %in% (trios_res_df_2 %>% dplyr::filter(Model_1 == "TRIPOD") %>% .$gene)){
      trios_res_df_2 %>% dplyr::filter(Model_1 == "TRIPOD" & gene == gene_) %>% dplyr::select(gene, peak) -> trios_res_df_2
      trios_res_df_2 = data.frame(
        chr = trios_res_df_2$peak %>% gsub("(.+)-(.+)-(.+)", "\\1", .), 
        start = trios_res_df_2$peak %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer(), 
        end = trios_res_df_2$peak %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer() 
      )      
    }else{
      data.frame(chr = character(), start = integer(), end = integer()) -> trios_res_df_2
    }
    
    if (gene_ %in% GSE155161_annotated_CHI_C$gene){
      GSE155161_annotated_CHI_C %>% dplyr::filter(gene == gene_) -> GSE155161_annotated_CHI_C
      GSE155161_annotated_CHI_C = data.frame(
        chr = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\1", .), 
        start = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\2", .) %>% as.integer(), 
        end = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\3", .) %>% as.integer() 
      )      
    }else{
      data.frame(chr = character(), start = integer(), end = integer()) -> GSE155161_annotated_CHI_C
    }
    
    rbind(trios_res_df_2, GSE155161_annotated_CHI_C) -> df_final
    if (nrow(df_final) >= 1){
      df_final %>% bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% 
        mutate(V4 = gene_) -> df_final      
    }else{
      df_final = data.frame(chr = character(), start = integer(), end = integer(), V4 = character())
    }
    return(df_final)
  }) %>% do.call(rbind, .) -> all_gene_peaks_linkage

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 6)

intersected_genes %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  #  .["Tmem14a"] %>% 
  future_map(function(x){
    gene_ = x
    print(gene_)
    
    trios_res_df_2 %>% filter(Model_1 == "TRIPOD" & gene == gene_) -> df_trios
    GenomicRanges::GRanges(gene = sprintf("peak_%s_%s", gene_, 1:nrow(df_trios)), 
                           seqnames = df_trios$peak %>% gsub("(.+)-.+-.+", "\\1", .), 
                           IRanges::IRanges(start = df_trios$peak %>% gsub("(.+)-(.+)-.+", "\\2", .) %>% as.integer(), 
                                            end = df_trios$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
                           ) 
    ) -> df_trios
    df_trios = unique(df_trios)
    
    GSE155161_annotated_CHI_C %>% filter(gene == gene_) -> df_CHI_C
    GenomicRanges::GRanges(gene = sprintf("CRE_%s_%s", gene_, 1:nrow(df_CHI_C)), 
                           seqnames = df_CHI_C$CRE %>% gsub("(.+):.+:.+", "\\1", .), 
                           IRanges::IRanges(start = df_CHI_C$CRE %>% gsub("(.+):(.+):.+", "\\2", .) %>% as.integer(), 
                                            end = df_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\3", .) %>% as.integer()
                           ) 
    ) -> df_CHI_C
    df_CHI_C = unique(df_CHI_C)
    
    ChIPpeakAnno::findOverlapsOfPeaks(df_trios, df_CHI_C, maxgap = 1) -> ol
    ol$venn_cnt -> venn_cnt
    venn_cnt[3, "count.df_trios"] -> trios_peak_unique
    venn_cnt[2, "count.df_CHI_C"] -> CHI_C_peak_unique
    venn_cnt[4, "count.df_trios"] -> trios_peak_intersect
    venn_cnt[4, "count.df_CHI_C"] -> CHI_C_peak_intersect
    data.frame(gene = gene_, trios_peak_unique = trios_peak_unique,
               CHI_C_peak_unique = CHI_C_peak_unique, trios_peak_intersect = trios_peak_intersect,
               CHI_C_peak_intersect = CHI_C_peak_intersect) -> df_final
    intersected_peaks = max(trios_peak_intersect, CHI_C_peak_intersect)
    A = intersected_peaks+trios_peak_unique
    B = intersected_peaks+CHI_C_peak_unique
    N = nrow(all_gene_peaks_linkage)
    k = intersected_peaks
    p_value <- phyper(k - 1, A, N - A, B, lower.tail = FALSE)
    df_final %>% mutate(p_value = p_value) -> df_final
  }) -> gene_peak_linkage_hypergeometric_test

gene_peak_linkage_hypergeometric_test %>% do.call(rbind, .) -> gene_peak_linkage_hypergeometric_test_df

save(all_gene_peaks_linkage, gene_peak_linkage_hypergeometric_test_df, file = "~/Dropbox/singulomics/github_rda/TRIPOD/hypergeometic_test_dat_new_scaled.rda")

trios_res_df_2 %>% 
  {
    df_ = .
    set.seed(123)
    sample(1:nrow(df_), nrow(df_), replace = F) -> se
    df_[["peak"]] = df_[["peak"]][se]
    df_
  } -> trios_res_df_2_shuffled

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 6)

intersected_genes %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  #  .["Tmem14a"] %>% 
  future_map(function(x){
    gene_ = x
    print(gene_)
    
    trios_res_df_2_shuffled %>% filter(Model_1 == "TRIPOD" & gene == gene_) -> df_trios
    GenomicRanges::GRanges(gene = sprintf("peak_%s_%s", gene_, 1:nrow(df_trios)), 
                           seqnames = df_trios$peak %>% gsub("(.+)-.+-.+", "\\1", .), 
                           IRanges::IRanges(start = df_trios$peak %>% gsub("(.+)-(.+)-.+", "\\2", .) %>% as.integer(), 
                                            end = df_trios$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
                           ) 
    ) -> df_trios
    df_trios = unique(df_trios)
    
    GSE155161_annotated_CHI_C %>% filter(gene == gene_) -> df_CHI_C
    GenomicRanges::GRanges(gene = sprintf("CRE_%s_%s", gene_, 1:nrow(df_CHI_C)), 
                           seqnames = df_CHI_C$CRE %>% gsub("(.+):.+:.+", "\\1", .), 
                           IRanges::IRanges(start = df_CHI_C$CRE %>% gsub("(.+):(.+):.+", "\\2", .) %>% as.integer(), 
                                            end = df_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\3", .) %>% as.integer()
                           ) 
    ) -> df_CHI_C
    df_CHI_C = unique(df_CHI_C)
    
    ChIPpeakAnno::findOverlapsOfPeaks(df_trios, df_CHI_C, maxgap = 1) -> ol
    ol$venn_cnt -> venn_cnt
    venn_cnt[3, "count.df_trios"] -> trios_peak_unique
    venn_cnt[2, "count.df_CHI_C"] -> CHI_C_peak_unique
    venn_cnt[4, "count.df_trios"] -> trios_peak_intersect
    venn_cnt[4, "count.df_CHI_C"] -> CHI_C_peak_intersect
    data.frame(gene = gene_, trios_peak_unique = trios_peak_unique,
               CHI_C_peak_unique = CHI_C_peak_unique, trios_peak_intersect = trios_peak_intersect,
               CHI_C_peak_intersect = CHI_C_peak_intersect) -> df_final
    intersected_peaks = max(trios_peak_intersect, CHI_C_peak_intersect)
    A = intersected_peaks+trios_peak_unique
    B = intersected_peaks+CHI_C_peak_unique
    N = nrow(all_gene_peaks_linkage)
    k = intersected_peaks
    p_value <- phyper(k - 1, A, N - A, B, lower.tail = FALSE)
    df_final %>% mutate(p_value = p_value) -> df_final
  }) -> gene_peak_linkage_hypergeometric_test_shuffled
gene_peak_linkage_hypergeometric_test_shuffled %>% do.call(rbind, .) -> gene_peak_linkage_hypergeometric_test_df_shuffled

gene_peak_linkage_hypergeometric_test_df %>% 
  mutate(p_value = case_when(
    p_value == 0 ~ .Machine$double.eps, 
    TRUE ~ p_value
  )) %>% mutate(group = "TRIPOD vs CHI-C") %>% 
  dplyr::filter(p_value < 0.05) %>% 
  mutate(p_value = -log(p_value, 10)) %>% 
  ggplot(aes(x = p_value)) + 
  geom_histogram(binwidth = 0.5, color = "black") + 
  xlab("Hypergeometric test -log10(p-value)") + 
  ggtitle("4241 sig. genes") + 
  theme_classic() -> p

gene_peak_linkage_hypergeometric_test_df %>% 
  mutate(p_value = case_when(
    p_value == 0 ~ .Machine$double.eps, 
    TRUE ~ p_value
  )) %>% mutate(group = "TRIPOD vs CHI-C") %>% 
  dplyr::filter(p_value < 0.05) %>% 
  {
    gene_ = .$gene
    gene_peak_linkage_hypergeometric_test_df_shuffled[gene_, ] %>% 
      mutate(p_value = case_when(
        p_value == 0 ~ .Machine$double.eps, 
        TRUE ~ p_value
      )) %>% mutate(group = "TRIPOD vs CHI-C") %>% 
      mutate(p_value = -log(p_value, 10)) %>% 
      ggplot(aes(x = p_value)) + 
      geom_histogram(binwidth = 0.5, color = "black") + 
      xlab("Hypergeometric test -log10(p-value)") + 
      ggtitle("4241 sig. genes") + 
      theme_classic()
  } -> p1

list(
  TRIPOD = p$data,
  random = p1$data
) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    df_ = x
    df_ %>% dplyr::select(gene, p_value) %>% 
      "colnames<-"(., c("gene", "-log10_P")) -> df_
    df_[["group"]] <- y
    df_
  }) %>% do.call(rbind, .) %>% 
  dplyr::arrange(gene) %>% 
  "rownames<-"(., 1:nrow(.)) -> df_hypergeometric_test

pvals = df_hypergeometric_test
tripod.p=10^(-pvals$`-log10_P`[pvals$group=='TRIPOD'])
random.p=10^(-pvals$`-log10_P`[pvals$group=='random'])

tripod.p=pmin(70, (pvals$`-log10_P`[pvals$group=='TRIPOD']))
random.p=(pvals$`-log10_P`[pvals$group=='random'])

pdf(file='~/Dropbox/singulomics/github_rda/HI_C_validation_new.pdf', width=5, height=4) #Fig_5F
hist(tripod.p, breaks=seq(0,70,2), col=adjustcolor("#E69F00", alpha.f=0.7), xlab='-log(p)', main='Distribution of -log(p) \nfrom hypergeometric test')
hist(random.p, add=TRUE, breaks=seq(0,70,2), col=adjustcolor("#999999", alpha.f = 0.7))
legend('topright', fill=c(adjustcolor("#999999", alpha.f=0.7), adjustcolor("#E69F00", alpha.f = 0.7)),
       legend = c('Random peaks', 'TRIPOD peaks'), bty='n')
dev.off()
####

## Process ChIP-seq data ##
##Process ChIP peak from GSE10115 (Rora and Rorg)
library(bedtoolsr)
options(bedtools.path = "/Users/chunyiptong/anaconda3/envs/bedtools_env/bin")
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/GSE101115_Rora_Rorg_ChIP", pattern = "narrowPeak|broadPeak", full.names = T) %>% 
  {.[grepl("narrow", .)]} %>% 
  {
  peak_files = .
  c("Rora", "Rorg") %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      gene_ = x
      peak_files %>% {.[grepl(sprintf("/%s", gene_), .)]} -> peak_files
      c("unique", "multiple") %>% 
        "names<-"(.,.) %>% 
        map(function(x){
          unique_ = x
          peak_files %>% {.[grepl(unique_, .)]} -> peak_files_
          peak_files_ %>% 
            map(function(x){
              print(x)
              read.table(x, header = F, sep = "\t", stringsAsFactors = F) %>% as_tibble() %>% .[,1:9] -> df_
              colnames(df_) = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
              print(min(df_$qValue))
              
              threshold_ = -log(0.1, 10)
              df_ %>% dplyr::filter(qValue > threshold_) -> df_
            }) %>% do.call(rbind, .) -> df_
          df_ %>% dplyr::select(chr, start, end) -> df_
          df_ %>% bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble() -> df_
        }) %>% do.call(rbind, .) %>% bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble() -> df_
    })
} -> GSE101115_peak_list
save(GSE101115_peak_list, file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE101115_Rora_Rorg_ChIP/GSE101115_peak_list.rda")
####

##Process ChIP peak from GSE26345 (Nr1d1)
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/GSE26345_Nr1d1_ChIP", pattern = "narrowPeak|broadPeak", full.names = T) %>% 
  {.[grepl("narrow", .)]} %>% 
  map(function(x){
    print(x)
    read.table(x, header = F, sep = "\t", stringsAsFactors = F) %>% as_tibble() %>% .[,1:9] -> df_
    colnames(df_) = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
    print(min(df_$qValue))
    
    threshold_ = -log(0.1, 10)
    df_ %>% dplyr::filter(qValue > threshold_) -> df_ 
  }) %>% do.call(rbind, .) %>% 
  dplyr::select(chr, start, end) %>% 
  bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble() -> GSE26345_peak_list
save(GSE26345_peak_list, file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE26345_Nr1d1_ChIP/GSE26345_peak_list.rda")
####

##Process ChIP peak from GSE59486 (Nfil3 and Rora)
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/GSE59486_Nfil3_Rora_ChIP", pattern = "narrowPeak|broadPeak", full.names = T) %>% 
  {.[grepl("narrow", .)]} %>% 
  {
    peak_files = .
    c("E4BP4", "ROR_alpha") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        gene_ = x
        peak_files %>% {.[grepl(sprintf("/%s", gene_), .)]} -> peak_files
        peak_files %>% 
          map(function(x){
            print(x)
            read.table(x, header = F, sep = "\t", stringsAsFactors = F) %>% as_tibble() %>% .[,1:9] -> df_
            colnames(df_) = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
            print(min(df_$qValue))
            
            threshold_ = -log(0.1, 10)
            df_ %>% dplyr::filter(qValue > threshold_) -> df_  
          }) %>% do.call(rbind, .) %>%
          dplyr::select(chr, start, end) %>% 
          bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble()
      })
  } -> GSE59486_peak_list
save(GSE59486_peak_list, file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE59486_Nfil3_Rora_ChIP/GSE59486_peak_list.rda")
####

##Process ChIP peak from PRJDB7796 (DBP and E4BP4) 
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/PRJDB7796_DBP_E4BP4_ChIP", pattern = "narrowPeak|broadPeak", full.names = T) %>% 
  {.[grepl("narrow", .)]} %>% 
  {
    peak_files = .
    c("DBP", "E4BP4") %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        gene_ = x
        peak_files %>% {.[grepl(sprintf("ZT\\d+_%s", gene_), .)]} -> peak_files
        peak_files %>% 
          map(function(x){
            print(x)
            read.table(x, header = F, sep = "\t", stringsAsFactors = F) %>% as_tibble() %>% .[,1:9] -> df_
            colnames(df_) = c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")
            print(min(df_$qValue))
            
            threshold_ = -log(0.1, 10)
            df_ %>% dplyr::filter(qValue > threshold_) -> df_  
          }) %>% do.call(rbind, .) %>%
          dplyr::select(chr, start, end) %>% 
          bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble() -> df_
      })
  } -> PRJDB7796_peak_list
save(PRJDB7796_peak_list, file = "~/Dropbox/singulomics/github_rda/TRIPOD/PRJDB7796_DBP_E4BP4_ChIP/PRJDB7796_peak_list.rda")
####

#Process ChIP peak from GSE39977 (Arntl, Npas2, Clock) 

read.table("~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/mapping.txt", header = F, sep = "\t", stringsAsFactors = F) %>% 
  "colnames<-"(., c("GSM", "sample")) -> mapping_df

 ####macs2 peak calling
mapping_df %>% mutate(group = gsub("Sample \\d+_(.+)", "\\1", sample)) %>% 
  filter(grepl("BMAL1|CLOCK|NPAS2", group)) %>% 
  mutate(group = gsub("Input DNA for ", "", group)) %>% 
  {
    df_ = .
    group_ = unique(df_$group) 
    group_ %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        group_ = x
        df_ %>% dplyr::filter(group == group_) -> df_
        df_treatment = df_ %>% filter(!grepl("Input", sample))
        if (nrow(df_treatment) > 6){
          df_treatment[1:(nrow(df_treatment)-1), ] -> df_treatment 
        }else{
          df_treatment = df_treatment
        }
        print(df_treatment)
        df_input = df_ %>% filter(grepl("Input", sample))
        if (nrow(df_input) > 1){
          input_GSM = c(rep(df_input$GSM[1], nrow(df_treatment)/2), rep(df_input$GSM[2], nrow(df_treatment)/2))
        }else{
          input_GSM = df_input$GSM
        }
        df_treatment %>% mutate(input = input_GSM) -> df_
        files = list.files("~/Downloads/GSE39977", pattern = "\\.bed$", full.names = T)
        
        df_$GSM %>% 
          map(function(x){
            files %>% {.[grepl(x, .)]} 
          }) %>% purrr::reduce(., c) -> treatment_bed
        df_$input %>% 
          map(function(x){
            files %>% {.[grepl(x, .)]} 
          }) %>% purrr::reduce(., c) -> input_bed
        CT_ = sprintf("CT_%s", seq(0, 4*(nrow(df_treatment)-1), 4))
        names_ = sprintf("%s_%s", group_, CT_) %>% gsub(" ", "_", .)
        output_dir = sprintf("~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/macs2_output/")
        command_ = sprintf("/Users/chunyiptong/anaconda3/envs/macs2_env/bin/macs2 callpeak --broad -t %s -c %s -f BED -g mm -n %s --outdir %s -p 0.05", 
                           treatment_bed, input_bed, names_, output_dir)
      }) %>% purrr::reduce(., c)
    } -> macs2_command_
for (x in macs2_command_){
  system(x)
}

list_GSE39977_peak_df = list()
list.files("~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/macs2_output", pattern = "\\.narrowPeak|\\.broadPeak",
           full.names = T) %>% 
  .[grepl("narrow", .)] %>% 
  gtools::mixedsort() %>% 
  as.data.frame() %>% "colnames<-"(., c("files")) %>% 
  mutate(group = gsub("/.+/(.+)_CT.+", "\\1", files)) %>% 
  {
    df_ = .
    group_ = unique(df_$group)
    group_ %>% 
#      .[1] %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        group_ = x
        df_ %>% filter(group == group_) -> df_
        df_$files %>% 
#          .[1] %>% 
          map(function(x){
            file_ = x
            print(file_)
            ZT_ = gsub(".+(CT_\\d+)_.+", "\\1", file_)
            print(ZT_)
            read.table(file_, header = F, sep = "\t", stringsAsFactors = F) %>% 
              "colnames<-"(., c("Chr", "Start", "End", "Peak_ID", "Score", "Strand", "SignalValue", "PValue", "QValue", "Peak")) -> df_peak
            
            c(0.05, 0.01, 0.001, 0.0001, 0.00001) %>% 
              "names<-"(., sprintf("p_%s", .)) %>% 
#              .[1] %>% 
              map2(.x=., .y = names(.), .f=function(x,y){
                df_peak
                df_peak %>% dplyr::filter(PValue > -log(x, 10)) -> df_peak
                list_GSE39977_peak_df[[group_]][[y]][[ZT_]] <<- df_peak
                data.frame(file = file_, group = group_) -> df_tmp
                df_tmp[[y]] = nrow(df_peak)
                return(df_tmp)
              }) -> list_p_threshold 
            c(0.1, 0.05, 0.01) %>% 
              "names<-"(., sprintf("q_%s", .)) %>% 
              #              .[1] %>% 
              map2(.x=., .y = names(.), .f=function(x,y){
                df_peak
                df_peak %>% dplyr::filter(QValue > -log(x, 10)) -> df_peak
                list_GSE39977_peak_df[[group_]][[y]][[ZT_]] <<- df_peak
                data.frame(file = file_, group = group_) -> df_tmp
                df_tmp[[y]] = nrow(df_peak)
                return(df_tmp)
              }) -> list_q_threshold 
            c(list_p_threshold, list_q_threshold) %>% 
            purrr::reduce(., left_join, by = c("file", "group"))
          }) %>% purrr::reduce(., rbind)
      }) %>% purrr::reduce(., rbind)
  } -> df_GSE39977_peak_summary

list_GSE39977_peak_df %>%
#  .[1] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    print(group_)
    selected_threshold = c("p_1e-04", "p_1e-05", "q_0.1", "q_0.05")
    x[selected_threshold] %>% 
#      .[1] %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        threshold_ = y
        print(threshold_)
        x %>% 
          map(function(x){
            df_ = x[,c(1:3)]
          }) %>% do.call(rbind, .) -> df_
        df_ %>% distinct() %>% bedtoolsr::bt.sort(.) %>% bedtoolsr::bt.merge(.) -> df_
      })
  }) -> list_GSE39977_peak_df_mm9

list_GSE39977_peak_df_mm9 %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    print(group_)
    x %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        threshold_ = y
        print(threshold_)
        df_ = x
        write.table(df_, "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed", quote = F, row.names = F, col.names = F, sep = "\t")
        input_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm9.bed" 
        liftOver_ = "~/Documents/UCSC_tools/liftOver"
        chain_ = "~/Dropbox/singulomics/github_rda/TRIPOD/mm9ToMm10.over.chain"
        output_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_mm10.bed"
        unlifted_ = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/tmp_unlifted.bed"
        command_ = sprintf("%s %s %s %s %s", liftOver_, input_, chain_, output_, unlifted_)
        system(command_)
        read.table(output_, header = F, sep = "\t", stringsAsFactors = F) -> df_
        #remove tmp files from the directory
        unlink(c(input_, output_, unlifted_))
        return(df_)
      })
  }) -> list_GSE39977_peak_df_mm10

c("p_1e-04", "p_1e-05", "q_0.1", "q_0.05") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    threshold_ = x
    print(threshold_)
    rbind(
      list_GSE39977_peak_df_mm10$`BMAL1_ChIP-seq`[[threshold_]], 
      list_GSE39977_peak_df_mm10$`BMAL1_ChIP-seq-2cycle`[[threshold_]] 
    ) -> df_
    df_ %>% distinct() %>% bedtoolsr::bt.sort(.) %>% bedtoolsr::bt.merge(.) -> df_
    list_GSE39977_peak_df_mm10[["BMAL1_All"]][[threshold_]] <<- df_
  })

save(list_GSE39977_peak_df, list_GSE39977_peak_df_mm9, list_GSE39977_peak_df_mm10, file = "~/Dropbox/singulomics/github_rda/TRIPOD/list_GSE39977_peak.RData")
load("~/Dropbox/singulomics/github_rda/TRIPOD/list_GSE39977_peak.RData")
rm(list_GSE39977_peak_df, list_GSE39977_peak_df_mm9)
####

## Plot ChIP-seq validation ##
load("~/Dropbox/singulomics/github_rda/TRIPOD/list_GSE39977_peak.RData")
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')

colnames(metacell.peak) %>% 
  {
    peak_ = .
    chr = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\1", peak_)
    start = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\2", peak_) %>% as.integer()
    end = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\3", peak_) %>% as.integer()
    GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end))
  } -> df_all_peaks

list_GSE39977_peak_df_mm10 %>% 
  #  .[1] %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    model_ = y
    TF_ = gsub("(^.+?)_.+", "\\1", model_)
    paste0(toupper(substr(TF_, 1, 1)), tolower(substr(TF_, 2, nchar(TF_)))) -> TF_
    if (TF_ == "Bmal1"){
      TF_ = "Arntl"
    }else{
      TF_ = TF_
    }
    x %>% 
      #      .[1] %>% 
      map2(.x=.,.y=names(.),.f=function(x,y){
        threshold_ = y
        df_ = x
        df_ %>% distinct() -> df_
        GenomicRanges::GRanges(seqnames = df_[,1], ranges = IRanges::IRanges(start = df_[,2], end = df_[,3])) -> df_ChIP
        
        trios_res_df_2 %>% dplyr::filter(Model_1 == "TRIPOD" & TF == TF_) %>% 
          .$peak %>% unique() %>% 
          {
            peak_ = .
            print(peak_[1:3])
            chr = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\1", peak_)
            start = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\2", peak_) %>% as.integer()
            end = gsub("(chr.+?)-(\\d+?)-(\\d+)", "\\3", peak_) %>% as.integer()
            #            data.frame(chr, start, end)[1:3, ]
            GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = start, end = end))
          } -> df_Trios
        
        set.seed(123)
        df_all_peaks[sample(1:length(df_all_peaks), nrow(df_Trios %>% as.data.frame()), replace = F), ] -> df_all_peaks_sampled
        
        c(0, 500, 1000, 1500, 2000) %>% 
          "names<-"(., sprintf("gap_%s", .)) %>% 
          #          .[1] %>% 
          map(function(x){
            gap_ = x
            ChIPpeakAnno::findOverlapsOfPeaks(df_ChIP, df_all_peaks, maxgap = gap_) -> ol_all_peaks
            ChIPpeakAnno::findOverlapsOfPeaks(df_ChIP, df_all_peaks_sampled, maxgap = gap_) -> ol_all_peaks_sampled
            ChIPpeakAnno::findOverlapsOfPeaks(df_ChIP, df_Trios, maxgap = gap_) -> ol
            
            ol_all_peaks$venn_cnt -> venn_cnt_all_peaks
            ol_all_peaks_sampled$venn_cnt -> venn_cnt_all_peaks_sampled
            ol$venn_cnt -> venn_cnt
            total_ChIP_peak = venn_cnt[,4] %>% sum()
            total_Trios_peak = venn_cnt[,5] %>% sum()
            Trios_unique = venn_cnt[2,5]
            Trios_intersect = venn_cnt[4,5]
            total_all_peaks = venn_cnt_all_peaks[,5] %>% sum()
            all_peaks_intersect = venn_cnt_all_peaks[4,5]
            total_sampled_peaks = venn_cnt_all_peaks_sampled[,5] %>% sum()
            sampled_peaks_intersect = venn_cnt_all_peaks_sampled[4,5]
            data.frame(
              model = model_,
              TF = TF_,
              macs2_threshold = threshold_,
              max_gap = gap_,
              total_ChIP_peak = total_ChIP_peak,
              total_Trios_peak = total_Trios_peak,
              Trios_unique = Trios_unique,
              Trios_intersect = Trios_intersect, 
              total_all_peaks = total_all_peaks, 
              all_peaks_intersect = all_peaks_intersect,
              total_sampled_peaks = total_sampled_peaks,
              sampled_peaks_intersect = sampled_peaks_intersect
            )
          }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) -> venn_list_df_

venn_list_df_ %>% filter(macs2_threshold == "q_0.1" & max_gap == 0) %>% 
  filter(!grepl("cycle|BMAL1_ChIP", model)) %>% 
  mutate(
    Trios = Trios_intersect / total_Trios_peak,
    All_peaks = all_peaks_intersect / total_all_peaks,
    Sampled_peaks = sampled_peaks_intersect / total_sampled_peaks
  ) %>% 
  dplyr::select(TF, Trios, All_peaks, Sampled_peaks) %>%
  pivot_longer(cols = -TF, names_to = "group", values_to = "ratio") %>% 
  ggplot(aes(x = TF, y = ratio, fill = group)) + 
  #  geom_bar(stat = "identity", width = 0.7) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  theme_classic() + 
  theme(legend.position = "top") +
  ylab("Fraction of peaks overlapping with ChIP-seq peaks") #Fig_5F

  ## Plot CRE phase vs Target gene phase (ARNTL and NR1D1) ##
  trios_res_df_2 %>% 
  dplyr::filter(Model_1 == "TRIPOD", gene_cauchy_BH.Q < 0.05, TF_cauchy_BH.Q < 0.05, peak_cauchy_BH.Q < 0.05) %>% 
  #  head() %>% 
  rowwise() %>% 
  dplyr::mutate(delta_phase_gene_TF = case_when(gene_phase - TF_phase > 12 ~ (gene_phase - TF_phase)-24, 
                                                gene_phase - TF_phase < -12 ~ (gene_phase - TF_phase)+24, 
                                                TRUE ~ gene_phase - TF_phase)) %>% 
  dplyr::mutate(delta_phase_peak_TF = case_when(peak_phase - TF_phase > 12 ~ (peak_phase-TF_phase)-24, 
                                                peak_phase - TF_phase < -12 ~ (peak_phase-TF_phase)+24, 
                                                TRUE ~ peak_phase - TF_phase)) %>% 
  dplyr::mutate(delta_phase_gene_peak = case_when(gene_phase - peak_phase > 12 ~ (gene_phase - peak_phase)-24, 
                                                  gene_phase - peak_phase < -12 ~ (gene_phase - peak_phase)+24,
                                                  TRUE ~ gene_phase - peak_phase)) %>% 
  ungroup() %>%
  #  dplyr::mutate(across(matches("delta_phase"), function(x){ifelse(x < -12, x+24, ifelse(x > 12, x-24, x))})) %>% 
  ungroup() -> trios_res_df_3


trios_res_df_3 %>% dplyr::filter(TF == "Arntl") %>% 
  dplyr::mutate(peak_phase = case_when(
    gene_phase - peak_phase > 12 ~ peak_phase + 24,
    gene_phase - peak_phase < -12 ~ peak_phase - 24,
    TRUE ~ peak_phase
  )) %>% 
  {
    df_ = .
    cor_ = cor(df_$gene_phase, df_$peak_phase, method = "pearson") %>% round(2)
    TF_phase = df_$TF_phase %>% unique()
    TF_phase - 12 -> TF_phase_12
    df_ %>% dplyr::filter(gene_phase>TF_phase-2.5, gene_phase<TF_phase+2.5) %>% .$gene %>% unique() ->> fast_regulated_genes
    df_ %>% dplyr::filter(gene_phase>TF_phase_12-2.5, gene_phase<TF_phase_12+2.5) %>% .$gene %>% unique() ->> latent_regulated_genes
    #    sprintf("Fast_regulated_genes: %s\nLatent_regulated_genes: %s", paste(fast_regulated_genes, collapse = ", "), paste(latent_regulated_genes, collapse = ", ")) %>% print()
    
    df_$group = ifelse(df_$gene %in% Core_clock_genes, "Core_clock_genes", "Other_genes")
    df_[, c("TF", "peak", "gene", "TF_phase", "peak_phase", "gene_phase", "group")] %>% distinct() -> df_
    print(df_ %>% dplyr::filter(gene == "Nr1d1"))
    df_ %>% 
      ggplot(aes(x = gene_phase, y = peak_phase)) + 
      #  geom_hex(bins = 50) + 
      geom_point(size = 0.5) + 
      geom_text_repel(data = subset(df_, group == "Core_clock_genes"), 
                      aes(label = gene), 
                      nudge_y = 5, 
                      direction = "y") + 
      geom_abline(color = "black", linetype = 2) + 
#      geom_vline(xintercept = TF_phase, color = "blue", linetype = 2) + 
#      geom_ribbon(aes(xmin = TF_phase-2.5, xmax = TF_phase+2.5), fill = "blue", alpha = 0.2) + 
#      geom_ribbon(aes(xmin = (TF_phase-12)-2.5, xmax = (TF_phase-12)+2.5), fill = "green", alpha = 0.2) + 
      theme_classic() + 
      ggtitle(sprintf("r=%s", cor_))
  } -> p_arntl#Fig_5E

trios_res_df_3 %>% dplyr::filter(TF == "Nr1d1") %>% 
  dplyr::mutate(peak_phase = case_when(
    gene_phase - peak_phase > 12 ~ peak_phase + 24,
    gene_phase - peak_phase < -12 ~ peak_phase - 24,
    TRUE ~ peak_phase
  )) %>% 
  {
    df_ = .
    cor_ = cor(df_$gene_phase, df_$peak_phase, method = "pearson") %>% round(2)
    TF_phase = df_$TF_phase %>% unique()
    TF_phase - 12 -> TF_phase_12
    df_ %>% dplyr::filter(gene_phase>TF_phase-2.5, gene_phase<TF_phase+2.5) %>% .$gene %>% unique() ->> fast_regulated_genes
    df_ %>% dplyr::filter(gene_phase>TF_phase_12-2.5, gene_phase<TF_phase_12+2.5) %>% .$gene %>% unique() ->> latent_regulated_genes
    #    sprintf("Fast_regulated_genes: %s\nLatent_regulated_genes: %s", paste(fast_regulated_genes, collapse = ", "), paste(latent_regulated_genes, collapse = ", ")) %>% print()
    
    df_$group = ifelse(df_$gene %in% Core_clock_genes, "Core_clock_genes", "Other_genes")
    df_[, c("TF", "peak", "gene", "TF_phase", "peak_phase", "gene_phase", "group")] %>% distinct() -> df_
    print(df_ %>% dplyr::filter(gene == "Nr1d1"))
    df_ %>% 
      ggplot(aes(x = gene_phase, y = peak_phase)) + 
      #  geom_hex(bins = 50) + 
      geom_point(size = 0.5) + 
      geom_text_repel(data = subset(df_, group == "Core_clock_genes"), 
                      aes(label = gene), 
                      nudge_y = 5, 
                      direction = "y") + 
      geom_abline(color = "black", linetype = 2) + 
#      geom_vline(xintercept = TF_phase, color = "blue", linetype = 2) + 
#      geom_ribbon(aes(xmin = TF_phase-2.5, xmax = TF_phase+2.5), fill = "blue", alpha = 0.2) + 
#      geom_ribbon(aes(xmin = (TF_phase+12)-2.5, xmax = (TF_phase+12)+2.5), fill = "green", alpha = 0.2) + 
      theme_classic() + 
      ggtitle(sprintf("r=%s", cor_))
  } -> p_nr1d1 #Fig_5E
patchwork::wrap_plots(p_arntl, p_nr1d1) #Fig_5E
####

## Add HC_C annotation ##
trios_res_df_2 %>% 
  {
    df_ = .
    GenomicRanges::GRanges(
      seqnames = df_$peak %>% gsub("(.+)-.+-.+", "\\1", .), 
      IRanges::IRanges(
        start = df_$peak %>% gsub("(.+)-(.+)-.+", "\\2", .) %>% as.integer(), 
        end = df_$peak %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer() 
      )
    ) %>% unique() 
  } -> trios_gr

GenomicRanges::GRanges(
  seqnames = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):.+:.+", "\\1", .),
  IRanges::IRanges(
    start = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):(.+):.+", "\\2", .) %>% as.integer(),
    end = GSE155161_annotated_CHI_C$CRE %>% gsub("(.+):(.+):(.+)", "\\3", .) %>% as.integer()
  ) 
) %>% unique() -> GSE155161_annotated_CHI_C_gr

ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, GSE155161_annotated_CHI_C_gr, maxgap = 1) -> ol
ol$overlappingPeaks[[1]][,c(2,3,4,8,9,10)] %>% as.data.frame() -> ol
data.frame(peak_trio = sprintf("%s-%s-%s", ol$seqnames, ol$start, ol$end), 
           peak_HI_C = sprintf("%s:%s:%s", ol$seqnames.1, ol$start.1, ol$end.1)) -> ol

GSE155161_annotated_CHI_C %>% dplyr::select(-ZT, -idx) %>% distinct() -> GSE155161_annotated_CHI_C_1

unique(ol$peak_trio) %>% 
  map(function(x){
    peak_trio_ = x
    ol %>% 
      dplyr::filter(peak_trio == peak_trio_) -> ol
    ol$peak_HI_C %>% 
      map(function(x){
        GSE155161_annotated_CHI_C_1 %>% 
          dplyr::filter(CRE == x) %>% 
          .$gene %>% unique() -> gene_
      }) %>% unlist() %>% unique() %>% {.[!is.na(.)]} %>% paste(., collapse = ";") -> gene_
    
    data.frame(peak_trio = peak_trio_, CHI_C_gene = gene_)
  }) %>% do.call(rbind, .) -> ol_1

trios_res_df_2 %>% 
  left_join(x = ., y = ol_1, by = c("peak"="peak_trio")) %>% 
  rowwise() %>% 
  dplyr::mutate(CHI_C_validate = case_when(
    any(gene %in% (CHI_C_gene %>% str_split(";") %>% unlist())) ~ TRUE, 
    TRUE ~ FALSE
  )) -> trios_res_df_2_1
####

#Add ChIP-seq annotation
load("~/Dropbox/singulomics/github_rda/TRIPOD/GSE101115_Rora_Rorg_ChIP/GSE101115_peak_list.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/GSE26345_Nr1d1_ChIP/GSE26345_peak_list.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/GSE59486_Nfil3_Rora_ChIP/GSE59486_peak_list.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/PRJDB7796_DBP_E4BP4_ChIP/PRJDB7796_peak_list.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/list_GSE39977_peak.RData")
rm(list_GSE39977_peak_df, list_GSE39977_peak_df_mm9)

list_GSE39977_peak_df_mm10 %>% 
  purrr::map(function(x){
    gene_ = x
    x[["q_0.1"]] -> df_
    return(df_)
  }) -> GSE39977_peak_list
GSE39977_peak_list %>% names()
GSE39977_peak_list[c("BMAL1_All", "CLOCK_ChIP-seq", "NPAS2_ChIP-seq")] -> GSE39977_peak_list
rm(list_GSE39977_peak_df_mm10)

names(GSE101115_peak_list) = c("Rora", "Rorc")
GSE26345_peak_list = list(Nr1d1 = GSE26345_peak_list)
names(GSE39977_peak_list) = c("Arntl", "Clock", "Npas2")
names(GSE59486_peak_list) = c("Nfil3", "Rora")
names(PRJDB7796_peak_list) = c("Dbp", "Nfil3")

ChIP_peak_list = list(
  Arntl = GSE39977_peak_list$Arntl,
  Clock = GSE39977_peak_list$Clock, 
  Npas2 = GSE39977_peak_list$Npas2,
  Rora = rbind(GSE59486_peak_list$Rora, GSE101115_peak_list$Rora),
  Rorc = GSE101115_peak_list$Rorc,
  Dbp = PRJDB7796_peak_list$Dbp,
  Nfil3 = rbind(PRJDB7796_peak_list$Nfil3, GSE59486_peak_list$Nfil3),
  Nr1d1 = GSE26345_peak_list$Nr1d1
)

ChIP_peak_list %>% 
  map(function(x){
    x %>% bedtoolsr::bt.sort() %>% bedtoolsr::bt.merge() %>% as_tibble() -> x
  }) -> ChIP_peak_list

ChIP_peak_list$Arntl %>% write.table(., file = "~/Downloads/Bmal1_ChIP_seq_peaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
ChIP_peak_list$Clock %>% write.table(., file = "~/Downloads/Clock_ChIP_seq_peaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)

c(names(ChIP_peak_list), "Others") %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  purrr::map(function(x){
    gene_ = x
    print(gene_)
    if (gene_ == "Others"){
      trios_res_df_2_1 %>% dplyr::filter(!(TF %in% names(ChIP_peak_list))) -> df_
      df_$ChIP_seq_validated = "FALSE"
      df_final = df_
    }else{
      trios_res_df_2_1 %>% dplyr::filter(TF == gene_) -> df_
      df_$peak %>% unique() -> peak_
      GenomicRanges::GRanges(seqnames = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\1", .), 
                             ranges = IRanges::IRanges(start = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer(), 
                                                       end = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer())) -> trios_gr
      ChIP_peak_list[[gene_]] %>% "colnames<-"(.,c("seqnames", "start", "end")) %>% 
        GenomicRanges::makeGRangesFromDataFrame() -> ChIP_gr
      ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, ChIP_gr, maxgap = 1) -> ol
      ol$overlappingPeaks[[1]] -> ol
      sprintf("%s-%s-%s", 
              ol[[2]] %>% as.character(), 
              ol[[3]], 
              ol[[4]]) -> validated_peak
      df_ %>% dplyr::mutate(ChIP_seq_validated = case_when(peak %in% validated_peak ~ "TRUE", TRUE ~ "FALSE")) -> df_final
    }
    return(df_final)
  }) %>% do.call(rbind, .) -> trios_res_df_2_1
trios_res_df_2_1 %>% as.data.frame() -> trios_res_df_2_1
####

## Make Supplementary Table 4 ##
trios_res_df_2_1 %>% 
  dplyr::select(gene, peak, TF, coef, pval, adj, model, sign, level, gene.pos, gene_strand, tf.pos, peak_to_tss_dist, peak_biotype, gene_phase, TF_cauchy_BH.Q, 
                TF_phase, peak_cauchy_BH.Q, peak_phase, chromatin_state_15, chromatin_state_18, CHI_C_validate, ChIP_seq_validated) %>% 
  dplyr::filter((CHI_C_validate == TRUE|ChIP_seq_validated==TRUE)) %>% 
  write.csv(x = ., file = "~/Downloads/00_Supp_Table_4_new.csv", quote = F, col.names = T, row.names = F) #Table_S4
####

## Plot ARNTL target gene trios ##
plot_trios = function(TF_, target_gene_, CRE_){
  
  metacell.motif_score[, TF_] %>% scales::rescale(., to = c(0,1)) %>% 
    as.data.frame() %>% rownames_to_column() %>% "colnames<-"(., c("metacell", "expr")) %>% 
    left_join(x = ., y = metacell_ZT, by = "metacell") -> df_
  df_$ZT %>% unique() %>% "names<-"(.,.) %>% 
    map(function(ZT_){
      df_ %>% dplyr::filter(ZT == ZT_) -> df_
      data.frame(ZT = ZT_, Mean = mean(df_$expr), SD = sd(df_$expr))
    }) %>% do.call(rbind, .) %>% dplyr::mutate(group = "TF") -> TF_motif_score_df
  
  metacell.rna[, target_gene_] %>% scales::rescale(., to = c(0,1)) %>% 
    as.data.frame() %>% rownames_to_column() %>% "colnames<-"(., c("metacell", "expr")) %>% 
    left_join(x = ., y = metacell_ZT, by = "metacell") -> df_
  df_$ZT %>% unique() %>% "names<-"(.,.) %>% 
    map(function(ZT_){
      df_ %>% dplyr::filter(ZT == ZT_) -> df_
      data.frame(ZT = ZT_, Mean = mean(df_$expr), SD = sd(df_$expr))
    }) %>% do.call(rbind, .) %>% dplyr::mutate(group = "Gene") -> RNA_expr_df
  
  metacell.peak[, CRE_] %>% scales::rescale(., to = c(0,1)) %>% 
    as.data.frame() %>% rownames_to_column() %>% "colnames<-"(., c("metacell", "expr")) %>% 
    left_join(x = ., y = metacell_ZT, by = "metacell") -> df_
  df_$ZT %>% unique() %>% "names<-"(.,.) %>% 
    map(function(ZT_){
      df_ %>% dplyr::filter(ZT == ZT_) -> df_
      data.frame(ZT = ZT_, Mean = mean(df_$expr), SD = sd(df_$expr))
    }) %>% do.call(rbind, .) %>% dplyr::mutate(group = "CRE") -> CRE_expr_df
  
  rbind(TF_motif_score_df, RNA_expr_df, CRE_expr_df) %>% dplyr::mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.numeric()) -> df_
  df_ %>% 
    ggplot(aes(x = ZT, y = Mean, group = group, color = group, fill = group)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Mean-SD, ymax = Mean+SD), alpha = 0.2, color = NA) + 
    scale_x_continuous(breaks = seq(2,22,4)) + 
    ggtitle(sprintf("TF: %s\nTarget gene: %s\nCRE: %s", TF_, target_gene_, CRE_)) 
}

trios_res_df_2_1 %>% dplyr::filter(TF == "Mef2d", gene == "Arntl")
plot_trios(TF_ = "Mef2d", target_gene_ = "Arntl", CRE_ = "chr7-113250569-113251463") + theme_classic() -> p1
trios_res_df_2_1 %>% dplyr::filter(TF == "Ppard", gene == "Arntl")
plot_trios(TF_ = "Ppard", target_gene_ = "Arntl", CRE_ = "chr7-113214355-113215153")
plot_trios(TF_ = "Ppard", target_gene_ = "Arntl", CRE_ = "chr7-113239502-113240409")
plot_trios(TF_ = "Ppard", target_gene_ = "Arntl", CRE_ = "chr7-113250569-113251463") + theme_classic() -> p2
trios_res_df_2_1 %>% dplyr::filter(TF == "Nr5a2", gene == "Arntl")
plot_trios(TF_ = "Nr5a2", target_gene_ = "Arntl", CRE_ = "chr7-113219262-113220256") + theme_classic() -> p3
trios_res_df_2_1 %>% dplyr::filter(TF == "Nfyb", gene == "Arntl")
plot_trios(TF_ = "Nfyb", target_gene_ = "Arntl", CRE_ = "chr7-113230019-113230865")
plot_trios(TF_ = "Nfyb", target_gene_ = "Arntl", CRE_ = "chr7-113242462-113243284") + theme_classic() -> p4
patchwork::wrap_plots(p1, p2, p3, p4, nrow = 2, guides = "collect") #Fig_S25

## Plot trios (E-box, D-box, RORE) ##
plot_trios(TF_ = "Clock", target_gene_ = "Per1", CRE_ = "chr11-69094571-69095496") -> p1
plot_trios(TF_ = "Rora", target_gene_ = "Npas2", CRE_ = "chr1-39192324-39193159") -> p2
plot_trios(TF_ = "Dbp", target_gene_ = "Per2", CRE_ = "chr1-91459045-91459906") -> p3
patchwork::wrap_plots(p1, p2, p3, nrow = 1, guides = "collect") #Fig_5B

####


## Plot case study, ARNTL/CLOCK regulate Nr1d1 ##
## load bigwig files ##
library(trackplot)
Sys.setenv(PATH = paste("/usr/local/mysql/bin", Sys.getenv("PATH"), sep = ":"))

list.files("~/Dropbox/singulomics/github_rda/bigwig/2014_Cell_Fang_Histone_ChIP", full.names = T) -> Histone_marker_bw
Histone_marker_bw = read_coldata(bws = Histone_marker_bw, build = "mm10")
Histone_marker_bw$bw_sample_names = c("Dnase-seq", "PolII-ChIP", "H3K4me1", "H3K27ac")
Histone_marker_bw[1, ] -> Dnase_seq_bw
Histone_marker_bw[-c(1:2), ] -> Histone_marker_bw

list.files("~/Dropbox/singulomics/github_rda/bigwig/2014_Science_Koike_Bmal1_ChIP", full.names = T) %>% 
  .[grepl("CT", .)] %>% gtools::mixedsort() -> Bmal1_bw
Bmal1_bw = read_coldata(bws = Bmal1_bw, build = "mm10")
Bmal1_bw$bw_sample_names = sprintf("Bmal1_ChIP_CT%s", seq(0,20,4))

list.files("~/Dropbox/singulomics/github_rda/bigwig/2014_Science_Koike_Clock_ChIP", full.names = T) %>% 
  gtools::mixedsort() -> Clock_bw
Clock_bw = read_coldata(bws = Clock_bw, build = "mm10")
Clock_bw$bw_sample_names = sprintf("Clock_ChIP_CT%s", seq(0,20,4))

list.files("~/Dropbox/singulomics/github_rda/bigwig/2014_Science_Koike_Npas2_ChIP", full.names = T) %>% 
  gtools::mixedsort() -> Npas2_bw
Npas2_bw = read_coldata(bws = Npas2_bw, build = "mm10")
Npas2_bw$bw_sample_names = sprintf("Npas2_ChIP_CT%s", seq(0,20,4))

list.files("~/Dropbox/singulomics/github_rda/bigwig/sc_multiomics_liver_ATAC", full.names = T) -> ATAC_bw
ATAC_bw = read_coldata(bws = ATAC_bw, build = "mm10")
ATAC_bw$bw_sample_names = sprintf("ATAC_ZT%s", seq(2,22,4))

list.files("~/Dropbox/singulomics/github_rda/bigwig/sc_multiomics_liver_RNA/", full.names = T) -> RNA_bw
RNA_bw = read_coldata(bws = RNA_bw, build = "mm10")
RNA_bw$bw_sample_names = sprintf("RNA_ZT%s", seq(2,22,4))
####

## Fig_5C ##
plot_trios(TF_ = "Clock", target_gene_ = "Nr1d1", CRE_ = "chr11-98747104-98747945") -> p1
plot_trios(TF_ = "Arntl", target_gene_ = "Nr1d1", CRE_ ="chr11-98775032-98775930") -> p2

c("chr11:98747104-98747945", "chr11:98775032-98775930") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  map(function(CRE_){
    c("Clock", "Bmal1") %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(TF_){
        if (TF_ == "Bmal1"){
          TF_CRE = track_extract(colData = Bmal1_bw, loci = CRE_) 
        }else{
          TF_CRE = track_extract(colData = Clock_bw, loci = CRE_) 
        }
        names(TF_CRE$data) %>% 
          map(function(timepoint_){
            TF_CRE$data[[timepoint_]] -> df_1
            df_1$max %>% sum(na.rm = T) %>% sum() -> measure_
            data.frame(TF = TF_, CRE = CRE_, measure = measure_, timepoint = timepoint_)
          }) %>% do.call(rbind, .) %>% dplyr::mutate(measure = scales::rescale(measure, to = c(0,1)))
      }) %>% do.call(rbind, .)
  })  %>% do.call(rbind, .) %>% 
  dplyr::mutate(timepoint = gsub(".+_CT(\\d+)", "\\1", timepoint) %>% as.numeric()) -> df_

df_ %>% 
  ggplot(aes(x = timepoint, y = measure, color = TF, linetype = CRE)) + 
  geom_line() + theme_classic() + 
  scale_x_continuous(breaks = seq(2,22,4)) -> p_chipseq

Plot_CHI_C = function(target_gene, Peak){
  GSE155161_annotated_CHI_C %>% dplyr::filter(gene == target_gene) -> df_
  data.frame(seqnames = gsub("(.+):.+:.+", "\\1", df_$CRE), 
             start = gsub("(.+):(.+):.+", "\\2", df_$CRE) %>% as.integer(), 
             end = gsub("(.+):(.+):(.+)", "\\3", df_$CRE) %>% as.integer(), 
             ZT = df_$ZT, 
             score = df_$score) -> df_
  
  data.frame(seqnames = gsub("(.+)-.+-.+", "\\1", Peak), 
             start = gsub("(.+)-(.+)-.+", "\\2", Peak) %>% as.integer(),
             end = gsub("(.+)-(.+)-(.+)", "\\3", Peak) %>% as.integer()) %>% 
    GenomicRanges::makeGRangesFromDataFrame(.) -> peak_gr
  
  df_$ZT %>% unique() %>% 
    gtools::mixedsort() %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      ZT_ = x
      df_ %>% dplyr::filter(ZT == ZT_) -> df_
      df_ %>% 
        group_by(seqnames, start, end) %>% 
        group_map(function(df_, y){
          df_ %>% dplyr::mutate(score = mean(score)) %>% distinct() -> df_
        }, .keep = T) %>% do.call(rbind, .) -> df_
      GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> CHI_C_gr
      ChIPpeakAnno::findOverlapsOfPeaks(peak_gr, CHI_C_gr, maxgap = 1) -> ol
      ol$overlappingPeaks[[1]] -> ol
      print(ol)
      ol$score %>% sum() -> interaction_score
      data.frame(Peaks = Peak, ZT = gsub("ZT(\\d+)", "\\1", ZT_) %>% as.integer(), measure = interaction_score)
      #      data.frame(group = "CHI_C", ZT = gsub("ZT(\\d+)", "\\1", ZT_) %>% as.integer(), measure = interaction_score)
      #      data.frame(group_ = "CHI_C", ZT_ = gsub("ZT(\\d+)", "\\1", ZT_) %>% as.integer(), measure_ = interaction_score)
    }) %>% do.call(rbind, .) -> CHI_C_dat
  
  CHI_C_dat %>% 
    ggplot(aes(x = ZT, y = measure)) + 
    geom_line() -> p
  return(p)
}

c("chr11-98747104-98747945", "chr11-98775032-98775930") %>% 
  "names<-"(.,.) %>% 
  map(function(CRE_){
    Plot_CHI_C(target_gene = "Nr1d1", Peak = CRE_) -> p
    p$data
  }) %>% do.call(rbind, .) -> df_
df_ %>% 
  ggplot(aes(x = ZT, y = measure, linetype = Peaks)) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(2,22,4)) + 
  theme_classic() -> p_HI_C

readxl::read_xlsx(path = "~/Dropbox/singulomics/MS22_Greenwell_Total_Nuclear_polysome_BMKO_v1.xlsx", sheet = 1) -> df_
df_[4:nrow(df_), ] %>% 
  {
    df_ = .
    df_[1,1] = "Gene"
    df_
  } %>% 
  "colnames<-"(., .[1,]) %>% 
  .[-1, ] %>% 
  dplyr::select(matches("Gene|KO\\d+ZT\\d+")) %>% 
  dplyr::mutate(across(matches("KO"), function(x){as.numeric(x)})) -> Bmal1_KO_df

df_[4:nrow(df_), ] %>% 
  {
    df_ = .
    df_[1,1] = "Gene"
    df_
  } %>% 
  "colnames<-"(., .[1,]) %>% 
  .[-1, ] %>% 
  dplyr::select(matches("Gene|NU\\d+ZT\\d+")) %>% 
  dplyr::mutate(across(matches("NU"), function(x){as.numeric(x)})) -> Bmal1_WT_df

source("~/Dropbox/singulomics/github/Calculate_HMP.R")
Bmal1_KO_df$Gene %>% table() %>% {.[. > 1]} %>% names() -> excluded_genes
Bmal1_KO_df %>% 
  dplyr::filter(!(Gene %in% excluded_genes)) %>% 
  {
    df_ = .
    gene_ = df_$Gene
    timepoint_ = colnames(df_)[-1] %>% gsub("KO\\d+ZT(\\d+)", "\\1", .) %>% as.integer()
    rep_ = colnames(df_)[-1] %>% gsub("KO(\\d+)ZT\\d+", "\\1", .) %>% as.integer()
    exp_mat_ = df_ %>% dplyr::select(-Gene)
    sprintf("ZT%s_REP%s", timepoint_, rep_) -> colnames(exp_mat_)
    colnames(exp_mat_)
    cyclic_HMP(exp_matrix = exp_mat_, gene = gene_, timepoint = timepoint_) -> res_
    res_
  } -> res_pval_KO
res_pval_KO %>% dplyr::select(matches("Gene|JTK|HR")) %>% 
  recal_cauchy_p_and_hmp() -> res_pval_KO

Bmal1_WT_df$Gene %>% table() %>% {.[. > 1]} %>% names() -> excluded_genes
Bmal1_WT_df %>% 
  dplyr::filter(!(Gene %in% excluded_genes)) %>% 
  {
    df_ = .
    gene_ = df_$Gene
    timepoint_ = colnames(df_)[-1] %>% gsub("NU\\d+ZT(\\d+)", "\\1", .) %>% as.integer()
    rep_ = colnames(df_)[-1] %>% gsub("NU(\\d+)ZT\\d+", "\\1", .) %>% as.integer()
    exp_mat_ = df_ %>% dplyr::select(-Gene)
    sprintf("ZT%s_REP%s", timepoint_, rep_) -> colnames(exp_mat_)
    colnames(exp_mat_)
    cyclic_HMP(exp_matrix = exp_mat_, gene = gene_, timepoint = timepoint_) -> res_
    res_
  } -> res_pval_WT
res_pval_WT %>% dplyr::select(matches("Gene|JTK|HR")) %>% 
  recal_cauchy_p_and_hmp() -> res_pval_WT

list(WT = Bmal1_WT_df, KO = Bmal1_KO_df) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    df_ = x
    df_ %>% dplyr::filter(Gene == "Nr1d1") %>% 
      pivot_longer(cols = -Gene, names_to = "ZT", values_to = "measure") %>% 
      dplyr::mutate(Rep = gsub(".+?(\\d+)ZT(\\d+)", "\\1", ZT)) %>% 
      dplyr::mutate(ZT = gsub(".+?(\\d+)ZT(\\d+)", "\\2", ZT) %>% as.integer()) %>% 
      group_by(Gene, ZT) %>% 
      summarise(expression = mean(measure), sd = sd(measure)) %>% 
      ungroup() %>% 
      dplyr::mutate(group = group_)
  }) %>% do.call(rbind, .) %>% 
  ggplot(aes(x = ZT, y = expression, fill = group, color = group)) +
  geom_line() + 
  geom_ribbon(aes(ymin = expression - sd, ymax = expression + sd), alpha = 0.3, color = NA) + 
  scale_x_continuous(breaks = seq(2,22,4)) +
  theme_classic() -> p_KO #Fig_5C

res_pval_KO %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Nr1d1")
res_pval_WT %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Nr1d1")


p1+p2+p_chipseq+p_HI_C+p_KO #Fig_5C

## Plot track plots Fig_5D##
extend_region = function(upstream, downstream, CRE){
 start = gsub(".+:(.+?)-(.+)", "\\1", CRE) %>% as.integer() 
 end = gsub(".+:(.+?)-(.+)", "\\2", CRE)  %>% as.integer()
 chr = gsub("(.+):(.+?)-(.+)", "\\1", CRE)
 new_start = start - upstream
 new_end = end + downstream
 new_CRE = sprintf("%s:%s-%s", chr, new_start, new_end)
 return(new_CRE)
}
range_to_bed = function(range){
  start = gsub(".+:(.+?)-(.+)", "\\1", range) %>% as.integer() 
  end = gsub(".+:(.+?)-(.+)", "\\2", range)  %>% as.integer()
  chr = gsub("(.+):(.+?)-(.+)", "\\1", range) 
  data.frame(seqnames = chr, start = start, end = end) -> df_
  return(df_)
}

CREs = c("chr11:98747104-98747945", "chr11:98775032-98775930")
CRE = "chr11:98747104-98747945"
Nr1d1 = "chr11:98765932-98777377"

extended_CRE = extend_region(upstream = 1000, downstream = 1000, CRE = CRE)
extended_gene = extend_region(upstream = 2000, downstream = 2000, CRE = Nr1d1)
range_to_bed(range = CREs) -> bed
write.table(x = bed, file = "~/Downloads/tmp.bed", col.names = F, row.names = F, quote = F, sep = "\t")

library(RColorBrewer)
colors <- brewer.pal(n = 6, name = "Set2")

#Histone
Histone_marker_t_CRE = track_extract(colData = Histone_marker_bw, loci = extended_CRE)
track_plot(summary_list = Histone_marker_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed")

Histone_marker_t_Nr1d1 = track_extract(colData = Histone_marker_bw, loci = extended_gene)
track_plot(summary_list = Histone_marker_t_Nr1d1, show_ideogram = F, peaks = "~/Downloads/tmp.bed")
###

#Dnase
Dnase_t_CRE = track_extract(colData = Dnase_seq_bw, loci = extended_CRE)
track_plot(summary_list = Dnase_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed")

Dnase_t_Nr1d1 = track_extract(colData = Dnase_seq_bw, loci = extended_gene)
track_plot(summary_list = Dnase_t_Nr1d1, show_ideogram = F, peaks = "~/Downloads/tmp.bed")
###

#Bmal1
Bmal1_t_CRE = track_extract(colData = Bmal1_bw, loci = extended_CRE)
track_plot(summary_list = Bmal1_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed", track_overlay = F, col = colors)

Bmal1_t_Nr1d1 = track_extract(colData = Bmal1_bw, loci = extended_gene)
track_plot(summary_list = Bmal1_t_Nr1d1, show_ideogram = F, track_overlay = F, col = colors)
####

#Clock
Clock_t_CRE = track_extract(colData = Clock_bw, loci = extended_CRE)
track_plot(summary_list = Clock_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed", track_overlay = F, col = colors)

Clock_t_Nr1d1 = track_extract(colData = Clock_bw, loci = extended_gene)
track_plot(summary_list = Clock_t_Nr1d1, show_ideogram = F, track_overlay = F, peaks = "~/Downloads/tmp.bed", col = colors)
####

#ATAC
ATAC_t_CRE = track_extract(colData = ATAC_bw, loci = extended_CRE)
track_plot(summary_list = ATAC_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed", track_overlay = F, col = colors)

ATAC_t_Nr1d1 = track_extract(colData = ATAC_bw, loci = extended_gene)
track_plot(summary_list = ATAC_t_Nr1d1, show_ideogram = F, track_overlay = F, col = colors, peaks = "~/Downloads/tmp.bed")
####

#RNA
RNA_t_CRE = track_extract(colData = RNA_bw, loci = extended_CRE)
track_plot(summary_list = RNA_t_CRE, show_ideogram = F, peaks = "~/Downloads/tmp.bed", track_overlay = F, col = colors)

RNA_t_Nr1d1 = track_extract(colData = RNA_bw, loci = extended_gene)
track_plot(summary_list = RNA_t_Nr1d1, show_ideogram = F, track_overlay = F, col = colors)