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

genes_ = sc@assays$SCT@counts %>% rownames()
ext.upstream <- ext.downstream <- 2e5
xymats.list <- bplapply(
  genes_,
  getXYMatrices,
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
save(xymats.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.list.rda')


xymats.m.list <- bplapply(
  xymats.list,
  fitModel,
  model.name = "marginal"
)
save(xymats.m.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.m.list.rda')
rm(xymats.m.list)

xymats.c.list <- bplapply(
  xymats.list,
  fitModel,
  model.name = "conditional"
)
save(xymats.c.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.c.list.rda')
rm(xymats.c.list)

xymats.i.list <- bplapply(
  xymats.list,
  fitModel,
  model.name = "interaction"
)
save(xymats.i.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.i.list.rda')
rm(xymats.i.list)

xymats.tripod.Xt.list <- bplapply(
  xymats.list,
  fitModel,
  model.name = "TRIPOD",
  match.by = "Xt"
)
save(xymats.tripod.Xt.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Xt.list.rda')

xymats.tripod.Yj.list <- bplapply(
  xymats.list,
  fitModel,
  model.name = "TRIPOD",
  match.by = "Yj"
)
save(xymats.tripod.Yj.list, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tripod.Yj.list.rda')

# set FDR < 0.01
fdr.thresh <- 0.01
# focus on positive sign
sign <- "positive"

xymats.tX1.pos.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 1
)
xymats.tX1.pos.df[order(xymats.tX1.pos.df$adj), ] -> xymats.tX1.pos.df 
save(xymats.tX1.pos.df, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tX1.pos.df.rda')

xymats.tX2.pos.df <- getTrios(
  xymats.list = xymats.tripod.Xt.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "TRIPOD",
  level = 2
)
xymats.tX2.pos.df[order(xymats.tX2.pos.df$adj), ] -> xymats.tX2.pos.df
save(xymats.tX2.pos.df, file = '~/Dropbox/singulomics/github_rda/TRIPOD/xymats.tX2.pos.df.rda')

fdr.thresh <- 0.01
list(
  marginal = xymats.m.list,
  conditional.on.Yj = xymats.c.list,
  conditional.on.Xt = xymats.c.list,
  interaction = xymats.i.list,
  TRIPOD.Yj = xymats.tripod.Yj.list,
  TRIPOD.Xt = xymats.tripod.Xt.list
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
list_trios_res = list_res ; rm(list_res)
save(list_trios_res, file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res.rda')
rm(list=ls())
####

library(Seurat)
library(Signac)
library(tidyverse)

load(file='~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')


readRDS(file="~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
gene.ref %>% as.data.frame() -> gene.ref

sc@meta.data[,c("ZT", "seurat_clusters")] %>% distinct() %>% dplyr::arrange(seurat_clusters) %>% 
  mutate(metacell = sprintf("metacell_%s", seurat_clusters)) %>% 
  "rownames<-"(., 1:nrow(.)) -> metacell_ZT

source("~/Dropbox/singulomics/github/Calculate_HMP.R")
list(rna = metacell.rna, peak = metacell.peak) %>% 
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
save(list_pval, file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_circadian_pval.rda')

list_trios_res %>% 
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
  }) %>% do.call(rbind, .) -> trios_res_df

list(trios_res_df = trios_res_df) %>% 
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
  }) %>% .[[1]] -> trios_res_df

trios_res_df %>% 
#  head() %>% 
  mutate(model_1 = case_when(
    grepl("TRIPOD", model) ~ "TRIPOD", 
    grepl("conditional", model) ~ "conditional", 
    TRUE ~ model
  )) -> trios_res_df

# Add CRE annotaiton (peak_biotype)
readRDS("~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
gene.ref[gene.ref$gene_biotype == "lincRNA"] -> gene.ref.lincRNA
trios_res_df %>% 
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

trios_res_df %>% 
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

c(lincRNA_map$peak, eRNA_map$peak) %>% length()
c(lincRNA_map$peak, eRNA_map$peak) %>% unique() %>% length()
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

full_join(x = trios_res_df, y = map_df, by = "peak") -> trios_res_df_1

load(file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_circadian_pval.rda')
list_pval$rna %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("gene", "gene_cauchy_BH.Q", "gene_phase")) -> gene_map
list_pval$rna %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("TF", "TF_cauchy_BH.Q", "TF_phase")) -> TF_map
list_pval$peak %>%
  dplyr::select(Gene, cauchy_BH.Q, HR_phi) %>% 
  dplyr::mutate(HR_phase = (HR_phi/(2*pi))*24) %>% 
  .[,c(1,2,4)] %>% 
  "colnames<-"(., c("peak", "peak_cauchy_BH.Q", "peak_phase")) -> peak_map
  
dim(trios_res_df_1)
gene_map %>% 
  left_join(x = trios_res_df_1, y = ., by = "gene") %>% 
  left_join(x = ., y = TF_map, by = "TF") %>% 
  left_join(x = ., y = peak_map, by = "peak") -> trios_res_df_1

rm(gene_map, TF_map, peak_map)

phase_diff_cal = function(phase_1, phase_2){
  phase_diff = min(abs(phase_1-phase_2), 24-abs(phase_1-phase_2))
  return(phase_diff)
}
trios_res_df_1 %>% 
#  head(100) %>% 
  rowwise() %>%
  mutate(TF_peak_phase_diff = phase_diff_cal(TF_phase, peak_phase), 
         peak_gene_phase_diff = phase_diff_cal(peak_phase, gene_phase), 
         TF_gene_phase_diff = phase_diff_cal(TF_phase, gene_phase)) -> trios_res_df_1
trios_res_df_1 %>% 
  mutate(regulation = case_when(
    (TF_gene_phase_diff<=1.5)&(TF_peak_phase_diff<=1.5)&(peak_gene_phase_diff<=1.5) ~ "activation",
    (TF_gene_phase_diff>=10.5)&(TF_peak_phase_diff<=1.5)&(peak_gene_phase_diff>=10.5) ~ "repression",
    (TF_gene_phase_diff>=10.5)&(TF_peak_phase_diff>=10.5)&(peak_gene_phase_diff<=1.5) ~ "repression",
    TRUE ~ "loop/indirect"
  )) %>% 
#  View()
  as.data.frame() -> trios_res_df_1

trios_res_df_1 %>% 
  mutate(peak_biotype = case_when(
    is.na(peak_biotype) ~ "No annotation", 
    TRUE ~ peak_biotype
  )) -> trios_res_df_1
save(trios_res_df_1, file = "~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1.RData")

# Add chromatin states annotations
load("~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1.RData")

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

trios_res_df_1 %>% 
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
     file = "~/Dropbox/singulomics/github_rda/TRIPOD/ChromHMM_dat.RData")

left_join(x = trios_res_df_1, y = mapping_df_list$`15-State`, by = "peak") %>% 
  "colnames<-"(., gsub("^chromatin_state$", "chromatin_state_15", colnames(.))) %>% 
  left_join(x = ., y = mapping_df_list$`18-State`, by = "peak") %>% 
  "colnames<-"(., gsub("^chromatin_state$", "chromatin_state_18", colnames(.))) -> trios_res_df_1

trios_res_df_1 %>% 
  dplyr::mutate(chromatin_state_15 = case_when(
    is.na(chromatin_state_15) ~ "no_annotation",
    TRUE ~ chromatin_state_15
  )) %>% 
  dplyr::mutate(chromatin_state_18 = case_when(
    is.na(chromatin_state_18) ~ "no_annotation",
    TRUE ~ chromatin_state_18
  )) -> trios_res_df_1

save(trios_res_df_1, file = "~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1.RData")

# Add gene and TF annotations
mm10_gtf = read.table("~/Dropbox/singulomics/github_rda/TRIPOD/gencode.vM10.annotation.gtf", head = F, sep = "\t", stringsAsFactors = F, skip = 5)

mm10_gtf %>% 
  dplyr::mutate(gene_name = gsub(".+gene_name (.+?); .+", "\\1", V9), .after = 5) -> mm10_gtf

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 7)

mm10_gtf %>% dplyr::select(V1,V4,V5,V3, gene_name) %>% 
  "colnames<-"(., c("chr", "start", "end", "annotation", "gene_name")) %>% 
  mutate(peak = sprintf("%s-%s-%s", chr, start, end)) %>% 
  {
    df_ = .
    peak_unique = df_$peak %>% unique()
    peak_unique %>% 
      "names<-"(.,.) %>% 
      future_map(function(x){
        peak_ = x
        df_ %>% dplyr::filter(peak == peak_) -> df_
        chr_ = df_$chr %>% unique()
        start_ = df_$start %>% unique()
        end_ = df_$end %>% unique()
        annotation_ = df_$annotation %>% unique() %>% paste(., collapse = ";")
        gene_name_ = df_$gene_name %>% unique() %>% paste(., collapse = ";")
        data.frame(chr = chr_, start = start_, end = end_, annotation = annotation_, gene_name = gene_name_)
      }) %>% do.call(rbind, .) -> df_
    GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> gr
    gr
  } -> mm10_gtf_gr

colnames(metacell.peak) %>% 
  {
    peak_ = .
    chr_ = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\1", .)
    start_ = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer()
    end_ = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer()
    data.frame(chr = chr_, start = start_, end = end_) -> df_
    GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> gr
    gr
  } -> All_peak_gr

ChIPpeakAnno::findOverlapsOfPeaks(All_peak_gr, mm10_gtf_gr) -> ol 

library(furrr)
2000*1024^2 -> max_size
options(future.globals.maxSize= max_size)
plan(multisession, workers = 7)
ol$overlappingPeaks[[1]] %>% 
  {
    df_ = .
    colnames(df_) -> col_
    col_[2:6] = sprintf("%s_1", col_[2:6])
    col_[7:12] = sprintf("%s_2", col_[7:12])
    colnames(df_) = col_ 
    
    overlapped_bp_ = function(start_1, end_1, start_2, end_2){
      sort(c(start_1, end_1, start_2, end_2)) %>% {.[3]-.[2]} -> overlapped_bp
      return(overlapped_bp)
    }
    
    df_ %>% 
      rowwise() %>% 
      mutate(overlapped_bp = overlapped_bp_(start_1, end_1, start_2, end_2)) -> df_
    
    unique_peak_ = df_$peaks1 %>% unique()
    unique_peak_ %>% 
      future_map(function(x){
        peak_ = x
        df_ %>% dplyr::filter(peaks1 == peak_) -> df_
        df_$overlapped_bp %>% max() -> max_overlapped_bp
        df_ %>% dplyr::filter(overlapped_bp == max_overlapped_bp) -> df_
        chr_ = df_$seqnames_1 %>% unique()
        start_ = df_$start_1 %>% unique()
        end_ = df_$end_1 %>% unique()
        annotation_ = df_$annotation %>% unique() %>% paste(., collapse = ";")
        gene_ = df_$gene_name %>% unique() %>% paste(., collapse = ";")
        data.frame(chr = chr_, start = start_, end = end_, gtf_annotation = annotation_, gtf_gene = gene_)
      }) %>% do.call(rbind ,.)
  } %>% dplyr::mutate(peak = sprintf("%s-%s-%s", chr, start, end)) %>% 
  dplyr::select(peak, gtf_annotation, gtf_gene) -> gtf_mapping_df

gtf_mapping_df %>% 
  rowwise() %>%
  dplyr::mutate(gtf_annotation = gtf_annotation %>% str_split(";") %>% unlist() %>% unique() %>% paste(., collapse = ";")) %>% 
  dplyr::mutate(gtf_gene = gtf_gene %>% str_split(";") %>% unlist() %>% unique() %>% paste(., collapse = ";")) -> gtf_mapping_df

save(mm10_gtf_gr, All_peak_gr, gtf_mapping_df, file = "~/Dropbox/singulomics/github_rda/TRIPOD/mm10_gtf_dat.rda")

load("~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1.RData")
left_join(x = trios_res_df_1, y = gtf_mapping_df, by = "peak") %>% 
  dplyr::mutate(gtf_annotation = case_when(
    is.na(gtf_annotation) ~ "intergenic", 
    TRUE ~ gtf_annotation
  )) %>% 
  dplyr::mutate(gtf_gene = case_when(
    is.na(gtf_gene) ~ "intergenic", 
    TRUE ~ gtf_gene
  )) -> trios_res_df_1
save(trios_res_df_1, file = "~/Dropbox/singulomics/github_rda/TRIPOD/trios_res_df_1.RData")
####

# Add promoter enriched HI-C data and annotation
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

trios_res_df_1 %>% 
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

ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, GSE155161_annotated_CHI_C_gr, maxgap = 1000) -> ol
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

trios_res_df_1 %>% 
  left_join(x = ., y = ol_1, by = c("peak"="peak_trio")) %>% 
  rowwise() %>% 
  dplyr::mutate(CHI_C_validate = case_when(
    any(gene %in% (CHI_C_gene %>% str_split(";") %>% unlist())) ~ TRUE, 
    TRUE ~ FALSE
  )) -> trios_res_df_2
trios_res_df_1 = trios_res_df_2; rm(trios_res_df_2)
save(trios_res_df_1, file = '~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res_processed.rda')

# Add ChIP seq data and annotation
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

#Process ChIP peak from GSE39977 (Arntl, Npas2, Clock) 

read.table("~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/mapping.txt", header = F, sep = "\t", stringsAsFactors = F) %>% 
  "colnames<-"(., c("GSM", "sample")) -> mapping_df

 ####macs2 peak calling ----
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

list_GSE39977_peak_df_mm10 %>% 
  purrr::map(function(x){
    gene_ = x
    x[["q_0.1"]] -> df_
    return(df_)
  }) -> GSE39977_peak_list
GSE39977_peak_list %>% names()
GSE39977_peak_list[c("BMAL1_All", "CLOCK_ChIP-seq", "NPAS2_ChIP-seq")] -> GSE39977_peak_list
rm(list_GSE39977_peak_df_mm10)

names(GSE101115_peak_list) #"Rora" "Rorg"
names(GSE26345_peak_list) #"Nr1d1"
names(GSE39977_peak_list) #"Bmal1" "Clock" "Npas2"
names(GSE59486_peak_list) # "E4bp4" "Rora"
names(PRJDB7796_peak_list) #"DBP" "E4BP4"

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

c(names(ChIP_peak_list), "Others") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  purrr::map(function(x){
    gene_ = x
    print(gene_)
    if (gene_ == "Others"){
      trios_res_df_1 %>% dplyr::filter(!(TF %in% names(ChIP_peak_list))) -> df_
      df_$ChIP_seq_validated = "FALSE"
      df_final = df_
    }else{
      trios_res_df_1 %>% dplyr::filter(TF == gene_) -> df_
      df_$peak %>% unique() -> peak_
      GenomicRanges::GRanges(seqnames = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\1", .), 
                             ranges = IRanges::IRanges(start = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\2", .) %>% as.integer(), 
                                                      end = peak_ %>% gsub("(.+)-(.+)-(.+)", "\\3", .) %>% as.integer())) -> trios_gr
      ChIP_peak_list[[gene_]] %>% "colnames<-"(.,c("seqnames", "start", "end")) %>% 
        GenomicRanges::makeGRangesFromDataFrame() -> ChIP_gr
      ChIPpeakAnno::findOverlapsOfPeaks(trios_gr, ChIP_gr, maxgap = 1000) -> ol
      ol$overlappingPeaks[[1]] -> ol
      sprintf("%s-%s-%s", 
              ol[[2]] %>% as.character(), 
              ol[[3]], 
              ol[[4]]) -> validated_peak
      df_ %>% dplyr::mutate(ChIP_seq_validated = case_when(peak %in% validated_peak ~ "TRUE", TRUE ~ "FALSE")) -> df_final
    }
    return(df_final)
  }) %>% do.call(rbind, .) -> trios_res_df_2
save(trios_res_df_2, file = "~/Downloads/trios_res_df_2.RData")
trios_res_df_2 %>% as.data.frame() -> trios_res_df_2
####