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

# 4. TRIPOD results (downstream analysis) ----
library(Seurat)
library(Signac)
library(tidyverse)

load("~/Dropbox/singulomics/github_rda/TRIPOD/GSE155161_CHI-C_annotated.rda")
load("~/Dropbox/singulomics/github_rda/TRIPOD/GSE39977/GSE39977_ChIP_timepoint_mat.rda")
load("~/Downloads/trios_res_df_2.RData")
trios_res_df_2 %>% as.data.frame() -> trios_res_df_2
trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD", (CHI_C_validate==TRUE|ChIP_seq_validated==TRUE)) %>% 
  dplyr::select(-c(peak_num, TF_num, model_1, peak_name, regulation, gtf_annotation, gtf_gene, CHI_C_gene)) %>% 
  {write.csv(., file="~/Downloads/trios_res_df_2.csv", row.names = F, quote = F)}

#load(file='~/Dropbox/singulomics/github_rda/TRIPOD/list_trios_res.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')
load('~/Dropbox/singulomics/github_rda/TRIPOD/list_circadian_pval.rda')


readRDS(file="~/Dropbox/singulomics/github_rda/gene.ref.rds") -> gene.ref
gene.ref %>% as.data.frame() -> gene.ref

sc@meta.data[,c("ZT", "seurat_clusters")] %>% distinct() %>% dplyr::arrange(seurat_clusters) %>% 
  mutate(metacell = sprintf("metacell_%s", seurat_clusters)) %>% 
  "rownames<-"(., 1:nrow(.)) -> metacell_ZT
rm(sc)

Core_clock_genes = c("Dbp", "Arntl", "Bhlhe40", "Bhlhe41", "Nfil3", "Rorc", "Rora", "Nr1d1", "Clock", "Npas2", "Cry1", "Ciart", "Per1", "Per2")
PlotGenePeakTF_custom = function(gene, TF, peak, peak_name = NULL){
  
  unique(metacell_ZT$ZT) %>% 
    as.character() %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      ZT_ = x
      metacell_ZT %>% dplyr::filter(ZT == ZT_) %>% .$metacell -> metacell_
      metacell.rna[metacell_, gene] %>% 
        as.data.frame() %>% 
        "colnames<-"(.,"gene_expression") %>% 
        mutate(ZT = ZT_)
    }) %>% do.call(rbind, .) %>% 
    mutate(ZT = factor(ZT)) %>% 
    ggplot(aes(x = ZT, y = gene_expression, fill = ZT)) + 
    geom_boxplot() -> p1
  p1 + theme_classic() + 
    theme(legend.position = "none") + 
    xlab(NULL) -> p1
  padj_ = list_pval[["rna"]] %>% dplyr::filter(Gene == gene) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
  p1 + ggtitle(sprintf("Target gene: %s\np.adj: %s", gene, padj_)) -> p1
  
  unique(metacell_ZT$ZT) %>% 
    as.character() %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      ZT_ = x
      metacell_ZT %>% dplyr::filter(ZT == ZT_) %>% .$metacell -> metacell_
      metacell.peak[metacell_, peak] %>% 
        as.data.frame() %>% 
        "colnames<-"(.,"peak_accessibility") %>% 
        mutate(ZT = ZT_)
    }) %>% do.call(rbind, .) %>% 
    mutate(ZT = factor(ZT)) %>% 
    ggplot(aes(x = ZT, y = peak_accessibility, fill = ZT)) + 
    geom_boxplot() -> p2
  p2 + theme_classic() + 
    theme(legend.position = "none") + 
    xlab(NULL) -> p2
  padj_ = list_pval[["peak"]] %>% dplyr::filter(Gene == peak) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
  if (is.null(peak_name)){
    p2 + ggtitle(sprintf("CRE: %s\np.adj: %s", peak, padj_)) -> p2    
  }else{
    p2 + ggtitle(sprintf("%s: %s\np.adj: %s", peak_name, peak, padj_)) -> p2    
  }
  
  
  unique(metacell_ZT$ZT) %>% 
    as.character() %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      ZT_ = x
      metacell_ZT %>% dplyr::filter(ZT == ZT_) %>% .$metacell -> metacell_
      metacell.rna[metacell_, TF] %>% 
        as.data.frame() %>% 
        "colnames<-"(.,"gene_expression") %>% 
        mutate(ZT = ZT_)
    }) %>% do.call(rbind, .) %>% 
    mutate(ZT = factor(ZT)) %>% 
    ggplot(aes(x = ZT, y = gene_expression, fill = ZT)) + 
    geom_boxplot() -> p3
  
  p3 + theme_classic() + 
    theme(legend.position = "none") + 
    xlab(NULL) -> p3
  padj_ = list_pval[["rna"]] %>% dplyr::filter(Gene == TF) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
  p3 + ggtitle(sprintf("TF: %s\np.adj: %s", TF, padj_)) -> p3
  
  patchwork::wrap_plots(p1, p2, p3, ncol = 3)
}
PlotGenePeakTF_custom_1 = function(plot){
  p = plot
  gene_df = p[[1]]$data
  gene_df %>% 
    dplyr::mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.integer()) %>% 
    mutate(group = "target_gene") -> gene_df
  p[[1]]$label$title %>% gsub(".+?: (.+)\n.+", "\\1", .) -> target_gene
  
  TF_df = p[[3]]$data
  TF_df %>% 
    dplyr::mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.integer()) %>% 
    mutate(group = "TF") -> TF_df
  p[[3]]$label$title %>% gsub(".+?: (.+)\n.+", "\\1", .) -> TF
  
  #    c(gene_df$gene_expression, TF_df$gene_expression) %>% max() -> gene_max
  
  peak_df = p[[2]]$data
  peak_df %>% 
    dplyr::mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.integer()) %>% 
    mutate(group = "peak") %>% 
    "colnames<-"(.,c("gene_expression", "ZT", "group")) -> peak_df
  p[[2]]$label$title %>% gsub(".+?: (.+)\n.+", "\\1", .) -> peak
  
  #    peak_df$gene_expression %>% max() -> peak_max
  #    scale_factor = gene_max/peak_max
  #    
  rbind(gene_df, TF_df, peak_df) -> df_
  df_ %>% 
    group_by(group, ZT) %>% 
    group_map(function(x,y){
      measure = x$gene_expression %>% mean()
      sd = x$gene_expression %>% sd()
      data.frame(group = y$group, ZT = y$ZT, measure = measure, sd = sd) -> df_
      
      return(df_)
    }) %>% do.call(rbind, .) -> df_
  gene_max = df_ %>% dplyr::filter(group != "peak") %>% .$measure %>% max()
  peak_max = df_ %>% dplyr::filter(group == "peak") %>% .$measure %>% max()
  scale_factor = gene_max/peak_max
  df_ %>% 
    mutate(measure = ifelse(group == "peak", measure*scale_factor, measure)) %>% 
    mutate(sd = ifelse(group == "peak", sd*scale_factor, sd)) -> df_
  
  df_ %>% 
#    ggplot(aes(x = ZT, y = measure, color = group, group = group, linetype = group)) + 
    ggplot(aes(x = ZT, y = measure, color = group, group = group)) + 
    geom_line() + 
    scale_linetype_manual(values = c("target_gene" = "solid", "TF" = "solid", "peak" = "dashed")) + 
    geom_ribbon(aes(ymin = measure - sd, ymax = measure + sd, fill = group), alpha = 0.2, color = NA) + 
    scale_y_continuous(name="Mean gene expression", sec.axis=sec_axis(~./scale_factor, name="Mean peak accessibility")) + 
    scale_x_continuous(breaks = seq(2, 22, 4)) + 
    ggtitle(sprintf("TF: %s\nTarget gene: %s\nCRE: %s", TF, target_gene, peak)) -> p
  return(p)
}
PlotGenePeakTF_custom_2 = function(plot){
  plot$data -> dat_
  plot$labels$title -> title_
  
  TF_ = gsub("TF: (.+?)\n.+", "\\1", title_)
  gene_ = gsub("TF: (.+?)\nTarget gene: (.+?)\nCRE: (.+)", "\\2", title_)
  CRE_ = gsub("TF: (.+?)\nTarget gene: (.+?)\nCRE: (.+)\n.+", "\\3", title_)
  validation_ = gsub("TF: (.+?)\nTarget gene: (.+?)\nCRE: (.+)\nValidation: (.+)", "\\4", title_)
#  sprintf("TF=%s, Target gene=%s, CRE=%s", TF_, gene_, CRE_) %>% print()
  
  if (TF_ == "Arntl"){
    TF__ = "Bmal1"
  }else{
    TF__  = TF_
  }
  
#  GSE39977_ChIP_timepoint_mat$multiple[[sprintf("%s_ChIP", TF__)]] %>% 
#  GSE39977_ChIP_timepoint_mat$unique[[sprintf("%s_ChIP", TF__)]] %>% 
  GSE39977_ChIP_timepoint_mat_1[[sprintf("%s_ChIP", TF__)]] %>%
    dplyr::mutate(across(matches("CT|Input"), function(x){(x/sum(x))*1e6})) %>% 
    dplyr::filter(Peaks == CRE_) %>% 
    dplyr::select(-Input) %>% 
    pivot_longer(-Peaks, names_to = "ZT", values_to = "measure") %>% 
    dplyr::mutate(ZT = gsub("CT(\\d+)", "\\1", ZT) %>% as.integer()) %>% 
    dplyr::mutate(sd = 0, group = "ChIP") %>% 
    dplyr::select(group, ZT, measure, sd) -> ChIP_df
#  print(ChIP_df)
  rbind(dat_, ChIP_df) -> dat_
#  print(dat_)
  
  gene_max = dat_ %>% dplyr::filter(group != "peak|ChIP") %>% .$measure %>% max()
  ChIP_max = dat_ %>% dplyr::filter(group == "ChIP") %>% .$measure %>% max()
  scale_factor = gene_max/ChIP_max
  
  dat_ %>% dplyr::mutate(measure = case_when(
    group == "ChIP" ~ measure*scale_factor,
    TRUE ~ measure
  )) -> dat_

  
  dat_$group = factor(dat_$group, levels = c("peak", "target_gene", "TF", "ChIP"))
  
  dat_ %>% 
    ggplot(aes(x = ZT, y = measure, color = group, group = group, linetype = group)) + 
    geom_line() + 
    scale_linetype_manual(values = c("target_gene" = "solid", "TF" = "solid", "peak" = "dashed", "ChIP" = "dashed")) + 
    geom_ribbon(aes(ymin = measure - sd, ymax = measure + sd, fill = group), alpha = 0.2, color = NA) + 
    scale_y_continuous(name="Gene/TF expression", sec.axis=sec_axis(~./scale_factor, name="ATAC peak accessibility/ChIP counts")) + 
    scale_x_continuous(breaks = seq(2, 22, 4)) + 
    scale_color_manual(values = c("peak" = "#F8766D", "TF" = "#619CFF", "target_gene" = "#00BA38", "ChIP" = "#BC2CFF")) + 
    scale_fill_manual(values = c("peak" = "#F8766D", "TF" = "#619CFF", "target_gene" = "#00BA38", "ChIP" = "#BC2CFF")) + 
    ggtitle(sprintf("TF: %s\nTarget gene: %s\nATAC/ChIP peak: %s\nValidation: %s", TF_, gene_, CRE_, validation_)) -> p
  return(p)
}
PlotGenePeakTF_custom_3 = function(gene_, peak_, plot_){
  GSE155161_annotated_CHI_C %>% dplyr::filter(gene == gene_) -> df_
  data.frame(seqnames = gsub("(.+):.+:.+", "\\1", df_$CRE), 
             start = gsub("(.+):(.+):.+", "\\2", df_$CRE) %>% as.integer(), 
             end = gsub("(.+):(.+):(.+)", "\\3", df_$CRE) %>% as.integer(), 
             ZT = df_$ZT, 
             score = df_$score) -> df_
  
  data.frame(seqnames = gsub("(.+)-.+-.+", "\\1", peak_), 
             start = gsub("(.+)-(.+)-.+", "\\2", peak_) %>% as.integer(),
             end = gsub("(.+)-(.+)-(.+)", "\\3", peak_) %>% as.integer()) %>% 
    GenomicRanges::makeGRangesFromDataFrame(.) -> peak_gr
  
  df_$ZT %>% unique() %>% 
    gtools::mixedsort() %>% 
    "names<-"(.,.) %>% 
    map(function(x){
      ZT_ = x
      df_ %>% dplyr::filter(ZT == ZT_) -> df_
      GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> CHI_C_gr
      ChIPpeakAnno::findOverlapsOfPeaks(peak_gr, CHI_C_gr, maxgap = 1000) -> ol
      ol$overlappingPeaks[[1]] -> ol
      ol$score %>% sum() -> interaction_score
      data.frame(group = "CHI_C", ZT = gsub("ZT(\\d+)", "\\1", ZT_) %>% as.integer(), measure = interaction_score, sd = 0)
#      data.frame(group_ = "CHI_C", ZT_ = gsub("ZT(\\d+)", "\\1", ZT_) %>% as.integer(), measure_ = interaction_score)
    }) %>% do.call(rbind, .) -> CHI_C_dat
  
#  scale_factor = max(plot_$data$measure)/max(CHI_C_dat$measure)
  CHI_C_dat$measure = scales::rescale(CHI_C_dat$measure, to = c(min(plot_$data$measure), max(plot_$data$measure)))
  plot_ + 
#    geom_line(data = CHI_C_dat, aes(x = ZT, y = measure, group = group), inherit.aes = F) -> p
    geom_point(data = CHI_C_dat, aes(x = ZT, y = measure, group = group), inherit.aes = F, shape = 8) -> p
  return(p)

#  return(ol)
#  return(peak_gr)
#  return(CHI_C_gr)
}
Plot_ChIP = function(TF, Peak){
  GSE39977_ChIP_timepoint_mat$multiple[[sprintf("%s_ChIP", TF)]] %>% 
    dplyr::mutate(across(matches("CT|Input"), function(x){(x/sum(x))*10^6})) %>% 
    dplyr::filter(Peaks == Peak) %>% 
    pivot_longer(cols = -Peaks, names_to = "ZT", values_to = "measure") %>% 
    dplyr::filter(ZT != "Input") %>%
    dplyr::mutate(ZT = gsub("CT(\\d+)", "\\1", ZT) %>% as.integer()) -> df_
  df_ %>% 
    ggplot(aes(x = ZT, y = measure)) + 
    geom_line() -> p
  return(p)
}
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
      GenomicRanges::makeGRangesFromDataFrame(df_, keep.extra.columns = T) -> CHI_C_gr
      ChIPpeakAnno::findOverlapsOfPeaks(peak_gr, CHI_C_gr, maxgap = 1000) -> ol
      ol$overlappingPeaks[[1]] -> ol
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
combine_ChIP_CHI_C = function(ChIP_plot, CHI_C_plot){
  ChIP_plot$data %>% 
    dplyr::mutate(group = "ChIP") -> ChIP_plot$data
  ChIP_max = max(ChIP_plot$data$measure)
  CHI_C_plot$data %>% 
    dplyr::mutate(group = "CHI_C") -> CHI_C_plot$data
  CHI_C_max = max(CHI_C_plot$data$measure)
  rbind(ChIP_plot$data, CHI_C_plot$data) -> df_
  
  scale_factor = ChIP_max/CHI_C_max
  df_ %>% dplyr::mutate(measure = case_when(
    group == "CHI_C" ~ measure*scale_factor,
    TRUE ~ measure
  )) -> df_
  
#  return(df_)
  ggplot(df_, aes(x = ZT, y = measure, color = group)) + 
    geom_line() -> p
  p + scale_y_continuous(name = "TF ChIP signal", sec.axis = sec_axis(~./scale_factor, name = "CHI-C interaction")) -> p
  return(p)
}

trios_res_df_2

trios_res_df_2 %>% 
  dplyr::filter(model_1 == "TRIPOD", TF %in% Core_clock_genes, gene %in% Core_clock_genes) %>% 
  dplyr::filter((CHI_C_validate == "TRUE")|(ChIP_seq_validated == "TRUE")) %>% 
  dplyr::arrange(TF) %>% 
  {
    df_ = .
    df_[,c("gene", "peak", "TF", "CHI_C_validate", "ChIP_seq_validated")] %>% distinct() -> df_
    df_ %>% dplyr::mutate(validation = case_when(
      (CHI_C_validate == "TRUE")&(ChIP_seq_validated != "TRUE") ~ "CHI_C",
      (ChIP_seq_validated == "TRUE")&(CHI_C_validate != "TRUE") ~ "ChIP-seq",
      (ChIP_seq_validated == "TRUE")&(CHI_C_validate == "TRUE") ~ "Both"
    )) -> df_
    df_$TF %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        TF_ = x
        df_ %>% dplyr::filter(TF == TF_) -> df_
        1:nrow(df_) %>% 
          map(function(i){
            df_[i, ] -> df_
            PlotGenePeakTF_custom(gene = df_$gene, TF = df_$TF, peak = df_$peak) -> p
            PlotGenePeakTF_custom_1(p) -> p
            p + ggtitle(sprintf("TF: %s\nTarget gene: %s\nCRE: %s\nValidation: %s", df_$TF, df_$gene, df_$peak, df_$validation)) -> p
            return(p)
          }) -> p_list
      }) -> p_list
    p_list
  } -> p_list


patchwork::wrap_plots(p_list$Arntl[[1]], p_list$Clock[[1]], p_list$Npas2[[7]], p_list$Rorc[[3]], ncol = 2, guides = "collect") # Fig 4B

#add motif score to the plot
library(Seurat)
library(Signac)
library(tidyverse)

load("~/Downloads/trios_res_df_2.RData")
trios_res_df_2 %>% as.data.frame() -> trios_res_df_2
trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$TF %>% unique() %>% toupper() -> TRIPOD_TF_

load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.rna.rda')
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')

load(file='~/Dropbox/singulomics/github_rda/TRIPOD/sc.rda')
sc@meta.data[,c("ZT", "seurat_clusters")] %>% distinct() %>% dplyr::arrange(seurat_clusters) %>% 
  mutate(metacell = sprintf("metacell_%s", seurat_clusters)) %>% 
  "rownames<-"(., 1:nrow(.)) -> metacell_ZT
metacell_ZT %>% dplyr::mutate(col_ = sprintf("%s_%s", ZT, metacell)) -> metacell_ZT

sc@assays$ATAC@motifs@motif.names %>% unlist() %>% toupper() -> TF_motif
colnames(metacell.rna) %>% toupper() -> RNA_genes
setdiff(TF_motif, RNA_genes)
intersect(TF_motif, RNA_genes) -> TF_
length(TF_)

TF_expr = metacell.rna %>% "colnames<-"(., toupper(colnames(.))) %>% .[,TF_]
rownames(TF_expr) = metacell_ZT$col_
TF_expr %>% t() %>% as.data.frame() -> TF_expr
TF_expr %>% dplyr::mutate(TRIPOD = case_when(rownames(.) %in% TRIPOD_TF_ ~ "yes", TRUE ~ "no")) -> TF_expr
dim(TF_expr)

metacell.peak %>% t() %>% as.data.frame() %>% "colnames<-"(., metacell_ZT$col_) -> CRE_expr
trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$peak %>% unique() -> TRIPOD_CRE_
CRE_expr %>% dplyr::mutate(TRIPOD = case_when(rownames(.) %in% TRIPOD_CRE_ ~ "yes", TRUE ~ "no")) -> CRE_expr

gene_expr = metacell.rna %>% "colnames<-"(., toupper(colnames(.)))
rownames(gene_expr) = metacell_ZT$col_
gene_expr %>% t() %>% as.data.frame() -> gene_expr
trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$gene %>% unique() %>% toupper() -> TRIPOD_gene_
gene_expr %>% dplyr::mutate(TRIPOD = case_when(rownames(.) %in% TRIPOD_gene_ ~ "yes", TRUE ~ "no")) -> gene_expr


sc@meta.data -> sc_meta_
rm(sc);gc()
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc
sc@assays$MOTIF@data -> motif_expr
motif_expr[,rownames(sc_meta_)] -> motif_expr
rownames(motif_expr) -> motif_names
sc@assays$ATAC@motifs@motif.names[motif_names] %>% unlist() %>% unname() %>% toupper() -> rownames(motif_expr)
sc_meta_ %>% 
  rownames_to_column("cells") %>% 
  group_by(seurat_clusters) %>% 
  group_map(function(meta_, cluster_){
    meta_$cells -> cells_
#    motif_expr[,cells_][1:5,1:5] %>% 
    motif_expr[,cells_] %>% 
      rowMeans() %>% 
      as.data.frame() %>% "colnames<-"(., sprintf("cluster_%s", cluster_$seurat_clusters)) -> df_
  }, .keep = T) %>% do.call(cbind, .) -> motif_expr
colnames(motif_expr) <- metacell_ZT$col_
dim(motif_expr)

trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$TF %>% unique() %>% toupper() -> TRIPOD_MOTIF_
motif_expr %>% dplyr::mutate(TRIPOD = case_when(rownames(.) %in% TRIPOD_MOTIF_ ~ "yes", TRUE ~ "no")) -> motif_expr

dim(TF_expr)
dim(gene_expr)
dim(CRE_expr)
dim(motif_expr)

write.csv(TF_expr, file = "~/Dropbox/singulomics/github_rda/TRIPOD/bayesian/TF_expression.csv", col.names = T, row.names = T, quote = F)
write.csv(gene_expr, file = "~/Dropbox/singulomics/github_rda/TRIPOD/bayesian/gene_expression.csv", col.names = T, row.names = T, quote = F)
write.csv(CRE_expr, file = "~/Dropbox/singulomics/github_rda/TRIPOD/bayesian/CRE_expression.csv", col.names = T, row.names = T, quote = F)
write.csv(motif_expr, file = "~/Dropbox/singulomics/github_rda/TRIPOD/bayesian/motif_expression.csv", col.names = T, row.names = T, quote = F)
###

read.csv("~/Dropbox/singulomics/github_rda/TRIPOD/bayesian/motif_expression.csv", row.names = 1, header = T, stringsAsFactors = F, sep = ",") -> motif_expr
list(ARNTL = p_list$Arntl[[1]], CLOCK = p_list$Clock[[1]], 
     RORC = p_list$Rorc[[6]], RORC = p_list$Rorc[[1]], 
     ARNTL = p_list$Arntl[[2]]) %>% 
  map2(.x=.,.y=names(.),.f=function(p_, TF_){
    p_$labels$title -> title_
    print(title_)
    p_$data -> df_
    df_ %>% 
      group_by(group) %>% 
      group_map(function(df_,y){
        max(df_$measure) - min(df_$measure) -> range_original
        df_ %>% 
          dplyr::mutate(measure = scales::rescale(measure, to = c(0,1))) %>% 
          dplyr::mutate(sd = sd/range_original) -> df_
      }, .keep = T) %>% do.call(rbind, .) -> df_
    motif_expr[TF_, ] %>% dplyr::select(-TRIPOD) %>% 
      pivot_longer(everything(), names_to = "group") %>% 
      dplyr::mutate(ZT = gsub("(ZT\\d+?)_.+", "\\1", group)) %>% 
      group_by(ZT) %>% 
      group_map(function(df_1, y){
        df_1$value %>% mean() -> mean_
        df_1$value %>% sd() -> sd_
        ZT_ = y$ZT %>% gsub("ZT(\\d+)", "\\1", .) %>% as.integer()
        data.frame(group = "motif", ZT = ZT_, measure = mean_, sd = sd_)
      }, .keep = T) %>% do.call(rbind, .) -> df_1
    
    max(df_1$measure) - min(df_1$measure) -> range_original
    df_1 %>% 
      dplyr::mutate(measure = scales::rescale(measure, to = c(0,1))) %>% 
      dplyr::mutate(sd = sd/range_original) -> df_1
    
    rbind(df_, df_1) -> df_
    df_ %>% 
      ggplot(aes(x = ZT, y = measure, group = group, color = group, fill = group)) + 
      geom_line() + 
      geom_ribbon(aes(ymin = measure-sd, ymax = measure+sd), alpha = 0.2, color = NA) + 
      scale_x_continuous(breaks = c(2,6,10,14,18,22)) + 
      ylab("relative value (TF expr, motif score, \nCRE accessibility, gene expr)") -> p
    p + ggtitle(title_) -> p
    p + scale_color_manual(values = c(peak = "#f3766e", target_gene = "#2ab34b", TF = "#7094cd", motif = "#a781ba")) -> p
    p + scale_fill_manual(values = c(peak = "#f3766e", target_gene = "#2ab34b", TF = "#7094cd", motif = "#a781ba")) -> p
  }) -> p_list_
patchwork::wrap_plots(p_list_[1:4], ncol = 4, nrow = 1, guides = "collect") #Fig 4B
p_list_[[5]] #Fig 4C


Plot_ChIP(TF = "Bmal1", Peak = "chr1-39102368-39103276") -> p1
Plot_CHI_C(target_gene = "Npas2", Peak = "chr1-39102368-39103276") -> p2
combine_ChIP_CHI_C(ChIP_plot = p1, CHI_C_plot = p2) -> p3 #Fig 4C

library(BSgenome.Mmusculus.UCSC.mm10)
gr = GRanges(seqnames = "chr1", ranges = IRanges(start = 39102368, end = 39103276))
genome <- BSgenome.Mmusculus.UCSC.mm10
dna_sequence <- getSeq(genome, gr)
dna_sequence <- dna_sequence[[1]]
ebox = "CACGTG"
motif_seq <- DNAString(ebox)
motif_pattern <- "CA..TG"
matches <- matchPattern(motif_seq, dna_sequence)
motif_df = as.data.frame(matches)
sequence_df <- data.frame(
  idx = 1:length(dna_sequence),
  position = as.data.frame(gr)$start:as.data.frame(gr)$end,
  nucleotide = as.character(dna_sequence) %>% str_split("") %>% .[[1]]
)
sequence_df = sequence_df[1:50,]
sequence_df$motif <- ifelse(sequence_df$idx %in% motif_df$start:(motif_df$end), "Motif", "Background")

ggplot(sequence_df, aes(x = position, y = 0, label = nucleotide, color = motif)) +
  geom_text(size = 5) + # Add nucleotide letters
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) + 
  ylab("") + 
  xlab("chr1") + 
  scale_x_continuous(breaks =  seq(sequence_df$position[1], sequence_df$position[nrow(sequence_df)], by = 10)) + 
  coord_fixed(ratio = 20) -> p_atac_peak #Fig 4C

#Plot Bmal1 KO vs WT (Npas2 and Cry1)
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

res_pval_KO %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Npas2")
res_pval_WT %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Npas2")

list(WT = Bmal1_WT_df, KO = Bmal1_KO_df) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    df_ = x
    df_ %>% dplyr::filter(Gene == "Npas2") %>% 
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
  geom_ribbon(aes(ymin = expression - sd, ymax = expression + sd), alpha = 0.3, color = NA) #Fig 4C

res_pval_KO %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Per1")
res_pval_WT %>% dplyr::mutate(cauchy_p = sprintf("%.2e", cauchy_p)) %>% dplyr::filter(Gene == "Per1")

list(WT = Bmal1_WT_df, KO = Bmal1_KO_df) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    group_ = y
    df_ = x
    df_ %>% dplyr::filter(Gene == "Cry1") %>% 
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
  geom_ribbon(aes(ymin = expression - sd, ymax = expression + sd), alpha = 0.3, color = NA) #Fig 4C

trios_res_df_2 %>% 
  dplyr::filter(model_1 == "TRIPOD", gene_cauchy_BH.Q < 0.05, TF_cauchy_BH.Q < 0.05, peak_cauchy_BH.Q < 0.05) %>% 
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

#Plot Arntl and Nr1d1 regulation
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
      geom_abline(color = "red", linetype = 2) + 
      geom_vline(xintercept = TF_phase, color = "blue", linetype = 2) + 
      geom_ribbon(aes(xmin = TF_phase-2.5, xmax = TF_phase+2.5), fill = "blue", alpha = 0.2) + 
      geom_ribbon(aes(xmin = (TF_phase-12)-2.5, xmax = (TF_phase-12)+2.5), fill = "green", alpha = 0.2) + 
      theme_classic() + 
      ggtitle(sprintf("r=%s", cor_))
  } #Fig 4D

trios_res_df_3 %>% dplyr::filter(TF == "Nr1d1") %>% 
  dplyr::filter(gene %in% Core_clock_genes) %>% View()

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
      geom_abline(color = "red", linetype = 2) + 
      geom_vline(xintercept = TF_phase, color = "blue", linetype = 2) + 
      geom_ribbon(aes(xmin = TF_phase-2.5, xmax = TF_phase+2.5), fill = "blue", alpha = 0.2) + 
      geom_ribbon(aes(xmin = (TF_phase+12)-2.5, xmax = (TF_phase+12)+2.5), fill = "green", alpha = 0.2) + 
      theme_classic() + 
      ggtitle(sprintf("r=%s", cor_))
  } #Fig 4D

#Plot Anrtl's TFs
Arntl_TF = c("Mef2c", "Klf1", "Nr5a2", "Mef2d", "Foxa3", "Esrra", "Hnf4a", "Mef2a", "Nfyb", "Ppard")
trios_res_df_2 %>% 
  dplyr::filter(model_1 == "TRIPOD", TF %in% Arntl_TF, gene == "Arntl") %>% 
  dplyr::filter((CHI_C_validate == "TRUE")) %>% 
  dplyr::arrange(TF) %>% 
  {
    df_ = .
    df_[,c("gene", "peak", "TF", "CHI_C_validate", "ChIP_seq_validated")] %>% distinct() -> df_
    df_ %>% dplyr::mutate(validation = case_when(
      (CHI_C_validate == "TRUE")&(ChIP_seq_validated != "TRUE") ~ "CHI_C",
      (ChIP_seq_validated == "TRUE")&(CHI_C_validate != "TRUE") ~ "ChIP-seq",
      (ChIP_seq_validated == "TRUE")&(CHI_C_validate == "TRUE") ~ "Both"
    )) -> df_
    df_$TF %>% 
      unique() %>% 
      "names<-"(.,.) %>% 
      map(function(x){
        TF_ = x
        df_ %>% dplyr::filter(TF == TF_) -> df_
        1:nrow(df_) %>% 
          map(function(i){
            df_[i, ] -> df_
            PlotGenePeakTF_custom(gene = df_$gene, TF = df_$TF, peak = df_$peak) -> p
            PlotGenePeakTF_custom_1(p) -> p
            p + ggtitle(sprintf("TF: %s\nTarget gene: %s\nCRE: %s\nValidation: %s", df_$TF, df_$gene, df_$peak, df_$validation)) -> p
            return(p)
          }) -> p_list
      }) -> p_list
    p_list
  } -> p_list #Supp Fig 22

p_list$Mef2a[[1]]
p_list$Mef2d[[1]]
p_list$Ppard[[2]]
p_list$Nr5a2[[1]]
p_list$Nfyb[[1]]

#Plot regulatory network of clock genes
Core_clock_genes = c("Dbp", "Arntl", "Bhlhe40", "Bhlhe41", "Nfil3", "Rorc", "Rora", "Nr1d1", "Clock", "Npas2", "Cry1", "Ciart", "Per1", "Per2")
TF_clock_genes = c("Dbp", "Arntl", "Bhlhe40", "Bhlhe41", "Nfil3", "Rora", "Rorc", "Nr1d1", "Clock", "Npas2")

trios_res_df_1 %>% 
  dplyr::filter(model_1 == "TRIPOD", gene %in% Core_clock_genes, 
                TF %in% Core_clock_genes,
                #                peak_cauchy_BH.Q < 0.05, 
                gene_cauchy_BH.Q < 0.05, 
                TF_cauchy_BH.Q < 0.05
  ) -> df_
head(df_)
dim(df_)

df_ %>% 
  group_by(gene, peak, TF) %>% 
  group_map(function(x, y){
    x$peak_biotype %>% unique() %>% paste(., collapse = ";") -> peak_biotype_
    x$CHI_C_validate %>% unique() %>% paste(., collapse = ";") -> CHI_C_validate_
    data.frame(gene = y$gene, peak = y$peak, TF = y$TF, peak_biotype = peak_biotype_, CHI_C_validate = CHI_C_validate_)
  }, .keep = T) %>% 
  do.call(rbind, .) -> df_
dim(df_)
View(df_)

Core_clock_genes %>% 
  as.data.frame() %>% 
  "colnames<-"(., "name") %>% 
  mutate(category = case_when(
    name %in% df_$TF ~ "TF/gene",
    TRUE ~ "gene"
  )) -> nodes
nodes %>% filter(name %in% c(df_$gene, df_$TF)) -> nodes
df_ %>% 
  dplyr::select(TF, gene, peak, peak_biotype, CHI_C_validate) %>% 
  "colnames<-"(.,c("from", "to", "by", "peak_biotype", "CHI_C_validate")) -> edges
edges %>% dplyr::mutate(CHI_C_validate = as.character(CHI_C_validate)) -> edges

g <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
#category_colors <- c("eRNA" = "red", "No annotation" = "blue")
#igraph::E(g)$color <- category_colors[edges$category]
igraph::E(g)$color <- "blue"
node_colors <- c("TF/gene" = "lightblue", "gene" = "orange")
igraph::V(g)$color <- node_colors[nodes$category]
line_types <- c("TRUE" = "solid", "FALSE" = "dashed")
igraph::E(g)$lty <- line_types[edges$CHI_C_validate]
#igraph::E(g)$width <- edges$weight

layout_kk <- igraph::layout_with_kk(g)
layout_circle <- igraph::layout.circle(g)

plot(g,
     #     layout = layout_kk,
     layout = layout_circle,
     edge.arrow.size = 0.4,        # Size of the arrow heads
     edge.curved = 0.3,            # Curvature of the edges
     vertex.size = 24,             # Size of the vertices
     vertex.label.color = "black", # Color of the vertex labels
     vertex.label.cex = 1,       # Font size of the vertex labels
     #     edge.width = igraph::E(g)$width,
     edge.lty = igraph::E(g)$lty,          # Linetypes of the edges
) #Supp Fig 23

#Plot Chromatin States validation
load("~/Dropbox/singulomics/github_rda/TRIPOD/ChromHMM_dat.RData")
load(file='~/Dropbox/singulomics/github_rda/TRIPOD/metacell.peak.rda')

c("chromatin_state_15", "chromatin_state_18") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    trios_res_df_2 %>% 
      dplyr::filter(model_1 == "TRIPOD") -> df_
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

trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$peak %>% unique() %>% length() -> N_TRIPOD_unique_peaks
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

write.csv(p18_state_df, file = "~/Downloads/chromatin_state_df.csv", row.names = F, quote = F)
read.csv("~/Downloads/chromatin_state_df.csv", header = T, stringsAsFactors = F) -> chromatin

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
  theme_minimal() + xlab("chromHMM states")+ scale_fill_manual(values=c('#999999','#E69F00')) #Fig 4E

# Plot CHI-C validation
load(file = "~/Dropbox/singulomics/github_rda/TRIPOD/GSE155161_CHI-C_annotated.rda")

intersect(
trios_res_df_2 %>% filter(model_1 == "TRIPOD") %>% .$gene %>% unique(),
GSE155161_annotated_CHI_C$gene %>% unique()
) -> intersected_genes
length(intersected_genes)

c(
  trios_res_df_2 %>% filter(model_1 == "TRIPOD") %>% .$gene %>% unique(),
  GSE155161_annotated_CHI_C$gene %>% unique()
) %>% unique() %>% 
  "names<-"(., .) %>% 
#  .[1] %>% 
  map(function(x){
    gene_ = x
    if (gene_ %in% (trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD") %>% .$gene)){
      trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD" & gene == gene_) %>% dplyr::select(gene, peak) -> trios_res_df_2
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
    
    trios_res_df_2 %>% filter(model_1 == "TRIPOD" & gene == gene_) -> df_trios
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
    
    ChIPpeakAnno::findOverlapsOfPeaks(df_trios, df_CHI_C, maxgap = 1000) -> ol
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

save(all_gene_peaks_linkage, gene_peak_linkage_hypergeometric_test_df, file = "~/Dropbox/singulomics/github_rda/TRIPOD/hypergeometic_test_dat.rda")

#Boxplot (hypergeometric test) -> TRIPOD vs CHI-C ----
c("p_val", "neg_log10_p_val") %>% 
  "names<-"(.,.) %>% 
  map(function(x){
    gene_peak_linkage_hypergeometric_test_df %>% 
      mutate(p_value = case_when(
        p_value == 0 ~ .Machine$double.eps, 
        TRUE ~ p_value
      )) %>% mutate(group = "TRIPOD vs CHI-C") %>% 
      dplyr::filter(p_value < 0.05) -> df_
    
    n_gene = df_$gene %>% unique() %>% length()
    n_sign = df_ %>% filter(p_value < 0.05) %>% nrow()
    median_ = df_$p_value %>% median() %>% sprintf("%.2e", .)
    
    if (x == "p_val"){
      df_ %>% 
        ggplot(aes(x = group, y = p_value, fill = group)) -> p
    }else{
      df_ %>% 
        ggplot(aes(x = group, y = -log10(p_value), fill = group)) -> p
    }

    p + 
      geom_boxplot(width = 0.5) + 
      theme_classic() + 
      xlab(NULL)  -> p
    
    if (x == "p_val"){
      p + ylab("Hypergeometric test p-value") -> p
    }else{
      p + ylab("Hypergeometric test -log10(p-value)") -> p
    }
    
    p + theme(legend.position = "none") +
      ggtitle(sprintf("%s sig. genes\nP-value median: %s", n_sign, median_)) -> p
  }) -> p_list

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

trios_res_df_2 %>% 
  {
    df_ = .
    set.seed(123)
    sample(1:nrow(df_), nrow(df_), replace = F) -> se
    df_[["peak"]] = df_[["peak"]][se]
    df_
  } -> trios_res_df_2_shuffled

intersected_genes %>% 
  "names<-"(.,.) %>% 
  #  .[1] %>% 
  #  .["Tmem14a"] %>% 
  future_map(function(x){
    gene_ = x
    print(gene_)
    
    trios_res_df_2_shuffled %>% filter(model_1 == "TRIPOD" & gene == gene_) -> df_trios
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
    
    ChIPpeakAnno::findOverlapsOfPeaks(df_trios, df_CHI_C, maxgap = 1000) -> ol
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
dim(gene_peak_linkage_hypergeometric_test_df_shuffled)
gene_peak_linkage_hypergeometric_test_df_shuffled$p_value %>% hist() %>% plot()

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

library(ggbreak)
rbind(
p$data %>% mutate(group = "TRIPOD predicted peaks") %>% dplyr::select(p_value, group),
p1$data %>% mutate(group = "random shuffled peaks") %>% dplyr::select(p_value, group)
) %>% 
  {
    dat_ = .
    ggplot() + 
#      geom_histogram(data=subset(dat_, group == "TRIPOD predicted peaks"), aes(x = p_value, y = ..count..), binwidth = 0.5, color = "#619CFF", fill = "#619CFF") + 
#      geom_histogram(data=subset(dat_, group == "random shuffled peaks"), aes(x = p_value, y = -..count..), binwidth = 0.5, color = "#00BA38", fill = "#00BA38") + 
      geom_histogram(data=subset(dat_, group == "TRIPOD predicted peaks"), aes(x = p_value, y = ..count.., fill = "TRIPOD peaks"), binwidth = 0.5) + 
      geom_histogram(data=subset(dat_, group == "random shuffled peaks"), aes(x = p_value, y = -..count.., fill = "Random peaks"), binwidth = 0.5) + 
      scale_y_continuous(breaks = seq(-4000, 300, 100)) + 
      scale_x_continuous(limits = c(-1, 100)) + 
      scale_y_break(c(-3800, -200), scales = 10) + 
      scale_fill_manual(values = c("TRIPOD peaks" = "#619CFF", "Random peaks" = "#00BA38"), breaks = c("TRIPOD peaks", "Random peaks")) +
      labs(fill = "Group") + 
      theme_classic() + 
      theme(legend.position = "top") + 
      xlab("Hypergeometric test -log10(p-value)") + 
      ggtitle("4241 sig. genes")
  }

p$data %>% head()
p1$data %>% head()

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

write.csv(df_hypergeometric_test, file = "~/Downloads/hypergeometric_test_df.csv", row.names = FALSE, quote = FALSE)

pvals=read.csv('~/Downloads/hypergeometric_test_df.csv')
tripod.p=10^(-pvals$X.log10_P[pvals$group=='TRIPOD'])
random.p=10^(-pvals$X.log10_P[pvals$group=='random'])

tripod.p=pmin(70, (pvals$X.log10_P[pvals$group=='TRIPOD']))
random.p=(pvals$X.log10_P[pvals$group=='random'])

pdf(file='pval.pdf', width=5, height=4) #Fig 5E
hist(tripod.p, breaks=seq(0,70,2), col=adjustcolor("#E69F00", alpha.f=0.7), xlab='-log(p)', main='Distribution of -log(p) \nfrom hypergeometric test')
hist(random.p, add=TRUE, breaks=seq(0,70,2), col=adjustcolor("#999999", alpha.f = 0.7))
legend('topright', fill=c(adjustcolor("#999999", alpha.f=0.7), adjustcolor("#E69F00", alpha.f = 0.7)),
       legend = c('Random peaks', 'TRIPOD peaks'), bty='n')
dev.off()


hist(random.p,  breaks=seq(0,70,2), col=adjustcolor("#999999", alpha.f = 0.7))

#Plot ChIP-seq validation

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
        
        trios_res_df_2 %>% dplyr::filter(model_1 == "TRIPOD" & TF == TF_) %>% 
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

venn_list_df_ %>% filter(macs2_threshold == "q_0.1" & max_gap == 1000) %>% 
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
  ylab("Fraction of peaks overlapping with ChIP-seq peaks") #Fig 5E