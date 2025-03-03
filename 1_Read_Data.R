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

# The 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("aggregate_analysis/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
rm(inputdata.10x)

dim(rna_counts)
dim(atac_counts)

# Create Seurat RNA object
sc <- CreateSeuratObject(counts = rna_counts)
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
rm(rna_counts)

# We'll only use peaks in standard chromosomes
peak.ref <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
seqlevelsStyle(peak.ref)='UCSC'
grange.use <- seqnames(peak.ref) %in% c(c(standardChromosomes(peak.ref),'chrM')) # chr1-19, chrX, chrY, chrM
atac_counts <- atac_counts[as.vector(grange.use), ]
peak.ref <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
peak.ref
rm(grange.use)

# Gene annotation for mm10
gene.ref <- genes(EnsDb.Mmusculus.v79)
genome(gene.ref)='mm10'
seqlevels(gene.ref) <- str_replace(string=paste("chr",seqlevels(gene.ref),sep=""), pattern="chrMT", replacement="chrM") # change to UCSC style
gene.ref=gene.ref[seqnames(gene.ref) %in% c(c(standardChromosomes(peak.ref),'chrM')),]

sc = sc[rownames(sc) %in% gene.ref$gene_name,]
gene.ref=gene.ref[match(rownames(sc), gene.ref$gene_name)]

dim(sc); length(gene.ref)
dim(atac_counts); length(peak.ref)

# Infer cell groups
temp=unlist(strsplit(colnames(sc),'-'))
temp=as.numeric(temp[seq(2, length(temp),2)])
table(temp)
group=factor(c('WT6','KO1','WT3','WT4','WT1','WT5','WT2','KO2')[temp])
ZT=factor(c('ZT22','ZT06','ZT10','ZT14','ZT02','ZT18','ZT06','ZT18')[temp])
rm(temp)
meta.data=data.frame(group, ZT)
rownames(meta.data)=colnames(sc)
sc=AddMetaData(sc, meta.data)
rm(group, meta.data, ZT)

# Now add in the ATAC-seq data
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
rm(chrom_assay); rm(atac_counts); rm(annotations)

table(sc$group)

DefaultAssay(sc)='ATAC'

sc <- NucleosomeSignal(sc)
sc <- TSSEnrichment(sc, fast = F)

sc$high.tss <- ifelse(sc$TSS.enrichment > 2, 'High', 'Low')
sc$nucleosome_group <- ifelse(sc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

save.image(file='rda/1_Read_Data.rda')


