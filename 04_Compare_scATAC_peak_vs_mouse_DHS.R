library(tidyverse)

# 1. Read Data ----
readRDS("~/Dropbox/singulomics/github_rda/integrated_sc_all_cell_types.rds") -> sc

# 2. Find overlapping peaks between scATAC and DHS peaks ----
options(bedtools.path = "/Users/chunyiptong/anaconda3/envs/bedtools_env/bin/bedtools")
library(bedtoolsr)
read.table("~/Dropbox/singulomics/jerome_chip_seq/mm10peaks_DHS_EncodeBenLiver_67443peaks.bed", header = F, sep = "\t", 
           stringsAsFactors = F) -> jerome_dhs
jerome_dhs %>% dplyr::arrange(V1,V2,V3) -> jerome_dhs
jerome_dhs[,1:3] -> jerome_dhs
colnames(jerome_dhs) = c("chrom", "chromStart", "chromEnd")
jerome_dhs = GenomicRanges::makeGRangesFromDataFrame(jerome_dhs)

sc@assays$ATAC@data %>% rownames() %>% 
  #  .[1:10] %>%
  map(function(x){
    strsplit(x, "-") %>% as.data.frame() %>% 
      t() %>% as.data.frame() %>% 
      "colnames<-"(., c("V1", "V2", "V3")) %>% 
      mutate(across(matches("2|3"), as.numeric))
  }) %>% do.call(rbind, .) %>% 
  "rownames<-"(., 1:nrow(.)) -> atac_peak
atac_peak %>% dplyr::arrange(V1, V2, V3) -> atac_peak
colnames(atac_peak) = c("chrom", "chromStart", "chromEnd")
atac_peak = GenomicRanges::makeGRangesFromDataFrame(atac_peak)

library(ChIPpeakAnno)
ChIPpeakAnno::findOverlapsOfPeaks(atac_peak, jerome_dhs, connectedPeaks = "keepAll") -> peak_overlap
#####

#Supp_Fig_6A ----
venn <- makeVennDiagram(peak_overlap, totalTest = 200000,
                        fill = c("#009E73", "#F0E442"),
                        col = c("#D55E00", "#0072B2"),
                        cat.col = c("#D55E00", "#0072B2"), 
                        connectedPeaks = "keepAll")

# 3. Calculate scATAC peaks intensitiy (ATAC unique, Mouse DHS unique, Overlapped and Random region) ----
library(Seurat)
library(Signac)

fpath = "/Users/chunyiptong/Dropbox/singulomics/aggregate_analysis/atac_fragments.tsv.gz"
cells_ = colnames(sc)
fragments_ <- CreateFragmentObject(fpath, cells = cells_)

peak_overlap$peaklist$`atac_peak///jerome_dhs` -> overlapped_peaks
overlapped_peaks_mat <- Signac::FeatureMatrix(fragments = fragments_, 
                                        features = overlapped_peaks, 
                                        cells = cells_)
dim(overlapped_peaks_mat)

peak_overlap$peaklist$jerome_dhs -> mouse_dhs
mouse_dhs_mat <- Signac::FeatureMatrix(fragments = fragments_, 
                                       features = mouse_dhs, 
                                       cells = cells_)
dim(mouse_dhs_mat)

peak_overlap$peaklist$atac_peak -> atac_peaks
atac_peaks_mat <- Signac::FeatureMatrix(fragments = fragments_, 
                                       features = atac_peaks, 
                                       cells = cells_)
dim(atac_peaks_mat)

#Make random Peak matrix
sc@assays$ATAC@counts %>% rownames() %>% 
  {list(peak = .)} %>% 
  map(function(x){
    gsub(".+?-(.+?)-.+", "\\1", x) %>% as.numeric() -> start
    gsub(".+?-.+?-(.+)", "\\1", x) %>% as.numeric() -> end
    len = end-start
    len = median(len)
    return(len)
  }) %>% .[[1]] -> peak_median_len

list(
  chr1=195471971,
  chr2=182113224,
  chr3=160039680,
  chr4=156508116,
  chr5=151834684,
  chr6=149736546,
  chr7=145441459,
  chr8=129401213,
  chr9=124595110,
  chr10=130694993,
  chr11=122082543,
  chr12=120129022,
  chr13=120421639,
  chr14=124902244,
  chr15=104043685,
  chr16=98207768,
  chr17=94987271,
  chr18=90702639,
  chr19=61431566
#  chrX=171031299,
#  chrY=91744698,
#  chrM=16299
) -> mm10_genome_size

mm10_genome_size %>% 
  map2(.x=.,.y=names(.),function(x,y){
    chr_ = y
    len = x - peak_median_len
    start = sample(1:len, round(10000/19), replace = F)
    end = start + peak_median_len
    data.frame(chrom=chr_, chromStart=start, chromEnd=end) -> df_
    df_ %>% dplyr::arrange(chromStart, chromEnd) -> df_
#    head(df_)
  }) %>% do.call(rbind, .) %>% 
  GenomicRanges::makeGRangesFromDataFrame(.) -> random_peaks

random_peaks_mat <- Signac::FeatureMatrix(fragments = fragments_, 
                                        features = random_peaks, 
                                        cells = cells_)

#Matrix normalization

list(ATAC = atac_peaks_mat, 
     Overlapped = overlapped_peaks_mat, 
     Mouse_DHS = mouse_dhs_mat, 
     random_peaks = random_peaks_mat) %>% 
  do.call(rbind, .) -> total_mat
total_mat %>% as.data.frame() %>% 
  mutate(across(everything(), function(x){(x/sum(x))*10^6})) -> total_mat

total_mat %>% rownames() %>% 
  list(peaks = .) %>% 
  map(function(x){
    gsub(".+?-(.+?)-.+", "\\1", x) %>% as.numeric() -> start
    gsub(".+?-.+?-(.+)", "\\1", x) %>% as.numeric() -> end
    len = end - start
    return(len)
  }) %>% .[[1]] -> len

dim(total_mat)
#[1] 164117  26068

list(ATAC=nrow(atac_peaks_mat), 
     Overlapped=nrow(overlapped_peaks_mat), 
     Mouse_DHS=nrow(mouse_dhs_mat), 
     Random_peaks = nrow(random_peaks_mat)) %>% 
  purrr::reduce(., c) %>% 
  cumsum() %>% {c(1,.)} -> total_mat_idx

rm(atac_peaks_mat, mouse_dhs_mat, overlapped_peaks_mat)
gc()

i = 0 
list(ATAC=1, Overlapped=2, Mouse_DHS=3, Random_peaks = 4) %>% 
  map2(.x=.,.y=names(.),.f=function(x,y){
    
    print(y)
    i + 1 ->> i
    
    if (i <= 4){
      total_mat_idx[x] -> start
      total_mat_idx[x+1] -> end      
    }
    
    if (i %in% 2:4){
      start = start+1
    }
    
    total_mat[start:end, ] -> df_
    
    df_ %>% rowMeans() %>% 
      as.data.frame() %>% 
      "colnames<-"(., "Mean") -> df_
    df_ %>% mutate(Group = y, .before=1) %>% 
      rownames_to_column("Peaks") -> df_
  }) %>% do.call(rbind, .) -> df_tmp

df_tmp %>% mutate(peak_len = (as.numeric(gsub(".+?-.+?-(.+)","\\1", Peaks))-as.numeric(gsub(".+?-(.+?)-.+","\\1", Peaks)))) %>% 
  mutate(Mean_by_len = (Mean/peak_len)*1000) -> df_tmp

##Add Histone Marks info
##Check histone mark ----
library(rtracklayer)

list.files("~/Dropbox/singulomics/github_rda/linkpeak/encode_mm10_histone_mark", pattern = "\\.bigwig", full.names = T) %>% 
  .[c(1,3)] %>% 
  "names<-"(., gsub("/.+/Encode_(.+?)_.+", "\\1", .)) %>% 
  map2(.x=.,.y=names(.),.f=function(bigwig_, bw_name_){
    print(bigwig_)
    bw <- import.bw(bigwig_)
    data.frame(
      seqnames = gsub("(.+)-(.+)-(.+)", "\\1", df_tmp$Peaks), 
      start = gsub("(.+)-(.+)-(.+)", "\\2", df_tmp$Peaks) %>% as.numeric(), 
      end = gsub("(.+)-(.+)-(.+)", "\\3", df_tmp$Peaks) %>% as.numeric()
    ) %>% GenomicRanges::makeGRangesFromDataFrame() -> peak_gr
    GenomicRanges::findOverlaps(peak_gr, bw) -> overlaps
    
    bw[subjectHits(overlaps)]@ranges@width -> bw_width
    peak_gr[queryHits(overlaps)]@ranges@width -> peak_width
    width_ <- pintersect(peak_gr[queryHits(overlaps)], bw[subjectHits(overlaps)])@ranges@width
    
    signal_values <- mcols(bw)$score[subjectHits(overlaps)]
    peaks_signal <- data.frame(
      peak_id = queryHits(overlaps),
      signal = signal_values,
      overlapping_width = width_,
      bw_width = bw_width, 
      peak_width = peak_width
    )
    
    mean_signal <- peaks_signal %>%
#      group_by(peak_id) %>%
      group_by(peak_id, peak_width) %>% 
      group_map(function(x,y){
        df_ = x
        df_ %>% dplyr::mutate(norm_signal = (signal*overlapping_width)) -> df_
        normalized_signal = sum(df_$norm_signal)/df_$peak_width[1]
        data.frame(peak_id = y$peak_id, mean_signal = normalized_signal)
      }, .keep = T) %>% do.call(rbind, .) %>% 
      "colnames<-"(.,c("peak_id" ,sprintf("%s_signal", bw_name_)))
  }) %>% purrr::reduce(., full_join, by = "peak_id") -> histone_mark_df
  
df_tmp %>% 
  dplyr::mutate(peak_id = 1:nrow(.)) %>% 
  dplyr::full_join(x = ., y = histone_mark_df, by = "peak_id") %>% 
  mutate(H3K27ac_signal = ifelse(is.na(H3K27ac_signal), 0, H3K27ac_signal)) %>% 
  mutate(H3K4me1_signal = ifelse(is.na(H3K4me1_signal), 0, H3K4me1_signal)) %>% 
  select(-peak_id) -> df_tmp

library(ggpubr)
library(ggbreak)

pairwise.t.test(df_tmp$H3K27ac_signal, df_tmp$Group, paired = F, alternative = "two.sided", pool.sd = F)
pairwise.wilcox.test(df_tmp$H3K27ac_signal, df_tmp$Group, paired = F, alternative = "two.sided", pool.sd = F)
pairwise.t.test(df_tmp$H3K4me1_signal, df_tmp$Group, paired = F, alternative = "two.sided", pool.sd = F)
pairwise.wilcox.test(df_tmp$H3K4me1_signal, df_tmp$Group, paired = F, alternative = "two.sided", pool.sd = F)

df_tmp %>% 
  ggplot(aes(x = Group, y = H3K27ac_signal, fill = Group)) + 
  geom_boxplot() + 
  scale_y_break(breaks = c(7,25), scales = 0.5) + 
  theme_classic() -> p1 #Supp_Fig_6B ----

df_tmp %>% 
  ggplot(aes(x = Group, y = H3K4me1_signal, fill = Group)) + 
  geom_boxplot() +
  scale_y_break(breaks = c(3,6), scales = 0.5) + 
  theme_classic() -> p2 #Supp_Fig_6B ----


my_comparisons = list(c("Overlapped", "Mouse_DHS"), 
                     c("Overlapped", "ATAC"), 
                     c("Overlapped", "Random_peaks"), 
                     c("ATAC", "Mouse_DHS"), 
                     c("ATAC", "Random_peaks"), 
                     c("Mouse_DHS", "Random_peaks"))
df_tmp %>% 
  mutate(Mean = sqrt(Mean)) %>%
  ggplot(aes(x = Group, y = Mean, fill = Group)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,20)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(13, 14, 15, 16, 17, 18)) + 
  ggpubr::stat_compare_means(method = "anova", label.y = c(20)) + 
  theme_classic() -> p3

df_tmp %>% 
  mutate(Mean_by_len = sqrt(Mean_by_len)) %>%
  ggplot(aes(x = Group, y = Mean_by_len, fill = Group)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,20)) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, label.y = c(13, 14, 15, 16, 17, 18)) + 
  ggpubr::stat_compare_means(method = "anova", label.y = c(20)) + 
  theme_classic() -> p4

df_tmp %>% 
  ggplot(aes(x = peak_len, color = Group, fill = Group)) + 
  geom_histogram(binwidth = 20) + 
  coord_cartesian(xlim = c(0,1100)) + 
  theme_classic() -> p5

#####

# 4. Draw coverage plot of the Mouse_DHS peaks and scATAC peaks (overlapping) -----

c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Dbp", "Nfil3", "Nr1d1", "Rora", "Rorc", 
  "Cry1", "Ciart", "Per1", "Per2") -> se_genes
se_genes %>% "names<-"(., .) %>% 
  map(function(x){
    CoveragePlot(
      object = sc,
      region = x,
      features = x,
      expression.assay = "SCT",
      group.by = "ZT",
      extend.upstream = 5000,
      extend.downstream = 5000, 
      ranges = jerome_dhs
    ) -> p
  }) -> p_list

#ggpubr::ggarrange(plotlist = p_list, ncol = 3, nrow = 5)
ggpubr::ggarrange(p_list$Arntl, p_list$Nr1d1, p_list$Dbp, ncol = 3) #Supp_Fig_6C ----