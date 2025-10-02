#Fig_S27
library(tidyverse)
library(bedtoolsr)
options(bedtools.path = "/Users/chunyiptong/anaconda3/envs/bedtools_env/bin")
setwd("~/Dropbox/singulomics")
load("./rda/3_transferred.rda")

gene.ref %>% width() %>% density() %>% plot()
gene.ref %>% width() %>% summary()
gene.ref %>% width() %>% {. > 2000} %>% sum()

gene.ref[width(gene.ref) > 2000, ] -> gene.ref_2000bp

upstream = 1000 
downstream = 1000
gene.tss.ref = promoters(gene.ref_2000bp, upstream = upstream, downstream = downstream) 
#subset gene.tss.ref to chr1 to chr19, chrM, chrX, chrY
gene.tss.ref = gene.tss.ref[seqnames(gene.tss.ref) %in% paste0("chr", c(1:19, "M", "X", "Y")),]
#convert gene.tss.ref to bed file format dataframe
gene.tss.ref.gtf = as.data.frame(gene.tss.ref) %>% 
  select(seqnames, start, end, gene_name) %>%
  mutate(source = "R", .after="seqnames") %>% 
  mutate(feature = "peaks", .after="source") %>% 
  mutate(score = ".", strand = ".", frame = ".", attribute = sprintf('gene_name "%s"', gene_name)) %>% 
  select(-gene_name)

#make a non-tss gene reference from gene.ref
1:length(gene.tss.ref) %>% 
#  .[1:1000] %>% 
  map(function(x){
    print(x)
    GenomicRanges::setdiff(gene.ref_2000bp[x, ], gene.tss.ref[x,], ignore.strand=TRUE)
  }) %>% purrr::reduce(., c) -> gene.non.tss.ref
#mcols(gene.non.tss.ref) <- mcols(gene.tss.ref[1:1000, ])
mcols(gene.non.tss.ref) <- mcols(gene.tss.ref)
gene.non.tss.ref

gene.non.tss.ref.gtf = as.data.frame(gene.non.tss.ref) %>% 
  select(seqnames, start, end, gene_name) %>%
  mutate(source = "R", .after="seqnames") %>% 
  mutate(feature = "peaks", .after="source") %>% 
  mutate(score = ".", strand = ".", frame = ".", attribute = sprintf('gene_name "%s"', gene_name)) %>% 
  select(-gene_name)

write.table(x = gene.tss.ref.gtf,
            file = "~/Dropbox/singulomics/rda/gene.tss.ref.gtf",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(x = gene.non.tss.ref.gtf,
            file = "~/Dropbox/singulomics/rda/gene.non.tss.ref.gtf",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#GSE35790 ----
##PolII ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Polr2b.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Polr2b.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

df_ %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  ggtitle("GSE35790 TSS counts/non TSS counts") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_ratio
p_tss_ratio_PolII_GSE35790 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("non TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_non_tss_counts

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)

list_GSE35790 = list()
list_GSE35790[["Polr2b"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)


#fig 6X12

##H3K36me3 ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("H3K36me3.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("H3K36me3.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

#df_ %>% 
list_GSE35790$H3K36me3$tss_ratio %>% 
  .[,-7] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("GSE35790 TSS counts/non TSS counts") -> p_tss_ratio
p_tss_ratio_H3K36me3_GSE35790 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("TSS counts normalized by lib_size") -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("non TSS counts normalized by lib_size") -> p_non_tss_counts

list_GSE35790[["H3K36me3"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)
#fig 6X12

##H3K4me3 ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("H3K4me3.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("H3K4me3.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

#df_ %>% 
list_GSE35790$H3K4me3$tss_ratio %>% 
  .[,-7] %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("GSE35790 TSS counts/non TSS counts") -> p_tss_ratio
p_tss_ratio_H3K4me3_GSE35790 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("TSS counts normalized by lib_size") -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle("non TSS counts normalized by lib_size") -> p_non_tss_counts

list_GSE35790[["H3K4me3"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)
#fig 6X12

##macs2 peakcalling ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/macs2", 
           pattern = "\\.broadPeak$", full.names = TRUE) %>% 
  list() %>% 
  map(function(x){
   group_ = gsub("/.+/(.+?)_ZT.+", "\\1", x) 
   data.frame(files_ = x, group_ = group_)
  }) %>% .[[1]] %>% group_by(group_) %>% 
  group_map(function(x,y){
    cbind(x,y) -> df_
    df_$files_ %>% 
      map(function(x_){
        read.table(x_, sep = "\t", header = FALSE, stringsAsFactors = F) -> df_
        df_ %>% filter(V9 > 2) %>% dplyr::select(V1, V2, V3) -> df_
      }) %>% do.call(rbind, .) %>% arrange(V1, V2) -> df_
    bedtoolsr::bt.merge(df_, d = 1000) -> df_
    assign(sprintf("GSE35790_%s_peaks", y$group_), value = df_, envir = .GlobalEnv)
    attr(df_, "group") <- y$group_
    return(df_)
  }) -> list_
map(list_, function(x){attributes(x)$group}) -> names_
names(list_) <- names_

list_$Polr2b %>% filter(!grepl("_", V1)) %>% .$V1 %>% table()

list_ %>% map(function(x){
  gene.ref = gene.ref[seqnames(gene.ref) %in% paste0("chr", c(1:19, "M", "X", "Y")),]
  x %>% filter(!grepl("_", V1)) -> x
  GRanges(seqnames = x$V1, ranges = IRanges(start = x$V2, end = x$V3)) -> gr_
  findOverlaps(query = gr_, subject = gene.ref) -> ol_
  se = ol_ %>% as.data.frame() %>% .$subjectHits %>% unique()
  gene.ref[se, ] -> gene.ref
  data.frame(gene_name = gene.ref$gene_name, biotype = gene.ref$gene_biotype) 
}) -> list_1

library(ggvenn)
ggvenn_list = list(PolII = list_1$Polr2b$gene_name, 
                   H3K4me3 = list_1$H3K4me3$gene_name,
                   H3K36me3 = list_1$H3K36me3$gene_name)
ggvenn(ggvenn_list, c("PolII", "H3K4me3", "H3K36me3"))

list_1 %>% map(function(x){
  x$gene_name
}) %>% purrr::reduce(intersect) %>% 
  {gene.ref[gene.ref$gene_name %in% ., ]} %>% 
  .$gene_biotype %>% table() %>% as.data.frame() %>% 
  "colnames<-"(c("biotype", "count")) 

list_1 %>% map(function(x){
  x$gene_name
}) %>% purrr::reduce(intersect) -> intersected_genes

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE35790/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("unique", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  list() %>% map(function(x){
    gsub("/.+/(.+?)_ZT.+", "\\1", x) -> group
    gsub("/.+/.+unique\\.(.+?)\\..+", "\\1", x) -> group_1
    data.frame(files_ = x, group_ = group, group_1 = group_1)
  }) %>% .[[1]] %>% group_by(group_) %>% 
  group_map(function(x,y){
    cbind(x,y) -> df_
    
    environment() -> env_
    names_ = c()
    
    df_ %>% group_by(group_1) %>% 
      group_map(function(x_,y_){
        names_ = c(names_, y_$group_1)
        assign("names_", names_, envir = env_)
        cbind(x_,y_) -> df_
      }) -> list_
    
    list_ %>% map(function(x){
      x$files_ %>% map(function(x){
        ZT = gsub("/.+/.+_(ZT\\d+?)_1.+", "\\1", x)
        read.table(x, sep = "\t", header = T, stringsAsFactors = F) -> df_
        df_[, c(1,7)] %>% 
          "colnames<-"(c("Geneid", ZT))
      }) %>% purrr::reduce(full_join, by = "Geneid")
    }) -> list_1
    names(list_1) = names_
    attr(list_1, "group") <- y$group_
    return(list_1)
  }) -> list_1
list_1 %>% map(function(x){attributes(x)$group}) -> names_
names(list_1) <- names_

list_peaks = list_
list_raw_counts = list_1
rm(list_, list_1)

## Plot intersected genes
list_raw_counts %>% map(function(x){
  x$tss %>% column_to_rownames("Geneid") %>%
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
  df_tss_norm[intersected_genes, ] -> df_tss_norm
  
  x$non %>% column_to_rownames("Geneid") %>%
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm
  df_non_tss_norm[intersected_genes, ] -> df_non_tss_norm
  
  x$tss %>% column_to_rownames("Geneid") %>% 
    {.[intersected_genes, ]} -> df_tss
  
  x$non %>% column_to_rownames("Geneid") %>% 
    {.[intersected_genes, ]} -> df_non_tss
  
  df_tss_ratio = (df_tss_norm + 1)/(df_non_tss_norm + 1)
  
  list_ = list()
  list_[["tss_count"]] = df_tss
  list_[["non_tss_count"]] = df_non_tss 
  list_[["tss_count_norm"]] = df_tss_norm
  list_[["non_tss_count_norm"]] = df_non_tss_norm
  list_[["tss_ratio"]] = df_tss_ratio
  return(list_)
}) -> list_interesected_genes

### PolII (intersected genes)----
list_interesected_genes$Polr2b[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Polr2b (intersected genes)", 
                color = "red", face = "bold", size = 14)) -> p
#4X10

### H3K4me3 (intersected genes)----
list_interesected_genes$H3K4me3[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("H4K4me3 (intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### H3K36me3 (intersected genes)----
list_interesected_genes$H3K36me3[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("H4K36me3 (intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

save.image(file = "./rda/5_PolII_data.rdata", compress = "xz")

rm(list = (ls() %>% .[-matches("^list_", vars = .)]))
list_interesected_genes -> list_interesected_genes_GSE35790
list_peaks -> list_peaks_GSE35790
list_raw_counts -> list_raw_counts_GSE35790
rm(list_interesected_genes, list_peaks, list_raw_counts)
save.image(file = "~/Documents/Tmp_5_PolII_data.rdata", compress = "xz")

load("./rda/5_test_resolution.rda")
rownames(list_interesected_genes_GSE35790$Polr2b$tss_count) -> intersected_genes_GSE35790
list(atac_gene = list_$res_1$JTK.atac$Hepatocytes %>% filter(JTK_BH.Q < 0.05) %>% .$CycID, 
     rna_gene = list_$res_1$JTK.rna$Hepatocytes %>% filter(JTK_BH.Q < 0.05) %>% .$CycID, 
     intersected_gene = intersected_genes_GSE35790) -> ggvenn_list
ggvenn(ggvenn_list, c("atac_gene", "intersected_gene", "rna_gene"))

purrr::reduce(ggvenn_list, intersect) -> intersected_genes_GSE35790_1

### PolII (intersected genes vs rna vs atac)----
list_interesected_genes_GSE35790$Polr2b[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE35790_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Polr2b (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### H3K4me3 (intersected genes vs rna vs atac)----
list_interesected_genes_GSE35790$H3K4me3[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE35790_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("H3K4me3 (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### H3K36me3 (intersected genes vs rna vs atac)----
list_interesected_genes_GSE35790$H3K36me3[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE35790_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("H3K36me3 (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

gc()
load("./rda/3_transferred.rda")
rm(sc)
gc()
gene.ref[gene.ref$gene_name %in% intersected_genes_GSE35790_1, ]$gene_biotype %>% table() %>% 
  as.data.frame() %>% "colnames<-"(c("biotype", "count"))

### run JTKcycle ----
list_interesected_genes_GSE35790 %>% 
  map2(.x = ., .y = names(list_interesected_genes_GSE35790), .f = function(x,y){
    names(x)
    x$tss_count_norm %>% drop_na() %>% rownames_to_column(var = "geneSymbol") %>%
      write.csv("circadian/temp.csv", row.names = F)
    JTK.i=meta2d(infile = 'circadian/temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = c(2,6,10,14,18,22,26), cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni')$meta
    file.remove('circadian/temp.csv')
    x[["JTK_tss_count_norm"]] = JTK.i
    return(x)
  }) -> list_tmp
list_interesected_genes_GSE35790 <- list_tmp
rm(list_tmp)

list_interesected_genes_GSE35790 %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Dbp", "Nfil3", "Nr1d1", "Rora", "Rorc", 
      "Cry1", "Ciart", "Per1", "Per2") -> gene_sele
    x$tss_count_norm[gene_sele, ] %>% rownames_to_column("Gene_name") %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "expression") %>% 
      mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.numeric()) -> df_
    gene_sele %>% map(function(x_){
      df_ %>% filter(Gene_name == x_) %>% 
        ggplot(., aes(x = ZT, y = expression)) + 
        geom_line() + 
        ggtitle(sprintf("%s %s: adj.pval=%s, adj.phase=%s", y, x_, 
                        x$JTK_tss_count_norm %>% as_tibble() %>% 
                          column_to_rownames("CycID") %>% {.[x_, "JTK_BH.Q"]} %>% round(., 3),
                        x$JTK_tss_count_norm %>% as_tibble() %>% 
                          column_to_rownames("CycID") %>% {.[x_, "JTK_adjphase"]} %>% round(., 3)))
    })
  }) -> p_list

names(p_list)
p_list[c("Polr2b", "H3K4me3", "H3K36me3")] -> p_list
purrr::reduce(p_list, c) -> p_list
ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 13, common.legend = T) -> p

#GSE96773 ----
##PolII ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Pol_II.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Pol_II.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

df_ %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  ggtitle("GSE96773 counts/non TSS counts") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_ratio
p_tss_ratio_PolII_GSE96773 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("non TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_non_tss_counts

list_GSE96773 = list()
list_GSE96773[["PolII"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)
#fig 6X12

##Nelf-A ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Nelf-A.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Nelf-A.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

df_ %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  ggtitle("GSE96773 TSS counts/non TSS counts") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_ratio
p_tss_ratio_NelfA_GSE96773 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("non TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_non_tss_counts

list_GSE96773[["Nelf-A"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)
#fig 6X12

##Tbp ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Tbp.+unique.+tss.+", vars = .)] %>% 
  .[-matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_tss
dim(df_tss)

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("Tbp.+unique.+tss.+", vars = .)] %>% 
  .[matches("non", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  map(function(x){
    samples_ = gsub(".+/(.+)_1_sorted.+", "\\1", x)
    read.table(x, sep = "\t", header = TRUE, stringsAsFactors = F) -> df_
    gene = df_$Geneid
    counts = df_[[7]]
    data.frame(gene = df_$Geneid) -> df_
    df_[[samples_]] <- counts
    df_
  }) %>% purrr::reduce(left_join, by = "gene") %>% 
  column_to_rownames("gene") -> df_non_tss
dim(df_non_tss)

df_ = (df_tss+1)/(df_non_tss+1)
df_

df_ %>% 
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
#  mutate(counts = log(counts, base = 2)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_count+1/non_tss_count+1") + 
  ggtitle("GSE96773 TSS counts/non TSS counts") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_ratio
p_tss_ratio_Tbp_GSE96773 = p_tss_ratio

df_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
df_non_tss %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm

df_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_tss_counts

df_non_tss_norm %>% rownames_to_column("gene") %>% 
  pivot_longer(!gene, names_to = "sample", values_to = "counts") %>%
  mutate(counts = log(counts, base = 10)) %>% 
  mutate(sample = factor(sample, levels = gtools::mixedsort(unique(sample)))) %>% 
  ggplot(aes(x = sample, y = counts, fill = sample)) + 
  geom_boxplot() + 
  ylab("log10((Counts+1)/lib_size*C)") + 
  ggtitle("non TSS counts normalized by lib_size") + 
  theme(axis.text.x = element_text(angle = 45)) -> p_non_tss_counts

list_GSE96773[["Tbp"]] = list(
  tss_ratio = df_,
  tss_raw_counts = df_tss,
  non_tss_raw_counts = df_non_tss,
  tss_norm_counts = df_tss_norm,
  non_tss_norm_counts = df_non_tss_norm
)

ggpubr::ggarrange(p_tss_ratio, p_tss_counts, p_non_tss_counts, ncol = 1, nrow = 3, common.legend = T)
#fig 6X12
ggpubr::ggarrange(p_tss_ratio_PolII_GSE96773, 
                  p_tss_ratio_PolII_GSE35790, 
                  p_tss_ratio_NelfA_GSE96773, 
                  p_tss_ratio_H3K4me3_GSE35790, 
                  p_tss_ratio_H3K36me3_GSE35790, 
                  p_tss_ratio_Tbp_GSE96773, 
                  ncol = 3, nrow = 2)

##macs2 peakcalling ----
list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/macs2", 
           pattern = "\\.broadPeak$", full.names = TRUE) %>% 
  list() %>% 
  map(function(x){
    group_ = gsub("/.+/(.+?)_ZT.+", "\\1", x) 
    data.frame(files_ = x, group_ = group_)
  }) %>% .[[1]] %>% group_by(group_) %>% 
  group_map(function(x,y){
    cbind(x,y) -> df_
    df_$files_ %>% 
      map(function(x_){
        env_ = environment()
        tryCatch(read.table(x_, sep = "\t", header = FALSE, stringsAsFactors = F) %>% 
                 filter(V9 > 2) %>% dplyr::select(V1, V2, V3) %>% 
                   {assign("df_1", value = ., envir = env_); print("success")}, 
                 error = function(e){
                   sprintf("%s have zero peak", x_) %>% print() ; 
                   data.frame(V1 = integer(0), V2 = integer(0), V3 = integer(0)) %>% 
                     {assign("df_1", value = ., envir = env_); print("fail")}
                 })
        print(sprintf("Dim: %s %s", nrow(df_1), ncol(df_1)))
        return(df_1)
#        read.table(x_, sep = "\t", header = FALSE, stringsAsFactors = F) -> df_
#        df_ %>% filter(V9 > 2) %>% dplyr::select(V1, V2, V3) -> df_
      }) %>% do.call(rbind, .) %>% arrange(V1, V2) -> df_
    bedtoolsr::bt.merge(df_, d = 1000) -> df_
    assign(sprintf("GSE96773_%s_peaks", y$group_), value = df_, envir = .GlobalEnv)
    attr(df_, "group") <- y$group_
    return(df_)
  }) -> list_peaks_GSE96773
map(list_peaks_GSE96773, function(x){attributes(x)$group}) -> names_
names(list_peaks_GSE96773) <- names_

list_peaks_GSE96773$Pol_II %>% filter(!grepl("_", V1)) %>% .$V1 %>% table()

list_peaks_GSE96773 %>% map(function(x){
  gene.ref = gene.ref[seqnames(gene.ref) %in% paste0("chr", c(1:19, "M", "X", "Y")),]
  x %>% filter(!grepl("_", V1)) -> x
  GRanges(seqnames = x$V1, ranges = IRanges(start = x$V2, end = x$V3)) -> gr_
  findOverlaps(query = gr_, subject = gene.ref) -> ol_
  se = ol_ %>% as.data.frame() %>% .$subjectHits %>% unique()
  gene.ref[se, ] -> gene.ref
  data.frame(gene_name = gene.ref$gene_name, biotype = gene.ref$gene_biotype) 
}) -> list_tmp

library(ggvenn)
ggvenn_list = list(PolII = list_tmp$Pol_II$gene_name, 
                   Tbp = list_tmp$Tbp$gene_name,
                   Nelf_A = list_tmp$`Nelf-A`$gene_name)
ggvenn(ggvenn_list, c("PolII", "Tbp", "Nelf_A"))

list_tmp %>% map(function(x){
  x$gene_name
}) %>% purrr::reduce(intersect) %>% 
  {gene.ref[gene.ref$gene_name %in% ., ]} %>% 
  .$gene_biotype %>% table() %>% as.data.frame() %>% 
  "colnames<-"(c("biotype", "count")) 

list_tmp %>% map(function(x){
  x$gene_name
}) %>% purrr::reduce(intersect) -> intersected_genes_GSE96773

list.files(path = "/Volumes/Crucial_X8/temp/GEO/GSE96773/processed_data",
           pattern = "\\.txt$", 
           full.names = TRUE) %>% 
  .[matches("unique", vars = .)] %>% 
  .[-matches("INPUT", vars = .)] %>% 
  gtools::mixedsort() %>% 
  list() %>% map(function(x){
    gsub("/.+/(.+?)_ZT.+", "\\1", x) -> group
    gsub("/.+/.+unique\\.(.+?)\\..+", "\\1", x) -> group_1
    data.frame(files_ = x, group_ = group, group_1 = group_1)
  }) %>% .[[1]] %>% group_by(group_) %>% 
  group_map(function(x,y){
    cbind(x,y) -> df_
    
    environment() -> env_
    names_ = c()
    
    df_ %>% group_by(group_1) %>% 
      group_map(function(x_,y_){
        names_ = c(names_, y_$group_1)
        assign("names_", names_, envir = env_)
        cbind(x_,y_) -> df_
      }) -> list_
    
    list_ %>% map(function(x){
      x$files_ %>% map(function(x){
        ZT = gsub("/.+/.+_(ZT\\d+?)_1.+", "\\1", x)
        read.table(x, sep = "\t", header = T, stringsAsFactors = F) -> df_
        df_[, c(1,7)] %>% 
          "colnames<-"(c("Geneid", ZT))
      }) %>% purrr::reduce(full_join, by = "Geneid")
    }) -> list_1
    names(list_1) = names_
    attr(list_1, "group") <- y$group_
    return(list_1)
  }) -> list_raw_counts_GSE96773
list_raw_counts_GSE96773 %>% 
  map(function(x){attributes(x)$group}) -> names_
names(list_raw_counts_GSE96773) <- names_

#list_peaks = list_
#list_raw_counts = list_1
#rm(list_, list_1)

## Plot intersected genes
list_raw_counts_GSE96773 %>% map(function(x){
  x$tss %>% column_to_rownames("Geneid") %>%
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_tss_norm
  df_tss_norm[intersected_genes_GSE96773, ] -> df_tss_norm
  
  x$non %>% column_to_rownames("Geneid") %>%
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_non_tss_norm
  df_non_tss_norm[intersected_genes_GSE96773, ] -> df_non_tss_norm
  
  x$tss %>% column_to_rownames("Geneid") %>% 
    {.[intersected_genes_GSE96773, ]} -> df_tss
  
  x$non %>% column_to_rownames("Geneid") %>% 
    {.[intersected_genes_GSE96773, ]} -> df_non_tss
  
  df_tss_ratio = (df_tss_norm + 1)/(df_non_tss_norm + 1)
  
  list_ = list()
  list_[["tss_count"]] = df_tss
  list_[["non_tss_count"]] = df_non_tss 
  list_[["tss_count_norm"]] = df_tss_norm
  list_[["non_tss_count_norm"]] = df_non_tss_norm
  list_[["tss_ratio"]] = df_tss_ratio
  return(list_)
}) -> list_interesected_genes_GSE96773

### PolII (intersected genes)----
list_interesected_genes_GSE96773$Pol_II[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Pol_II (intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### Nelf-A (intersected genes)----
list_interesected_genes_GSE96773$`Nelf-A`[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Nelf-A (intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### Tbp (intersected genes)----
list_interesected_genes_GSE96773$Tbp[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% rownames_to_column("Gene_name") %>% drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) %>% 
        mutate(ZT = factor(ZT, levels = gtools::mixedsort(unique(ZT)))) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Tbp (intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

list(atac_gene = list_$res_1$JTK.atac$Hepatocytes %>% filter(JTK_BH.Q < 0.05) %>% .$CycID, 
     rna_gene = list_$res_1$JTK.rna$Hepatocytes %>% filter(JTK_BH.Q < 0.05) %>% .$CycID, 
     intersected_gene = intersected_genes_GSE96773) -> ggvenn_list
ggvenn(ggvenn_list, c("atac_gene", "intersected_gene", "rna_gene"))

purrr::reduce(ggvenn_list, intersect) -> intersected_genes_GSE96773_1

gene.ref[gene.ref$gene_name %in% intersected_genes_GSE96773_1, ]$gene_biotype %>% table() %>% 
  as.data.frame() %>% "colnames<-"(c("biotype", "count"))

### PolII (intersected genes vs rna vs atac)----
list_interesected_genes_GSE96773$Pol_II[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE96773_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Pol_II (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### Nelf-A (intersected genes vs rna vs atac)----
list_interesected_genes_GSE96773$`Nelf-A`[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE96773_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Nelf-A (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### Tbp (intersected genes vs rna vs atac)----
list_interesected_genes_GSE96773$Tbp[3:5] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% {.[intersected_genes_GSE96773_1, ]} %>% rownames_to_column("Gene_name") %>% 
      drop_na() %>% 
      pivot_longer(!Gene_name, names_to = "ZT", values_to = "counts") -> df_
    
    if (grepl("ratio", y)){
      df_ %>% mutate(counts = log(counts, base = 2)) -> df_
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log2(tss_count+1/non_tss_count+1)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle("TSS counts/non TSS counts") -> p
    }else{
      df_ %>% mutate(counts = log(counts, base = 10)) -> df_ 
      ggplot(df_, aes(x = ZT, y = counts, fill = ZT)) + 
        geom_boxplot() + 
        ylab("log10((Counts+1)/lib_size*C)") + 
        theme(axis.text.x = element_text(angle = 45)) + 
        ggtitle(y) -> p
    }
    
    return(p)
  }) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = T) -> p
ggpubr::annotate_figure(p, top = ggpubr::text_grob("Tbp (new intersected genes)", 
                                                   color = "red", face = "bold", size = 14)) -> p
#4X10

### run JTKcycle ----
list_interesected_genes_GSE96773 %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    names(x)
    x$tss_count_norm %>% drop_na() %>% rownames_to_column(var = "geneSymbol") %>%
      write.csv("circadian/temp.csv", row.names = F)
    JTK.i=meta2d(infile = 'circadian/temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = c(2,6,10,14,18,22), cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni')$meta
    file.remove('circadian/temp.csv')
    x[["JTK_tss_count_norm"]] = JTK.i
    return(x)
  }) -> list_tmp
list_interesected_genes_GSE96773 <- list_tmp
rm(list_tmp)

list_interesected_genes_GSE96773 %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    c("Arntl", "Bhlhe40", "Bhlhe41", "Clock", "Dbp", "Nfil3", "Nr1d1", "Rora", "Rorc", 
      "Cry1", "Ciart", "Per1", "Per2") -> gene_sele
    list_raw_counts_GSE96773[[y]][["tss"]] %>% column_to_rownames("Geneid") -> df_raw
    df_raw %>% mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) -> df_norm
    
    df_norm %>% rownames_to_column("Geneid") %>% 
      write.csv(., file = "circadian/temp.csv", row.names = F)
    JTK.i=meta2d(infile = 'circadian/temp.csv',filestyle = "csv",
                 outputFile = FALSE, timepoints = c(2,6,10,14,18,22), cycMethod = c("JTK"),
                 maxper = 24, minper = 24, combinePvalue = 'bonferroni')$meta
    file.remove('circadian/temp.csv')
    
    df_norm %>% rownames_to_column("Geneid") %>%
      pivot_longer(!Geneid, names_to = "ZT", values_to = "expression") %>% 
      mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.numeric()) -> df_
    
    gene_sele %>% map(function(x_){
      df_ %>% filter(Geneid == x_) %>% 
        ggplot(., aes(x = ZT, y = expression)) + 
        geom_line() + 
        ggtitle(sprintf("%s %s: adj.pval=%s, adj.phase=%s", y, x_, 
                        JTK.i %>% as_tibble() %>%
                          column_to_rownames("CycID") %>% {.[x_, "JTK_BH.Q"]} %>% round(., 3),
                        JTK.i %>% as_tibble() %>%
                          column_to_rownames("CycID") %>% {.[x_, "JTK_adjphase"]} %>% round(., 3)))
    })
  }) -> p_list

names(p_list)
p_list[c("Pol_II", "Nelf-A", "Tbp")] -> p_list
purrr::reduce(p_list, c) -> p_list
ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 13, common.legend = T) -> p

save.image(file = "./rda/5_PolII_data.rdata", compress = "xz")

#pheatmap----
setwd("~/Dropbox/singulomics/")
library(pheatmap)
library(tidyverse)
library(patchwork)
##GSE96773----

### TSS norm counts ----
list_raw_counts_GSE96773$Pol_II$tss
list_raw_counts_GSE96773 %>% map(function(x){
  x$tss %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) %>% 
    apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_
  return(df_)
}) -> list_norm_counts_rescale_GSE96773

list_norm_counts_rescale_GSE96773 %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list

pheatmap_list %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

pheatmap_list %>% map2(.x = ., .y = names(.), .f = function(x,y){
  cutree(x$cutree$tree_2, k = 2) -> tree_
  table(tree_) %>% sort(decreasing = F)
  (tree_ == 2) %>% tree_[.] %>% names() -> genes_
}) -> ggvenn_list

ggvenn::ggvenn(ggvenn_list, c("Pol_II", "Nelf-A", "Tbp"))
ggvenn::ggvenn(ggvenn_list, c("Pol_II", "Tbp"))
#test plotting
#p_list$Tbp$tree_2$gtable %>% ggplotify::as.ggplot()
#p_list$Tbp$tree_2$gtable %>% ggpubr::as_ggplot()
#p_list$Tbp$tree_2$gtable %>% ggplotify::as.grob()
#p_list$Tbp$tree_2$gtable %>% cowplot::as_grob()
#p_list$Pol_II$tree_2

### non TSS norm counts ----
list_raw_counts_GSE96773 %>% map(function(x){
  x$non %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) %>% 
    apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_
  return(df_)
}) -> list_norm_counts_rescale_GSE96773_non_tss

list_norm_counts_rescale_GSE96773_non_tss %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_non_tss

pheatmap_list_non_tss %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

pheatmap_list %>% map2(.x = ., .y = names(.), .f = function(x,y){
  cutree(x$cutree$tree_2, k = 2) -> tree_
  table(tree_) %>% sort(decreasing = F)
  (tree_ == 2) %>% tree_[.] %>% names() -> genes_
}) -> ggvenn_list

ggvenn::ggvenn(ggvenn_list, c("Pol_II", "Nelf-A", "Tbp"))
ggvenn::ggvenn(ggvenn_list, c("Pol_II", "Tbp"))

### TSS/non_tss counts ----
list_raw_counts_GSE96773 %>% map(function(x){
  x$tss %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) ->> df_tss
  x$non %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) -> df_non_tss
  df_tss/df_non_tss -> df_ratio
#  df_ratio %>% head()
    apply(df_ratio, 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_ratio
    return(df_ratio)
}) -> list_norm_counts_rescale_GSE96773_tss_ratio

list_norm_counts_rescale_GSE96773_tss_ratio %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_tss_ratio
pheatmap_list_tss_ratio -> pheatmap_list_tss_ratio_GSE96773

pheatmap_list_tss_ratio %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

## GSE35790----
### TSS norm counts ----
list_raw_counts_GSE35790$Polr2$tss
list_raw_counts_GSE35790 %>% map(function(x){
  x$tss %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) %>% 
    apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_
  return(df_)
}) -> list_norm_counts_rescale_GSE35790

list_norm_counts_rescale_GSE35790 %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list

pheatmap_list %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

### non TSS norm counts ----
list_raw_counts_GSE35790 %>% map(function(x){
  x$non %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){((x+1)/sum(x))*10^6})) %>% 
    apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_
  return(df_)
}) -> list_norm_counts_rescale_GSE35790_non_tss

list_norm_counts_rescale_GSE35790_non_tss %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_non_tss

pheatmap_list_non_tss %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

list_raw_counts_GSE96773 %>% map(function(x){
  x$tss %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) -> df_tss
  x$non %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) -> df_non_tss
  df_tss/df_non_tss -> df_ratio
#  df_ratio %>% head()
    apply(df_ratio, 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_ratio
    return(df_ratio)
}) -> list_norm_counts_rescale_GSE96773_tss_ratio

list_norm_counts_rescale_GSE96773_tss_ratio %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_tss_ratio
pheatmap_list_tss_ratio -> pheatmap_list_tss_ratio_GSE96773
rm(pheatmap_list_tss_ratio)

pheatmap_list_tss_ratio %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

### TSS/non_tss counts ----
list_raw_counts_GSE35790 %>% map(function(x){
  x$tss %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) -> df_tss
  x$non %>% 
    column_to_rownames("Geneid") %>% 
    mutate(across(everything(), function(x){x+1})) -> df_non_tss
  df_tss/df_non_tss -> df_ratio
  #  df_ratio %>% head()
  apply(df_ratio, 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
    t() %>% base::as.data.frame() -> df_ratio
  return(df_ratio)
}) -> list_norm_counts_rescale_GSE35790_tss_ratio

list_norm_counts_rescale_GSE35790_tss_ratio %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_tss_ratio
pheatmap_list_tss_ratio -> pheatmap_list_tss_ratio_GSE35790

pheatmap_list_tss_ratio %>% map2(.x = ., .y = names(.), .f = function(x,y){
  x$plot$tree_2 + 
    ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  return(p)
}) -> p_list

ggpubr::ggarrange(plotlist = p_list, nrow = 3, ncol = 1, common.legend = F) -> p
rm(p)

#sc ATAC ----
readRDS("~/Desktop/output.tss.raw.list.rds") -> output.tss.raw.list
readRDS("~/Desktop/output.non.tss.raw.list.rds") -> output.non.tss.raw.list

output.tss.raw.list$Sox17$celltype %>% unique() %>% 
  map(function(x){
    i = 0
    output.tss.raw.list %>% 
      map2(.x = ., .y = names(.), .f = function(x_, y_){
        i + 1 ->> i
        sprintf("cell type: %s, gene no: %s, TSS", x, i) %>% print()
      x_ %>% filter(celltype == x) %>% 
          group_by(ZT) %>% summarise(mean = mean(mean)) -> df_
       list_ = list() 
       list_[[y_]] = df_$mean
       list_ %>% as.data.frame() %>% "rownames<-"(., sprintf("ZT%s", seq(2,22,4))) %>% 
         t() %>% as.data.frame()
    }) %>% do.call(rbind, .) -> df_tss
    
    i = 0
    output.non.tss.raw.list %>% 
      map2(.x = ., .y = names(.), .f = function(x_, y_){
        i + 1 ->> i
        sprintf("cell type: %s, gene no: %s, non TSS", x, i) %>% print()
        x_ %>% filter(celltype == x) %>% 
          group_by(ZT) %>% summarise(mean = mean(mean)) -> df_
        list_ = list() 
        list_[[y_]] = df_$mean
        list_ %>% as.data.frame() %>% "rownames<-"(., sprintf("ZT%s", seq(2,22,4))) %>% 
          t() %>% as.data.frame()
      }) %>% do.call(rbind, .) -> df_non_tss
    
    intersect(rownames(df_tss), rownames(df_non_tss)) -> common_genes
    df_tss[common_genes, ]/df_non_tss[common_genes, ] -> df_ratio
    df_ratio %>% filter_all(all_vars(!is.infinite(.))) %>% 
      drop_na() -> df_final
    attr(df_final, "celltype") <- x
    return(df_final)
  }) -> list_atac_tss_ratio
map(list_atac_tss_ratio, function(x){attributes(x)$celltype}) -> names(list_atac_tss_ratio)

list_atac_tss_ratio %>% 
  map(function(x){
    apply(x, 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
      t() %>% base::as.data.frame() -> df_ratio
    return(df_ratio)
  }) -> list_atac_tss_ratio_rescale

list_atac_tss_ratio_rescale %>% 
  map(function(x){
    tree_2 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 2) 
    tree_3 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 3) 
    tree_4 = pheatmap(x, cluster_rows = T, cluster_cols = F, 
                      show_rownames = F, cutree_rows = 4) 
    
    p_list = list(tree_2 = ggplotify::as.ggplot(tree_2$gtable), 
                  tree_3 = ggplotify::as.ggplot(tree_3$gtable),
                  tree_4 = ggplotify::as.ggplot(tree_4$gtable))
    
    cutree_list = list(tree_2 = tree_2$tree_row, 
                       tree_3 = tree_3$tree_row,
                       tree_4 = tree_4$tree_row)
    
    list(cutree = cutree_list, plot = p_list) -> list_
    return(list_)
  }) -> pheatmap_list_atac_tss_ratio

pheatmap_list_atac_tss_ratio %>% 
  map2(.x = ., .y = names(.), function(x,y){
    x$plot$tree_2 + ggtitle(y) + 
    theme(plot.title = element_text(hjust = 0.5)) -> p
  }) -> p_list
ggpubr::ggarrange(p_list[[1]], 
                  p_list[[2]],
                  p_list[[3]], 
                  p_list[[4]], 
                  nrow = 2, ncol = 2, common.legend = F) -> p

cutree(pheatmap_list_atac_tss_ratio$Hepatocytes$cutree$tree_2, k = 2) %>% 
  sort(decreasing = T) %>% names() %>% 
  {list_atac_tss_ratio_rescale$Hepatocytes[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)}

#load("~/Desktop/tmp.rda")
intersect(pheatmap_list_tss_ratio_GSE96773$Pol_II$cutree$tree_2$labels, 
          pheatmap_list_atac_tss_ratio$Hepatocytes$cutree$tree_2$labels) -> common_genes 

cutree(pheatmap_list_atac_tss_ratio$Hepatocytes$cutree$tree_2, k = 2) %>% .[common_genes] %>%
  sort(decreasing = T) %>% names() -> se_genes
se_genes %>%
  {list_atac_tss_ratio_rescale$Hepatocytes[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>% 
  {ggplotify::as.ggplot(.) + ggtitle("ATAC hepatocytpes TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_1

se_genes %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$Pol_II[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 PolII Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_2

se_genes %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$Polr2b[., ]} %>% 
  dplyr::select(-ZT26) %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 PolII Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_3

se_genes %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$`Nelf-A`[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Nelf-A Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_4

se_genes %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$Tbp[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Tbp Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_5

se_genes %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$H3K4me3[., ]} %>% 
  dplyr::select(-ZT26) %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 H3K4me3 Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_6

readRDS("~/Desktop/otuput.rna.list.rds") -> output.rna.list
output.rna.list %>% 
#  .[1:10] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% filter(celltype == "Hepatocytes") %>% 
      group_by(ZT) %>% 
      summarise(mean = mean(mean)) -> df_
    list(mean = df_$mean) %>% 
      as.data.frame() %>% t() %>% 
      as.data.frame() %>% 
      "rownames<-"(., y) %>% 
      "colnames<-"(., sprintf("ZT%s", seq(2,22,4)))
  }) %>% do.call(rbind, .) %>%
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  {.[se_genes, ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("scRNA gene expression") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_7

se_genes %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$H3K36me3[., ]} %>% 
  dplyr::select(-ZT26) %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 H3K36me3 Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_8

#ggpubr::ggarrange(p_1, p_2, p_3, p_4, p_5, p_6, p_7, nrow = 3, ncol = 3, common.legend = F) -> p
ggpubr::ggarrange(p_4, p_1, p_2, p_7, p_5, p_8, p_6, nrow = 4, ncol = 2) -> p

intersect(
list_interesected_genes_GSE96773$Pol_II$tss_ratio %>% rownames(), 
list_interesected_genes_GSE35790$Polr2b$tss_ratio %>% rownames()
) -> common_genes_1

(common_genes_1 %in% rownames(list_atac_tss_ratio_rescale$Hepatocytes)) %>% 
  common_genes_1[.] -> common_genes_1

common_genes_1 %>%
  {list_atac_tss_ratio_rescale$Hepatocytes[., ]} %>% 
  {pheatmap(., cluster_rows = T, cluster_cols = F, 
            show_rownames = F, cutree_rows = 1)} -> p_1

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_atac_tss_ratio_rescale$Hepatocytes[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>% 
{ggplotify::as.ggplot(.) + ggtitle("ATAC hepatocytpes TSS ratio") + 
    theme(plot.title = element_text(hjust = 0.5))} -> p_2

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$Pol_II[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 PolII Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_3

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$Polr2b[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 PolII Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_4

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$`Nelf-A`[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Nelf-A Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_5

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_norm_counts_rescale_GSE96773_tss_ratio$Tbp[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Tbp Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_6

p_1$tree_row$labels[p_1$tree_row$order] %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$H3K4me3[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 H3K4me3 Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p_7

#Use highly variable genes
readRDS("~/Downloads/rds/highly_variable_genes.rds") -> highly_variable_genes
list_atac_tss_ratio_rescale$Hepatocytes[highly_variable_genes, ] %>% drop_na() %>% 
  pheatmap(., cluster_rows = T, cluster_cols = F, 
           show_rownames = F, cutree_rows = 4) -> p1
c(1,3,2,4) %>% 
  map(function(x){
    cutree(p1$tree_row, k = 4) %>% .[. == x] %>% 
      names()
  }) %>% unlist() -> highly_variable_genes
list_atac_tss_ratio_rescale$Hepatocytes[highly_variable_genes, ] %>% drop_na() %>% 
  pheatmap(., cluster_rows = F, cluster_cols = F, 
           show_rownames = F) %>% 
  {ggplotify::as.ggplot(.) + ggtitle("ATAC hepatocytes TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p2

highly_variable_genes %>% 
{list_norm_counts_rescale_GSE96773_tss_ratio$`Nelf-A`[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Nelf-A Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p3

highly_variable_genes %>%
{list_norm_counts_rescale_GSE96773_tss_ratio$Pol_II[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 PolII Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p4

highly_variable_genes %>%
{list_norm_counts_rescale_GSE96773_tss_ratio$Tbp[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE96773 Tbp Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p5

highly_variable_genes %>% 
{list_norm_counts_rescale_GSE35790_tss_ratio$H3K4me3[., ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 H3K4me3 Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p6

output.rna.list %>% 
  #  .[1:10] %>% 
  map2(.x = ., .y = names(.), .f = function(x,y){
    x %>% filter(celltype == "Hepatocytes") %>% 
      group_by(ZT) %>% 
      summarise(mean = mean(mean)) -> df_
    list(mean = df_$mean) %>% 
      as.data.frame() %>% t() %>% 
      as.data.frame() %>% 
      "rownames<-"(., y) %>% 
      "colnames<-"(., sprintf("ZT%s", seq(2,22,4)))
  }) %>% do.call(rbind, .) %>%
  apply(., 1, function(x){scales::rescale(x, to = c(0,1))}) %>% 
  t() %>% base::as.data.frame() %>% 
  {.[highly_variable_genes, ]} %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("scRNA gene expression") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p7

highly_variable_genes %>% 
  {list_norm_counts_rescale_GSE35790_tss_ratio$H3K36me3[., ]} %>% 
  dplyr::select(-ZT26) %>% 
  {pheatmap(., cluster_rows = F, cluster_cols = F, 
            show_rownames = F)} %>%
  {ggplotify::as.ggplot(.) + ggtitle("GSE35790 H3K36me3 Chip TSS ratio") + 
      theme(plot.title = element_text(hjust = 0.5))} -> p8

ggpubr::ggarrange(p3, p2, p4, p7, p5, p8, p6, ncol = 2, nrow = 4)
length(highly_variable_genes)

save.image(file = "./rda/5_PolII_data.rdata", compress = "xz")

list_atac_tss_ratio_rescale$Hepatocytes %>% 
  pivot_longer(everything(), names_to = "sample", values_to = "value") %>% 
  mutate(sample = factor(sample, levels = c("ZT2", "ZT6", "ZT10", "ZT14", "ZT18", "ZT22"))) %>%
  ggplot(., aes(x = sample, y = value, fill = sample)) + 
  geom_boxplot() + 
  ylab("tss_counts/non_tss_counts") -> p_atac

list(
  p_atac + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)),
  
  p_tss_ratio_PolII_GSE96773 + 
    ggtitle(NULL) +
    scale_fill_discrete(name = "sample", labels = sprintf("ZT%s", seq(2,22,4))) + 
    scale_x_discrete(labels=sprintf("ZT%s", seq(2,22,4))) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)),
  
  p_tss_ratio_Tbp_GSE96773 + 
    ggtitle(NULL) +
    scale_fill_discrete(name = "sample", labels = sprintf("ZT%s", seq(2,22,4))) + 
    scale_x_discrete(labels=sprintf("ZT%s", seq(2,22,4))) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)),
  
  p_tss_ratio_NelfA_GSE96773 + 
    ggtitle(NULL) +
    scale_fill_discrete(name = "sample", labels = sprintf("ZT%s", seq(2,22,4))) + 
    scale_x_discrete(labels=sprintf("ZT%s", seq(2,22,4))) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)),
  
#  p_tss_ratio_H3K4me3_GSE35790 + 
#    ggtitle(NULL) +
#    scale_fill_discrete(name = "sample", labels = sprintf("ZT%s", seq(2,22,4))) + 
#    scale_x_discrete(labels=sprintf("ZT%s", seq(2,22,4))) + 
#    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)),
  
#  p_tss_ratio_H3K36me3_GSE35790 + 
#    ggtitle(NULL) +
#    scale_fill_discrete(name = "sample", labels = sprintf("ZT%s", seq(2,22,4))) + 
#    scale_x_discrete(labels=sprintf("ZT%s", seq(2,22,4))) + 
#    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))

  p_tss_ratio_H3K4me3_GSE35790$data %>% 
    mutate(sample = sprintf("ZT%s", gsub(".+_ZT(\\d+)", "\\1", sample) %>% as.numeric())) %>% 
    filter(sample %in% sprintf("ZT%s", seq(2,22,4))) %>% 
    mutate(sample = factor(sample, levels = sprintf("ZT%s", seq(2,22,4)))) %>%
    ggplot(aes(x = sample, y = counts, fill = sample)) + 
    geom_boxplot() + 
    ylab("tss_count+1/non_tss_count+1"), 

  p_tss_ratio_H3K36me3_GSE35790$data %>% 
    mutate(sample = sprintf("ZT%s", gsub(".+_ZT(\\d+)", "\\1", sample) %>% as.numeric())) %>% 
    filter(sample %in% sprintf("ZT%s", seq(2,22,4))) %>% 
    mutate(sample = factor(sample, levels = sprintf("ZT%s", seq(2,22,4)))) %>%
    ggplot(aes(x = sample, y = counts, fill = sample)) + 
    geom_boxplot() + 
    ylab("tss_count+1/non_tss_count+1") 
) -> p_list

ggpubr::ggarrange(plotlist = p_list, ncol = 3, nrow = 2, common.legend = T, legend = "right")
