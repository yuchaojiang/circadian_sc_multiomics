library(shiny)
library(bslib)
library(DT)
library(tidyverse)

getwd()

#Prepare_data_table ----
readxl::read_xlsx("./Data/00_Supp_Table_1.xlsx", sheet = "Supplementary_Table_1", skip = 1) -> df_

df_[,1:16] -> df_Hep_RNA
colnames(df_Hep_RNA) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Hep_RNA)
head(df_Hep_RNA)

df_[,18:33] -> df_Hep_ATAC
colnames(df_Hep_ATAC) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Hep_ATAC)
head(df_Hep_ATAC)

df_[,35:50] -> df_Endo_RNA
colnames(df_Endo_RNA) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Endo_RNA)
head(df_Endo_RNA)

df_[,52:67] -> df_Endo_ATAC
colnames(df_Endo_ATAC) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Endo_ATAC)
head(df_Endo_ATAC)

df_[,69:84] -> df_Fibro_RNA
colnames(df_Fibro_RNA) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Fibro_RNA)
head(df_Fibro_RNA)

df_[,86:101] -> df_Fibro_ATAC
colnames(df_Fibro_ATAC) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Fibro_ATAC)
head(df_Fibro_ATAC)

df_[,103:118] -> df_Kupffer_RNA
colnames(df_Kupffer_RNA) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Kupffer_RNA)
head(df_Kupffer_RNA)

df_[,120:135] -> df_Kupffer_ATAC
colnames(df_Kupffer_ATAC) %>% gsub("(.+?)\\.\\..+", "\\1", .) -> colnames(df_Kupffer_ATAC)
head(df_Kupffer_ATAC)

list(Hep_RNA=df_Hep_RNA, Hep_ATAC=df_Hep_ATAC,
     Endo_RNA=df_Endo_RNA, Endo_ATAC=df_Endo_ATAC,
     Fibro_RNA=df_Fibro_RNA, Fibro_ATAC=df_Fibro_ATAC,
     Kupffer_RNA=df_Kupffer_RNA, Kupffer_ATAC=df_Kupffer_ATAC) %>% 
  purrr::map2(.x=.,.y=names(.),.f=function(df_, group_){
    group_ %>% gsub("(.+)_(.+)", "\\1", .) -> celltype_
    c(Hep = "Hepatocytes", Endo = "Endothelial cells", Fibro = "Fibroblasts", Kupffer = "Kupffer cells")[celltype_] %>% unname() -> celltype_
    group_ %>% gsub("(.+)_(.+)", "\\2", .) -> assay_
    df_ %>% dplyr::mutate(Assay = assay_, .before = 1) %>% 
      dplyr::mutate(CellType = celltype_, .after = 1) -> df_
  }) %>% do.call(rbind, .) -> df_all
celltype_specific_analysis = df_all
save(celltype_specific_analysis, file = "./Data/celltype_specific_analysis.RData")

readxl::read_xlsx("./Data/00_Supp_Table_3.xlsx", sheet = "Sheet1") -> df_
CRE_linkage_analysis = as.data.frame(df_)
save(CRE_linkage_analysis, file = "./Data/CRE_linkage_analysis.RData")

readxl::read_xlsx("./Data/00_Supp_Table_4.xlsx", sheet = "Sheet1") -> df_
TRIPOD_analysis = as.data.frame(df_)
save(TRIPOD_analysis, file = "./Data/TRIPOD_analysis.RData")

#Make preloaded plots ----
load("./Data/celltype_specific_analysis.RData")
colnames(celltype_specific_analysis)[c(3, 1:2, 4:18)] %>% celltype_specific_analysis[,.] -> celltype_specific_analysis
celltype_specific_analysis = as.data.frame(celltype_specific_analysis)
celltype_specific_analysis %>% dplyr::mutate(Plot = 1:nrow(.) %>% sprintf("P_%s", .), .before=1) -> celltype_specific_analysis

c("RNA", "gene_activity") %>% 
  "names<-"(.,.) %>% 
#  .[1] %>% 
  map(function(assay_){
    print(assay_)
    c("Hepatocytes", "Endothelial_cells", "Fibroblasts", "Kupffer_cells") %>% 
      "names<-"(.,.) %>% 
#      .[1] %>% 
      map(function(celltype_){
        print(celltype_)
        list.files("~/Dropbox/singulomics/github_rda/output/Celltype_specific/Raw_data_downsampled/replicate_10", full.names = T) %>% 
          {.[grepl(assay_, .)]} %>% {.[grepl(celltype_, .)]} %>% gtools::mixedsort() -> files_       
#        seeds_ = c()
        files_ %>% 
#          .[1:3] %>% 
          map(function(file_){
            print(file_)
            seed_ = gsub("/.+/.+_(seed_\\d+?)_.+", "\\1", file_)
#            seeds_ <<- c(seeds_, seed_)
            read.csv(file_, header = T, stringsAsFactors = F) -> df_
            df_ %>% pivot_longer(-Gene, names_to = "Group", values_to = "Exprs") -> df_
            df_ %>% dplyr::mutate(ZT = gsub("ZT(\\d+)_.+", "\\1", Group) %>% as.numeric()) %>% 
              dplyr::mutate(REP = gsub(".+_REP(\\d+)", "\\1", Group) %>% as.numeric()) -> df_
            df_ %>% dplyr::select(-Group) -> df_
            df_ %>% dplyr::mutate(Celltype = celltype_, Assay = assay_, Seed = seed_) -> df_
            df_[,c("Celltype", "Assay", "Gene", "ZT", "REP", "Seed", "Exprs")] -> df_
          }) %>% do.call(rbind, .) -> df_
        df_ %>% group_by(Celltype, Assay, Gene, ZT, REP) %>% 
          summarise(Exprs = mean(Exprs), .groups = "drop") %>% 
          ungroup() -> df_
        df_ %>% group_by(Celltype, Assay, Gene, ZT) %>% 
          summarize(Mean_exprs = mean(Exprs), SD_exprs = sd(Exprs), .groups = "drop") -> df_
      }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) -> df_plot
saveRDS(df_plot, file = "./Data/celltype_specific_analysis_plot_df.RData")

1:nrow(celltype_specific_analysis) %>% 
  "names<-"(., sprintf("P_%s", .)) %>% 
#  .[1:100] %>% 
  purrr::map(function(i_){
    print(i_)
    celltype_specific_analysis[i_, ] -> df_
    gene_ = df_$Gene
    assay_ = df_$Assay
    celltype_ = df_$CellType
    BH.Q = df_$cauchy_BH.Q %>% sprintf("%.2e", .)
    
    df_plot %>% dplyr::filter(Assay == assay_, Gene == gene_, Celltype == celltype_) -> df_plot
    
    sprintf("%s,%s,%s", assay_, celltype_, gene_)
    if (assay_ == "RNA"){
      color_ = "blue"
      ylab_ = "RNA expression"
    }else{
      color_ = "red"
      ylab_ = "ATAC activity"
    }
    
    df_plot %>% 
      ggplot(aes(x = ZT, y = Mean_exprs)) +
      geom_line(color = color_) + 
      geom_ribbon(aes(ymax = Mean_exprs+SD_exprs, ymin = Mean_exprs-SD_exprs), alpha = 0.2, color = NA, fill = color_) + 
      ylab(ylab_) +
      ggtitle(sprintf("%s (%s)\nBH.Q=%s", gene_, celltype_, BH.Q)) + 
      theme_minimal()
  }) -> p_list

load("./Data/CRE_linkage_analysis.RData")
CRE_linkage_analysis = as.data.frame(CRE_linkage_analysis)
CRE_linkage_analysis %>% dplyr::mutate(Plot = 1:nrow(.) %>% sprintf("P_%s", .), .before=1) -> CRE_linkage_analysis

readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_rna.rds") -> output.rna.list
readRDS("~/Dropbox/singulomics/github_rda/output/Hepatocytes_transient/normalized_data/res_2.5_gene_activity.rds") -> output.atac.list

1:nrow(CRE_linkage_analysis) %>% 
  "names<-"(., sprintf("P_%s", .)) %>% 
  .[1:10] %>% 
  purrr::map(function(i_){
    CRE_linkage_analysis[i_, ] -> df_
    df_$Gene -> gene_
    output.rna.list[[gene_]]
  })

# Shiny app ----
ui <- page_sidebar(
  title = "Single-Cell Multiomic Analysis of Circadian Rhythmicity in Mouse Liver",
  sidebar = sidebar(
    accordion(
      id = "nav",                   # <- give the accordion an id
      accordion_panel("Cell-type specific analysis", value = "cell"),
      accordion_panel("Hepatocytes CRE linkage analysis (Spatial & Circadian)",   value = "spatiotemp"),
      accordion_panel("Hepatocytes TRIPOD analysis",        value = "tripod"),
      multiple = FALSE               # one open at a time (set TRUE if you prefer)
    )
  ),
  # MAIN CONTENT: single table output that we’ll feed different data
  card(
    card_header(textOutput("title")),
    DTOutput("table")
  )
)

load("./Data/celltype_specific_analysis.RData")
celltype_specific_analysis = as.data.frame(celltype_specific_analysis)
load("./Data/CRE_linkage_analysis.RData")
CRE_linkage_analysis = as.data.frame(CRE_linkage_analysis)
load("./Data/TRIPOD_analysis.RData")
TRIPOD_analysis = as.data.frame(TRIPOD_analysis)

server <- function(input, output, session) {
  # example datasets you want to show
  tbl_cell      <- celltype_specific_analysis
  tbl_spatiotemp<- CRE_linkage_analysis
  tbl_tripod    <- TRIPOD_analysis
  
  current_data <- reactive({
    switch(input$nav,
           "cell"      = tbl_cell,
           "beyond"    = tbl_beyond,
           "spatiotemp"= tbl_spatiotemp,
           "tripod"    = tbl_tripod,
           tbl_cell)  # default
  })
  
  output$title <- renderText({
    c(cell="Cell-type specific analysis",
      spatiotemp="Hepatocytes CRE linkage analysis (Spatial & Circadian)",
      tripod="Hepatocytes TRIPOD analysis")[input$nav %||% "cell"]
  })
  
  output$table <- renderDataTable({
    datatable(
      current_data(),
      filter = list(position = "top", clear = FALSE),  # << add column filters
      options = list(
        pageLength = 10,
        autoWidth  = TRUE,
        dom = 'lrtip'   # hide global search ('f') but keep column filters
      )
    )
  })
}

shinyApp(ui, server)
