library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(htmltools)
library(ggplot2)
library(patchwork)

# Load functions ----
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
    ggtitle(sprintf("TF: %s\nTarget gene: %s\nCRE: %s", TF_, target_gene_, CRE_)) + 
    theme_minimal()
}

#Load your table data (must include a "Plot" id column like "P_1") ----
load("./Data/celltype_specific_analysis.RData")
colnames(celltype_specific_analysis)[c(3, 1:2, 4:18)] %>% celltype_specific_analysis[,.] -> celltype_specific_analysis
celltype_specific_analysis = as.data.frame(celltype_specific_analysis)
celltype_specific_analysis %>% dplyr::mutate(Plot = 1:nrow(.) %>% sprintf("P_%s", .), .before=1) -> celltype_specific_analysis

load("./Data/CRE_linkage_analysis.RData")
CRE_linkage_analysis = as.data.frame(CRE_linkage_analysis)
CRE_linkage_analysis %>% dplyr::mutate(Plot = 1:nrow(.) %>% sprintf("P_%s", .), .before=1) -> CRE_linkage_analysis


load("./Data/TRIPOD_analysis.RData")
colnames(TRIPOD_analysis) %>% stringr::str_replace(., "^\\w", toupper) -> colnames(TRIPOD_analysis)
TRIPOD_analysis = as.data.frame(TRIPOD_analysis)
TRIPOD_analysis %>% dplyr::mutate(Plot = 1:nrow(.) %>% sprintf("P_%s", .), .before=1) -> TRIPOD_analysis

#load plot data ----
readRDS("./Data/celltype_specific_analysis_plot_df.RData") -> celltype_specific_analysis_plot_df
celltype_specific_analysis_plot_df %>% dplyr::mutate(Assay = case_when(
  Assay == "RNA" ~ "RNA",
  TRUE ~ "ATAC"
)) -> celltype_specific_analysis_plot_df
celltype_specific_analysis_plot_df %>% dplyr::mutate(Celltype = gsub("_", " ", Celltype)) -> celltype_specific_analysis_plot_df
celltype_specific_analysis %>% dplyr::mutate(group = sprintf("%s_%s_%s", Assay, CellType, Gene)) %>% dplyr::select(group, Plot) -> df_tmp
celltype_specific_analysis_plot_df %>% dplyr::mutate(group = sprintf("%s_%s_%s", Assay, Celltype, Gene), .before = 1) -> celltype_specific_analysis_plot_df
left_join(x = celltype_specific_analysis_plot_df, y = df_tmp, by = "group") %>% dplyr::select(-group) -> celltype_specific_analysis_plot_df
celltype_specific_analysis_plot_df %>% dplyr::filter(!is.na(Plot)) -> celltype_specific_analysis_plot_df
rm(df_tmp)

load("./Data/CRE_linkage_analysis_plot_dat.RData")
load("./Data/CRE_linkage_analysis_cor_plot_dat.RData")

load("./Data/TRIPOD_analysis_plot_dat.RData")

#Server ----
server <- function(input, output, session) {
  
  # ---- Which table to show in the main area ----
  current_data <- reactive({
    switch(input$nav,
           "cell"       = celltype_specific_analysis,
           "spatiotemp" = CRE_linkage_analysis,
           "tripod"     = TRIPOD_analysis,
           celltype_specific_analysis
    )
  })
  
  titles <- c(
    cell       = "Cell-type specific analysis",
    spatiotemp = "Hepatocytes CRE linkage analysis (Spatial & Circadian)",
    tripod     = "Hepatocytes TRIPOD analysis"
  )
  
  output$title <- renderText({
    nav <- input$nav %||% "cell"     # shiny’s “or” helper
    titles[[nav]]
  })
  
  # ---- Render the table; first column is a clickable link that sends (ds, id) to the server ----
  output$table <- renderDT({
    dat <- current_data()
    ds  <- input$nav %||% "cell"
    
    if ("Plot" %in% names(dat)) {
      id <- htmlEscape(dat$Plot)
      dat$Plot <- sprintf(
        '<a href="#" class="plot-link" data-id="%s" data-ds="%s">%s</a>',
        id, ds, id
      )
    }
    
    datatable(
      dat,
      escape  = FALSE,                                # keep HTML in Plot column
      filter  = list(position = "top", clear = FALSE),
      options = list(pageLength = 10, autoWidth = TRUE, dom = "lrtip"),
      callback = DT::JS(
        # Use DataTables click -> send to Shiny
        "table.on('click', 'a.plot-link', function(e){",
        "  e.preventDefault();",
        "  var id = $(this).data('id');",
        "  var ds = $(this).data('ds');",
        "  Shiny.setInputValue('plot_click', {id:id, ds:ds}, {priority:'event'});",
        "});"
      )
    )
  })
  
  # ================== Plot builders (one per dataset) ==================
  
  build_plot_cell <- function(id) {
    df_ <- celltype_specific_analysis_plot_df[celltype_specific_analysis_plot_df$Plot == id, , drop = FALSE]
    #    print(sprintf("id=%s", id))
    df_$Celltype %>% unique() -> celltype_
    df_$Assay %>% unique() -> assay_
    df_$Gene %>% unique() -> gene_
    BH.Q = celltype_specific_analysis %>% dplyr::filter(Plot == id) %>% .$cauchy_BH.Q %>% sprintf("%.2e", .)
    
    if (assay_ == "RNA"){
      color_ = "blue"
      ylab_ = "RNA expression"
    }else{
      color_ = "red"
      ylab_ = "ATAC activity"
    }
    
    print(df_)
    df_ %>% 
      ggplot(aes(x = ZT, y = Mean_exprs)) +
      geom_line(color = color_) + 
      geom_ribbon(aes(ymax = Mean_exprs+SD_exprs, ymin = Mean_exprs-SD_exprs), alpha = 0.2, color = NA, fill = color_) + 
      ylab(ylab_) +
      ggtitle(sprintf("%s (%s)\nBH.Q=%s", gene_, celltype_, BH.Q)) + 
      theme_minimal() -> p 
    return(p)
  }
  
  build_plot_spatiotemp <- function(id) {
    CRE_linkage_analysis %>% dplyr::filter(Plot == id) -> df_
    gene_ = df_$Gene
    peak_ = df_$Peak
    Spatial_r = df_$Spatial_r %>% round(., 2)
    Temporal_r = df_$Temporal_r %>% round(., 2)
    pval_list_CRE_linkage %>% dplyr::filter(Gene == gene_) -> df_pval
    Circadian_p_adj = df_pval %>% .$Circadian_p.adj %>% sprintf("%.2e", .)
    Transient_p_adj = df_pval %>% .$Transient_p.adj %>% sprintf("%.2e", .)
    Circadian_transient_p_adj = df_pval %>% .$Circadian_transient_p.adj %>% sprintf("%.2e", .)
    
    plot_list_CRE_linkage[[gene_]] -> df_plot
    df_plot %>% dplyr::mutate(ZT = gsub("ZT(\\d+)", "\\1", ZT) %>% as.numeric()) -> df_plot
    df_plot %>% 
      ggplot(aes(x = ZT, y = Mean_expression, group = cluster, color = pseudotime)) + 
      geom_point() +
      geom_line() + 
      theme_classic() + 
      scale_x_continuous(breaks = seq(0, 24, 4)) -> p 
    
    p + 
      ggtitle(sprintf("%s\nC p.adj=%s\nS p.adj=%s\nS&C p.adj=%s", gene_, Circadian_p_adj, Transient_p_adj, Circadian_transient_p_adj)) -> p
    p + 
      scale_x_continuous(breaks = seq(2, 22, 4)) -> p_1
    
    CRE_linkage_cor_plot_list[[gene_]][["temporal"]][, c("Gene_expr", peak_)] %>% as.data.frame() -> df_temporal
    colnames(df_temporal) = c("RNA_expr", "Peak_accessibility")
    df_temporal[["ZT"]] = rownames(df_temporal) %>% gsub("ZT(\\d+)_.+", "\\1", .) %>% as.numeric() %>% as.factor()
    df_temporal %>% 
      ggplot(aes(x = RNA_expr, y = Peak_accessibility, color = ZT, group = ZT)) + 
      geom_point() + 
      ylab("ATAC peak accessibility (Circadian)") + 
      xlab("RNA expression (Circadian)") + 
      ggtitle(sprintf("Circadian r=%s", Temporal_r)) + 
      theme_minimal() -> p_2
    
    CRE_linkage_cor_plot_list[[gene_]][["spatial"]][, c("Gene_expr", peak_)] %>% as.data.frame() -> df_spatial
    colnames(df_spatial) = c("RNA_expr", "Peak_accessibility")
    df_spatial[["Zone"]] = rownames(df_spatial) %>% gsub("zone(\\d+)_.+", "\\1", .) %>% as.numeric() %>% as.factor()
    df_spatial %>% 
      ggplot(aes(x = RNA_expr, y = Peak_accessibility, group = Zone, color = Zone)) + 
      geom_point() + 
      ylab("ATAC peak accessibility (Spatial)") + 
      xlab("RNA expression (Spatial)") + 
      ggtitle(sprintf("Spatial r=%s", Spatial_r)) + 
      theme_minimal() -> p_3
    
    patchwork::wrap_plots(p_2,p_3, ncol = 1) -> p_2_3
    patchwork::wrap_plots(p_1, p_2_3, ncol = 2) -> p
    return(p)
  }
  
  build_plot_tripod <- function(id) {
    TRIPOD_analysis %>% dplyr::filter(Plot == id) -> df_
    df_$TF -> TF_
    df_$Gene -> Target_gene_
    df_$Peak -> CRE_
    plot_trios(TF_ = TF_, CRE_ = CRE_, target_gene_ = Target_gene_) -> p
    return(p)
  }
  
  # ---- Dispatcher: pick the right builder for (dataset, plot id) ----
  make_plot <- function(ds, id) {
    switch(ds,
           "cell"       = build_plot_cell(id),
           "spatiotemp" = build_plot_spatiotemp(id),
           "tripod"     = build_plot_tripod(id),
           build_plot_cell(id)
    )
  }
  
  # ================== Click -> modal with plot (cached) ==================
  observeEvent(input$plot_click, {
    req(input$plot_click$id, input$plot_click$ds)
    
    showModal(modalDialog(
      title = paste("Plot:", input$plot_click$id, "—", toupper(substr(input$plot_click$ds,1,1)), substring(input$plot_click$ds,2)),
      plotOutput("plot_modal", height = 520),
      easyClose = TRUE, size = "l"
    ))
    
    output$plot_modal <- renderPlot({
      make_plot(input$plot_click$ds, input$plot_click$id)
    }) %>% bindCache(input$plot_click$ds, input$plot_click$id)  # cache by (dataset, id)
  })
}