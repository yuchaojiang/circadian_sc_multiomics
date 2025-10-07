library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(htmltools)
library(ggplot2)
library(patchwork)

# ================== UI ==================
ui <- page_sidebar(
  title = "Single-Cell Multiomic Analysis of Circadian Rhythmicity in Mouse Liver",
  sidebar = sidebar(
    accordion(
      id = "nav",
      accordion_panel("Cell-type specific analysis", value = "cell"),
      accordion_panel("Hepatocytes CRE linkage analysis (Spatial & Circadian)", value = "spatiotemp"),
      accordion_panel("Hepatocytes TRIPOD analysis", value = "tripod"),
      multiple = FALSE
    )
  ),
  card(
    card_header(textOutput("title")),
    DTOutput("table")
  )
)