shiny_panel_tcga <- fluidPage(
  #Adding box formatting
           box(width = 6,
          #Adding hyperlink to TCGA main page
          tags$a(href = "https://gdc.cancer.gov/resources-tcga-
           users/tcga-code-tables/tcga-study-abbreviations",
                      "Full tumor name list"),
    uiOutput(outputId = "tcga_tumor"),
    #Adding help tooltips
    bsTooltip("tcga_tumor",
              "Select tumors you want to download from this list",
              placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("import_tcga", "Import"),
    #Adding hyperlink to TCGA main page
    tags$a(href = "https://gdc.cancer.gov/resources-tcga-
           users/tcga-code-tables/tcga-study-abbreviations",
           "Full tumor name list")
  )
)
