shinyPanelTCGA <- fluidPage(
  #Adding box formatting
  fluidRow(hr(),
           box(width = 6,
    uiOutput(outputId = "tcga_tumor"),
    #Adding help tooltips
    bsTooltip("tcga_tumor", "Select tumors you want to download from this list", placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("import_tcga", "Import"),
    #Adding hyperlink to TCGA main page
    tags$a(href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations", "Full tumor name list"),
    #div(dataTableOutput("tcga_contents"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;")
  )
  )
)
