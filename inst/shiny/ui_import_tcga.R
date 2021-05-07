shinyPanelTCGA <- fluidPage(
  fluidRow(hr(),
           box(width = 12,
    useShinyalert(),
    add_busy_spinner(spin = "fading-circle"),
    uiOutput(outputId = "tcga_tumor"),
    bsTooltip("tcga_tumor", "Select tumors you want to download from this list", placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("import_tcga", "Import"),
    tags$a(href="https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations", "Full tumor name list"),
    #div(dataTableOutput("tcga_contents"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;")
  )
  )
)
