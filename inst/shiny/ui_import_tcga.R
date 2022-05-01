shiny_panel_tcga <- fluidPage(
  #Adding box formatting
           box(width = 6,
    uiOutput(outputId = "tcga_tumor"),
    #Adding help tooltips
    bsTooltip("tcga_tumor",
              "Select tumors you want to download from this list",
              placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("import_tcga", "Import"),
  )
)
