shinyPanelTCGA <- fluidPage(
  fluidRow(box(
    shinyalert::useShinyalert(),
    shinybusy::add_busy_spinner(spin = "fading-circle"),
    uiOutput(outputId = "tcga_tumor"),
    actionButton("import_tcga", "Import")
  )
  )
)