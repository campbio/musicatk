shinyPanelTables <- fluidPage(
  box(
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96" = 1, "SBS192" = 2,
                               "DBS" = 3, "Indel" = 4,
                               "Custom" = 5),
                selected = 1),
    fluidRow(
      column(width = 12,
        textInput("SelectGenome", h3("Select Genome"), value = "hg38")
      )
    ),
    textInput("TableName", h3("Table Name")),
    actionButton("AddTable", h3("Add Table"))
  )
)