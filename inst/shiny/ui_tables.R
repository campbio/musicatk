shinyPanelTables <- fluidPage(
  box(
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96", "SBS192",
                               "DBS", "Indel",
                               "Custom"),
                selected = 1),
    fluidRow(
      column(width = 12,
        textInput("Genome", h3("Select Genome"), value = "hg38")
      )
    ),
    #textInput("TableName", h3("Table Name")),
    actionButton("AddTable", h3("Add Table"))
  )
)