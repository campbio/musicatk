shinyPanelTables <- fluidPage(
  box(
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96", "SBS192",
                               "DBS", "Indel",
                               "Custom"),
                selected = 1),
    hidden(selectInput("StrandType", h3("Strand Type"), 
                  choices = list("", "Transcript_Strand", "Replication_Strand"),
                  selected = NULL)),
    radioButtons("OverwriteTable", h3("Overwrite Existing Tables"), 
                 choices = list("No" = F, "Yes" = T)),
    #textInput("TableName", h3("Table Name")),
    actionButton("AddTable", h3("Add Table"))
  )
)
