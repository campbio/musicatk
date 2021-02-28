shinyPanelTables <- fluidPage(
  box(
    # uiOutput("TableMusicaList"),
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96", "SBS192 - Transcript_Strand", 
                "SBS192 - Replication_Strand", "DBS", "Indel"),
                selected = 1),
    uiOutput("AllowTable")
    )
)
