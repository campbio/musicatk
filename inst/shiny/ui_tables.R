shinyPanelTables <- fluidPage(
  box(width = 12,
    # uiOutput("TableMusicaList"),
    helpText("Use this tab if you have not created a mutational count table,
    or if you wish to add additional count tables."),
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96", "SBS192 - Transcript_Strand", 
                "SBS192 - Replication_Strand", "DBS", "Indel"),
                selected = 1),
    uiOutput("AllowTable"),
    shinybusy::use_busy_spinner(spin = "double-bounce")
    
  ),
  uiOutput("CombineTables"),
  shinybusy::use_busy_spinner(spin = "double-bounce")
  
)
