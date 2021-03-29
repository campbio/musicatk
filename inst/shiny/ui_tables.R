shinyPanelTables <- fluidPage(
  box(width = 12,
    # uiOutput("TableMusicaList"),
    helpText("Use this tab if your musica object has no mutational count table,
    or if you wish to add additional count tables."),
    selectInput("SelectTable", "Select Count Table",
                choices = list("SBS96", "SBS192 - Transcript_Strand", 
                "SBS192 - Replication_Strand", "DBS", "Indel"),
                selected = 1),
    uiOutput("TableGenomeList"),
    textOutput("TableGenomeWarning"),
    uiOutput("AllowTable"),
    shinybusy::use_busy_spinner(spin = "double-bounce"),
    bsTooltip("SelectTable",
              "Name of the standard table to build.", 
              placement = "bottom", trigger = "hover", options = NULL)
  ),
  uiOutput("CombineTable"),
  shinybusy::use_busy_spinner(spin = "double-bounce")
  
)
