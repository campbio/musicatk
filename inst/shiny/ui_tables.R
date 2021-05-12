shinyPanelTables <- fluidPage(
  box(width = 6,
    selectInput("SelectTable", "Select Count Table",
                choices = list("SBS96", "SBS192 - Transcript_Strand", 
                "SBS192 - Replication_Strand", "DBS", "Indel"),
                selected = 1),
    selectInput("TableGenomeList", "Reference genome:",
                list("hg38","hg19", "mm9","mm10"), 
                width ='100%'),
    textOutput("TableGenomeWarning"),
    uiOutput("AllowTable"),
    shinybusy::use_busy_spinner(spin = "double-bounce"),
    bsTooltip("SelectTable",
              "Name of the standard table to build.", 
              placement = "right", trigger = "hover", options = NULL)
  ),
  uiOutput("CombineTable"),
  shinybusy::use_busy_spinner(spin = "double-bounce")
  
)
