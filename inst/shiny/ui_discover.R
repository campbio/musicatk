shinyPanelDiscover <- fluidPage(
  box(
    selectInput("SelectTable", h3("Select Count Table"),
                choices = list("SBS96" = 1, "SBS192" = 2, 
                               "DBS" = 3, "Indel" = 4), 
                selected = 1),
    textInput("NumberOfSignatures", h3("Number of signatures")),
    
    radioButtons("Method", h3("Method"),
                 choices = list("LDA" = 1, "NMF" = 2)),
    textInput("Seed", h3("Seed")),
    textInput("nStart", h3("Number of random starts")),
    actionButton("MusicaResults", h3("Discover Signatures"))
  )
)
