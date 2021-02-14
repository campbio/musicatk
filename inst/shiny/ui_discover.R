shinyPanelDiscover <- fluidPage(
  box(
    # selectInput("SelectTable", h3("Select Count Table"),
    #             choices = list("SBS96" = 1, "SBS192" = 2, 
    #                            "DBS" = 3, "Indel" = 4,
    #                            "Custom" = 5), selected = 1),
    uiOutput("DiscoverTable"),
    hidden(textInput("GetTableName", h3("Custon Table Name"))),
    textInput("NumberOfSignatures", h3("Number of signatures"), value = 5),
    radioButtons("Method", h3("Method"),
                 choices = list("lda", "nmf")),
    #textInput("Seed", h3("Seed")),
    textInput("nStart", h3("Number of random starts"), value = 10),
    textInput("MusicaResultName", h3("Name for musica result object")),
    actionButton("DiscoverSignatures", h3("Discover Signatures"))
    
  )
)

