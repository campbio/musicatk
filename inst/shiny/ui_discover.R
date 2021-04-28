shinyPanelDiscover <- fluidPage(
  box(width = 12,
    uiOutput("DiscoverTable"),
    hidden(textInput("GetTableName", "Custon Table Name")),
    textInput("NumberOfSignatures", "Number of signatures", value = 5),
    radioButtons("Method", "Method",
                 choices = list("Latent Dirichlet Allocation (lda)" = "lda", 
                                "Negative Matrix Factorization (nmf)" = "nmf")),
    #textInput("Seed", h3("Seed")),
    textInput("nStart", "Number of random starts", value = 10),
    # textInput("MusicaResultName", "Name for musica result object", value = "SBS96-Result"),
    uiOutput("DiscoverResultName"),
    textOutput("DiscoverWarning"),
    actionButton("DiscoverSignatures", "Discover Signatures"),
    bsTooltip("DiscoverSignatures",
              "Create a musica result object that contains the signatures and exposures of each sample.", 
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("Method",
              "Method to use for mutational signature discovery.", 
              placement = "left", trigger = "hover", options = NULL),
    bsTooltip("NumberOfSignatures",
              "Number of signatures to discover.", 
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("GetTableName",
              "Name of the table to use for signature discovery.", 
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("nStart",
              "Number of independent random starts used in the mutatinoal signature algorithms.", 
              placement = "bottom", trigger = "hover", options = NULL),
    shinybusy::use_busy_spinner(spin = "double-bounce")
    
  )
)

