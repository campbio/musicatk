shiny_panel_discover <- fluidPage(
  box(width = 6,
    # UI to choose counts table
    uiOutput("discover_table"),
    hidden(textInput("get_table_name", "Custon Table Name")),
    textInput("number_of_signatures", "Number of signatures", value = 5),
    radioButtons("Method", "Method",
                 choices = list("Latent Dirichlet Allocation (lda)" = "lda",
                                "Negative Matrix Factorization (nmf)" = "nmf")),
    textInput("n_start", "Number of random starts", value = 10),
    # UI to set result name. Default name of the counts table.
    uiOutput("discover_result_name"),
    textOutput("discover_warning"),
    actionButton("discover_signatures", "Discover Signatures"),
    bsTooltip("discover_signatures",
              "Create a musica result object that contains the signatures
              and exposures of each sample.",
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("Method",
              "Method to use for mutational signature discovery.",
              placement = "left", trigger = "hover", options = NULL),
    bsTooltip("number_of_signatures",
              "Number of signatures to discover.",
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("get_table_name",
              "Name of the table to use for signature discovery.",
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("n_start",
              "Number of independent random starts used in the mutatinoal
              signature algorithms.",
              placement = "bottom", trigger = "hover", options = NULL),
    shinybusy::use_busy_spinner(spin = "double-bounce")
  )
)
