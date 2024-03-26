shiny_panel_discover <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),  
  div(style = "position: relative;", box(width = 6,
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
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "Mutational signatures and exposures are discovered using a Latent Dirichlet Allocation (LDA) or a Non-Negative Matrix Factorization (NMF) algorithms. 
           These algorithms will deconvolute the mutation count matrix into two matrices: 1) a “signature” matrix containing the probability of each mutation type in each sample and 
           2) an “exposure” matrix containing the estimated counts for each signature in each sample. Select a count table, algorithm, number of signatures, and specify the name of 
           the result and then click “Discover signatures”.",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
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
  ))
)
