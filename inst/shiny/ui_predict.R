shiny_panel_predict <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),  
  div(style = "position: relative;", box(width = 6,
    uiOutput("predicted_result"),
    uiOutput("precited_signatures"),
    uiOutput("predict_table"),
    radioButtons("predict_algorithm", "Algorithm",
                 choices = list("Latent Dirichlet Allocation (lda)" = "lda",
                                "decompTumor2Sig",
                                "deconstructSigs")),
    hidden(selectInput("predict_genome_list", "Reference genome:",
                       list("hg38", "hg19", "mm9", "mm10"),
                       width = "100%")),
    uiOutput("predict_result_name"),
    textOutput("predict_warning"),
    actionButton("predict_sigs", "Predict Signature Exposures"),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "Exposures for samples can be predicted using an existing set of signatures stored in a musica result object. Algorithms available for prediction include “lda”, “decompTumor2Sig”, and “deconstructSigs”.
           Note that “deconstructSigs” can only work with variants stored in hg19 format. The “Signatures to Predict” dropdown contains all the signatures in the result object selected from the “Result to Predict” dropdown. 
           You can search this dropdown and select multiple signatures.",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
    bsTooltip("cosmic_count_table",
              "Select the Cosmic object containing mutational
              signatures to predict.",
              placement = "right", trigger = "hover", options = NULL),
    bsTooltip("cosmic_SBS_Sigs",
              "Which Cosmic signature to use",
              placement = "right", trigger = "hover", options = NULL),
    bsTooltip("cosmic_DBS_sigs",
              "Which Cosmic signature to use",
              placement = "right", trigger = "hover", options = NULL),
    bsTooltip("cosmic_INDEL_Sigs",
              "Which Cosmic signature to use",
              placement = "right", trigger = "hover", options = NULL),
    bsTooltip("predict_algorithm",
              "Algorithm to use for prediction of exposures",
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("predict_cosmic",
              "Create a musica result object that contains the predicted
              exposures for samples using an existing set of signatures.",
              placement = "bottom", trigger = "hover", options = NULL)
  ))
)
