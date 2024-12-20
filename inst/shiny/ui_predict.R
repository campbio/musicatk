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
    uiOutput("predict_table"),                                     
    uiOutput("predicted_result"),
    uiOutput("predicted_signatures"),
    radioButtons("predict_algorithm", "Algorithm",
                 choices = list("Latent Dirichlet Allocation (lda)" = "lda",
                                "decompTumor2Sig")),
    #hidden(selectInput("predict_genome_list", "Reference genome:",
    #                   list("hg38", "hg19", "mm9", "mm10"),
    #                   width = "100%")),
    uiOutput("predict_result_name"),
    uiOutput("predict_model_name"),
    textOutput("predict_warning"),
    actionButton("predict_sigs", "Predict Signature Exposures"),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "Exposures for samples can be predicted using an existing set of signatures stored in a result_model object. Algorithms available for prediction are “lda” and “decompTumor2Sig”.
           The “Signatures to Predict” dropdown contains all the signatures in the result_model object selected from the “Result to Predict” dropdown. 
           You can search this dropdown and select multiple signatures.",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
    #bsTooltip("predict_table",
    #          "Select the modality.",
    #          placement = "right", trigger = "hover", options = NULL),
    #bsTooltip("cosmic_SBS_Sigs",
    #          "Which Cosmic signature to use",
    #          placement = "right", trigger = "hover", options = NULL),
    #bsTooltip("cosmic_DBS_sigs",
    #          "Which Cosmic signature to use",
    #          placement = "right", trigger = "hover", options = NULL),
    #bsTooltip("cosmic_INDEL_Sigs",
    #          "Which Cosmic signature to use",
    #          placement = "right", trigger = "hover", options = NULL),
    #bsTooltip("predict_algorithm",
    #          "Algorithm to use for prediction of exposures",
    #          placement = "bottom", trigger = "hover", options = NULL),
    #bsTooltip("predict_cosmic",
    #          "Update the musica object tho contain the predicted
    #          exposures for samples using an existing set of signatures.",
    #          placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("predict_algorithm",
              "Algorithm for prediction.",
              placement = "right", trigger = "hover", options = NULL)
  ))
)
