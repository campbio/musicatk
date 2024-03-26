shiny_panel_compare <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),
  div(style = "position: relative;", box(width = 6,
    uiOutput("compare_result_a"),
    uiOutput("compare_result_b"),
    # textInput("Threshold", "Threshold", value = "0.9"),
    textInput("Threshold", "Threshold", value = ".5"),
    radioButtons("compare_metric", "Similarity Metric",
                 choices = c("Cosine" = "cosine",
                             "Jensen-Shannon Divergence (jsd)" = "jsd")),
    uiOutput("compare_validate"),
    actionButton("compare_results", "Compare Results"),
    shinybusy::use_busy_spinner(spin = "double-bounce"),
    bsTooltip("Threshold",
              "Treshold for similarity",
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("compare_metric",
              "Method for calculating similarity.",
              placement = "left", trigger = "hover", options = NULL),
    bsTooltip("compare_results",
              "Compare two result objects to find similar signatures.",
              placement = "bottom", trigger = "hover", options = NULL),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "A common analysis is to compare discovered signatures estimated in a dataset to those generated in other datasets or to those in the COSMIC database. 
           The threshold acts as a cutoff similarity score and can be any value between 0 and 1. Results will populate in a downloadable data table.",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;")   
    
  )),
  box(width = 12,
       uiOutput("download_comparison"),
       dataTableOutput("compare_table")
       )
)
