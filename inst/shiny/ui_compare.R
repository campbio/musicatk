shinyPanelCompare <- fluidPage(
  box( width = 6,
    uiOutput("CompareResultA"),
    uiOutput("CompareResultB"),
    # textInput("Threshold", "Threshold", value = "0.9"),
    textInput("Threshold", "Threshold", value = ".5"),
    radioButtons("CompareMetric", "Similarity Metric", 
                 choices = c("Cosine" = "cosine", 
                             "Jensen-Shannon Divergence (jsd)" = "jsd")),
    uiOutput("CompareValidate"),
    actionButton("CompareResults", "Compare Results"), 
    shinybusy::use_busy_spinner(spin = "double-bounce"),
    bsTooltip("Threshold",
              "Treshold for similarity", 
              placement = "bottom", trigger = "hover", options = NULL),
    bsTooltip("CompareMetric",
              "Method for calculating similarity.", 
              placement = "left", trigger = "hover", options = NULL),
    bsTooltip("CompareResults",
              "Compare two result objects to find similar signatures.", 
              placement = "bottom", trigger = "hover", options = NULL)
  ),
  box( width = 12,
       uiOutput("DownloadComparison"),
       dataTableOutput("CompareTable")
       )
)
