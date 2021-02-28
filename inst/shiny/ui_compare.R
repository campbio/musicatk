shinyPanelCompare <- fluidPage(
  box( width = 12,
    uiOutput("CompareResultA"),
    uiOutput("CompareResultB"),
    # textInput("Threshold", "Threshold", value = "0.9"),
    textInput("Threshold", "Threshold", value = ".5"),
    radioButtons("CompareMetric", h2("Similarity Metric"), 
                 choices = c("Cosine" = "cosine", 
                             "Jensen-Shannon Divergence (jsd)" = "jsd")),
    uiOutput("CompareValidate"),
    actionButton("CompareResults", "CompareResults"),
    #uiOutput("ComparisonTable")
    plotOutput("ComparePlot")
    #uiOutput("CompareTable")
  )
)
