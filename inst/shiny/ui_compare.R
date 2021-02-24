shinyPanelCompare <- fluidPage(
  box( 
    uiOutput("CompareResultA"),
    uiOutput("CompareResultB"),
    # textInput("Threshold", "Threshold", value = "0.9"),
    sliderInput("Threshold", "Threshold", min = 0, max = 1, .5),
    radioButtons("CompareMetric", h2("Similarity Metric"), 
                 choices = c("Cosine" = "cosine", 
                             "Jensen-Shannon Divergence (jsd)" = "jsd")),
    actionButton("CompareResults", "CompareResults")
    # plotOutput("ComparePlot")
  )
)