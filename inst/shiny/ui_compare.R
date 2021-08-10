shiny_panel_compare <- fluidPage(
  box(width = 6,
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
              placement = "bottom", trigger = "hover", options = NULL)
  ),
  box(width = 12,
       uiOutput("download_comparison"),
       dataTableOutput("compare_table")
       )
)
