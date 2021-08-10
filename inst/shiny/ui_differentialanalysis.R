shiny_panel_diffanal <- fluidPage(
  box(width = 12,
      uiOutput("diff_anal_result"),
      uiOutput("diff_anal_annot"),
      radioButtons("diff_method", label = "Method",
                   choices = list("Wilcoxon Rank Sum Test" = "wilcox",
                                  "Kruskal-Wallis Rank Sum Test" = "kruskal",
                                  "Negative Binomial Regression" = "glm.nb")),
      uiOutput("diff_anal_groups"),
      textOutput("diff_error"),
      actionButton("run_diff_anal", "Run Differential Analysis")
      ),
  box(width = 12,
      downloadButton("download_diff", "Download"),
      bsTooltip("download_diff",
                "Download the differential exposure table",
                placement = "bottom", trigger = "hover", options = NULL),
      dataTableOutput("diff_table")
      )
)
