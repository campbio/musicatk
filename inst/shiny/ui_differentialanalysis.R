shinyPanelDifferentialAnalysis <- fluidPage(
  box(width = 12,
      helpText("Use this tab to run differential analysis on the 
               signature exposures of annotated samples within the 
               musica_result object"),
      uiOutput("DiffAnalResult"),
      uiOutput("DiffAnalAnnot"),
      radioButtons("DiffMethod", label = "Method", 
                   choices = list("Wilcoxon Rank Sum Test" = "wilcox",
                                  "Kruskal-Wallis Rank Sum Test" = "kruskal",
                                  "Negative Binomial Regression" = "glm.nb")),
      uiOutput("DiffAnalGroups"),
      textOutput("DiffError"),
      actionButton("RunDiffAnal", "Run Differential Analysis")
      ),
  box(width = 12,
      downloadButton("DownloadDiff", "Download"),
      bsTooltip("DownloadDiff",
                "Download the differential exposure table",
                placement = "bottom", trigger = "hover", options = NULL),
      dataTableOutput("DiffTable")
      )
)