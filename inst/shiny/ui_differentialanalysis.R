shinyPanelDifferentialAnalysis <- fluidPage(
  box(width = 12,
      helpText("Use this tab to run differential analysis on the 
               signature exposures of annotated samples within the 
               musica_result object"),
      uiOutput("DiffAnalResult"),
      uiOutput("DiffAnalAnnot"),
      radioButtons("DiffMethod", label = "Method", 
                   choices = c("wilcoxon", "kruskal", "glm.nb")),
      uiOutput("DiffAnalGroups"),
      actionButton("RunDiffAnal", "Run Differential Analysis"),
      dataTableOutput("DiffTable")
      )
)