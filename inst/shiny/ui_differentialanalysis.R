shiny_panel_diffanal <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),
  div(box(width = 12,
      uiOutput("diff_anal_result"),
      uiOutput("diff_anal_annot"),
      radioButtons("diff_method", label = "Method",
                   choices = list("Wilcoxon Rank Sum Test" = "wilcox",
                                  "Kruskal-Wallis Rank Sum Test" = "kruskal",
                                  "Negative Binomial Regression" = "glm.nb")),
      uiOutput("diff_anal_groups"),
      textOutput("diff_error"),
      actionButton("run_diff_anal", "Run Differential Analysis"),
      tags$a(href = "#", 
             tags$i(class = "fas fa-question-circle"),
             title = "Need help?", 
             `data-toggle` = "popover", 
             `data-trigger` = "focus", 
             `data-content` = "The “Exposure Differential Analysis” tab is used to run 
             differential analysis on the signature exposures of annotated samples. There 
             are 3 methods to perform the differential analysis: Wilcoxon Rank Sum Test,
             Kruskal-Wallis Rank Sum Test, and a negative binomial regression (glm).",
             `data-html` = "true",
             `data-placement` = "left",
             style = "position: absolute; top: 5px; right: 5px; cursor: pointer;")        
      )),
  box(width = 12,
      downloadButton("download_diff", "Download"),
      bsTooltip("download_diff",
                "Download the differential exposure table",
                placement = "bottom", trigger = "hover", options = NULL),
      dataTableOutput("diff_table")
      )
)
