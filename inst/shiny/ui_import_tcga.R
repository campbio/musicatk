shiny_panel_tcga <- fluidPage(
            tags$head(
              tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
              tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
            ),  
  #Adding box formatting
           div(style = "position: relative;",
             box(p("Select desired TCGA dataset"),
             width = 6,
    uiOutput(outputId = "tcga_tumor"),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "Select desired TCGA dataset, click import button to download them
                            It may take some time to download and process all of the selected datasets. 
                            Once finished, a completion notification will appear at the bottom",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
    #Adding help tooltips
    bsTooltip("tcga_tumor",
              "Select tumors you want to download from this list",
              placement = "bottom", trigger = "hover",
              options = NULL),
    actionButton("import_tcga", "Import"),
    hr(),
    div(DT::DTOutput("tcga_contents"),
        style = "height:500px;
                  overflow-y: scroll;overflow-x: scroll;"),
  )
  )
)

