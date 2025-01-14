shiny_panel_import <- fluidPage(
  #fluidRow(box(width = 6, actionButton("example", "Import example data"))),
          fluidRow(
            tags$head(
              tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
              tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
            ),
          div(style = "position: relative;",
            box(
              width = 6,
              div(fileInput("file", "Select file: ",
                  multiple = TRUE,
                  accept = c(".maf", ".vcf", ".txt"))),
                  tags$a(href = "#", 
                     tags$i(class = "fas fa-question-circle"),
                     title = "Need help?", 
                     `data-toggle` = "popover", 
                     `data-trigger` = "focus", 
                     `data-content` = "Select .maf, .vcf, or .txt files to extract variants. You can select multiple files.
                                        Clicking the ‘Add Samples’ button will add the files to the list of samples to be imported",
                     `data-html` = "true",
                     `data-placement` = "left",
                     style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),
          br(),
          actionButton("upload", "Add samples")))),
          bsTooltip("upload", "Press button to add your uploaded files to Sample List",
                    placement = "bottom", trigger = "hover", options = NULL),
          bsTooltip("example",
                    "Press button to add an example file",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          fluidRow(box(width = 6, HTML(paste0("<b>","Added files:","</b>")),
          tags$div(id = "file_id", DT::DTOutput("dtable"),
                   style = "overflow-y: scroll;overflow-x: scroll;"),
          uiOutput("undo_ui"),
          hr(),
          uiOutput("spinner"))),
          #Adding the download button
          fluidRow(tags$head(
                  tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
                  tags$script(HTML('
                  $(document).ready(function(){
                  $("[data-toggle=\'popover\']").popover();
                  });
                  '))
                    ),
                    div(box(p("Import data and Download variants"),

                       width = 6, actionButton("import", "Import"),
                       downloadButton("download_musica", "Download Variants"),
                       tags$a(href = "#", 
                              tags$i(class = "fas fa-question-circle"),
                              title = "Need help?", 
                              `data-toggle` = "popover", 
                              `data-trigger` = "focus", 
                              `data-content` = "the ‘Import’ button can be clicked to start the importing process.
                                                variant table can be downloaded by pressing the ‘Download Variants’ button",
                              `data-html` = "true",
                              `data-placement` = "left",
                              style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),
                       hr(),
                       div(DT::DTOutput("musica_contents"),
                           style = "height:500px;
                           overflow-y: scroll;overflow-x: scroll;"),
                       ))),
          #Adding help tootltips

          bsTooltip("file_id", "Table of files that have been added by you",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("import",
                    "Press button to import the files in the Sample List",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("undo_ui",
                    "Press button to undo any deletions
                    you made in the Sample List",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          fluidRow(
            box(
              width = 12,
              HTML('<p> For more information, visit <a href="https://camplab.net/musicatk/current/articles/articles/tutorial_tcga_ui.html">this link</a>.</p>')
            )
          )
)
