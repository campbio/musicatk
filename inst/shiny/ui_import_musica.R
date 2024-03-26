shiny_panel_result <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ), 
  #Adding box formatting and creating radio buttons
  fluidRow(div(style = "position: relative;",box(width = 6,
              radioButtons("musica_button",
                           "Upload Musica Result or Object",
                           c("Musica Result" = "result",
                            "Musica Object" = "object")),
              tags$a(href = "#", 
                     tags$i(class = "fas fa-question-circle"),
                     title = "Need help?", 
                     `data-toggle` = "popover", 
                     `data-trigger` = "focus", 
                     `data-content` = " import a user’s previously generated musica variant object or a musica result object in .rda or .rds format.The browse button in the musica result or musica result object can be selected for selecting files. 
                     Once selected, a bar saying file upload is completed will appear. By default, the musica result objects will be named after the file’s name but it can be changed in the text box labeled “Name your musica result object”. 
                     Once finished, a notification message in the bottom right will appear along with a table of variants. If a musica variant object was uploaded then the users can skip the Create Musica tab and move to the build tables tab directly. 
                     If a musica result object was uploaded, then all the downstream tabs will be accessible.",
                     `data-html` = "true",
                     `data-placement` = "left",
                     style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),              
              #Adding file upload feature
              fileInput("musica_file",
                        "Upload Musica Result or Object:",
                        multiple = TRUE,
                        accept = c(".rda", "rds", ".RDS")),
              hr(),
              uiOutput("musica_result_name"),
              actionButton(inputId = "upload_musica", label = "Upload"),
              hr(),
              #adding download feature
              downloadButton("download_musica_result", "Download Variants"),
              div(textOutput("musica_upload_summary")),
              hr(),
              div(dataTableOutput("musica_upload"),
                  style = "height:500px;
                  overflow-y: scroll;overflow-x: scroll;"),
              #Adding help tooltips
              bsTooltip("upload_musica",
                        "Press button to upload your own Musica Result Object.",
                        placement = "bottom", trigger = "hover",
              options = NULL),
              bsTooltip("download_musica_result",
                        "Press button to download the musica variants.",
                        placement = "bottom", trigger = "hover",
              options = NULL)
))))
