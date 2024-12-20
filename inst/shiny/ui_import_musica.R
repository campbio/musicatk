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
  fluidRow(box(width = 6, actionButton("example", "Import example data"))),
  fluidRow(div(style = "position: relative;",box(width = 6,
              #radioButtons("musica_button",
              #            "Upload Musica Result or Object",
              #             c("Musica Result" = "result",
              #              "Musica Object" = "object")),
              tags$a(href = "#", 
                     tags$i(class = "fas fa-question-circle"),
                     title = "Need help?", 
                     `data-toggle` = "popover", 
                     `data-trigger` = "focus", 
                     `data-content` = " import a user’s previously generated musica object in .rda or .rds format.The browse button can be selected for choosing files. 
                     Once selected, a bar saying file upload is completed will appear. 
                     Once finished, a notification message in the bottom right will appear along with a table of variants. Users can then skip the Create Musica tab and move to the build tables tab directly.",
                     `data-html` = "true",
                     `data-placement` = "left",
                     style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),              
              #Adding file upload feature
              fileInput("musica_file",
                        "Upload Musica Object:",
                        multiple = TRUE,
                        accept = c(".rda", "rds", ".RDS")),
              hr(),
              #uiOutput("musica_result_name"),
              actionButton(inputId = "upload_musica", label = "Upload"),
              hr(),
              #adding download feature
              downloadButton("download_musica_result", "Download Variants"),
              div(textOutput("musica_upload_summary")),
              hr(),
              div(DT::DTOutput("musica_upload"),
                  style = "height:500px;
                  overflow-y: scroll;overflow-x: scroll;"),
              #Adding help tooltips
              bsTooltip("upload_musica",
                        "Press button to upload your own Musica object.",
                        placement = "bottom", trigger = "hover",
              options = NULL),
              bsTooltip("download_musica_result",
                        "Press button to download the musica variants.",
                        placement = "bottom", trigger = "hover",
              options = NULL)
))))
