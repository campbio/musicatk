shiny_panel_result <- fluidPage(
  #Adding box formatting and creating radio buttons
  fluidRow(box(width = 6,
              radioButtons("musica_button",
                           "Upload Musica Result or Object",
                           c("Musica Result" = "result",
                            "Musica Object" = "object")),
              #Adding file upload feature
              fileInput("musica_file",
                        "Upload Musica Result or Object:",
                        multiple = TRUE,
                        accept = c(".rda", "rds")),
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
)))
