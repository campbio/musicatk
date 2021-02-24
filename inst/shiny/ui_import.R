shinyPanelImport <- fluidPage(
          fluidRow(
            column(width = 12,
                    h3("Step 1: Select File")
            )
          ),
          
          fileInput("file", "File:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          actionButton("upload", "Upload"),
          h3("Samples"),
          #tags$div(id = "file_id",dataTableOutput("my_file_name"),style =  "font-size:40%"),
          uiOutput('undoUI'),
          tags$div(id = "file_id",DT::dataTableOutput("dtable")),
          
          hr(),
          actionButton("import", "Import"),
          uiOutput("spinner"),
          
          
          downloadButton("download_musica", "Download Musica Variants"),
          actionButton(inputId = "reset", label = "Clear Musica Summary"),
          div(textOutput("musica_contents_summary")),
          hr(),
          div(dataTableOutput("musica_contents"))
)
          
