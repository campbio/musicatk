shinyPanelImport <- fluidPage(
          fluidRow(
            column(width = 12,
                    h3("Step 1: Select File")
            )
          ),
          
          fileInput("file", "File:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          h3("Samples"),
          #dataTableOutput("deletable"),
          tags$div(id = "file_id",dataTableOutput("my_file_name"),style =  "font-size:40%"),
          actionButton("import", "Import"),
          uiOutput("spinner"),
          hr(),
          
          wellPanel(id = "musica_data",
                    #h1("Musica Object Summary"),
                    fileInput("musica_file", "Upload Musica Object:",
                              multiple = TRUE,
                              accept = c(".rda","rds")),
                    hr(),
                    downloadButton("download_musica", "Download Musica Variants"),
                    actionButton(inputId = "reset", label = "Clear Musica Summary"),
                    div(textOutput("musica_contents_summary"))),
                    hr(),
                    div(dataTableOutput("musica_contents")),
                    div(tableOutput("musica_upload"))
          
          #actionButton("get_musica", "Get Musica Data"),
          #actionButton("MusicaResults","Get Musica Result Object")
          
        
)
