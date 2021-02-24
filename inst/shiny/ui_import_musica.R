shinyPanelResult <- fluidPage(
  fluidRow(
    wellPanel(id = "musica_data",
              #h1("Musica Object Summary"),
              fileInput("musica_file", "Upload Musica Result Object:",
                        multiple = TRUE,
                        accept = c(".rda","rds")),
              helpText("Use this option only if you don't want to upload your own files in Step 1 for generating
                             a musica object."),
              hr(),
              downloadButton("download_musica_result", "Download Musica Variants"),
              actionButton(inputId = "reset_musica", label = "Clear Musica Summary"), 
              div(textOutput("musica_upload_summary")),
              hr(),
              div(dataTableOutput("musica_upload"))
    
    #actionButton("get_musica", "Get Musica Data"),
    #actionButton("MusicaResults","Get Musica Result Object")
    
    
  )
  ))
  
