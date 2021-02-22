shinyPanelMusica <- fluidPage(
  fluidRow(
    column(width = 12,
            h3("Step 3: Create Musica Object"),
           actionButton("get_musica_object", "Get Musica Object")
    )),
    
    wellPanel(id = "musica_data",
            h3("Musica Object Summary"),
            fileInput("musica_file", "Upload Musica Object:",
                      multiple = TRUE,
                      accept = ".rda"),
            downloadButton("download_musica", "Download Musica Variants"),
            tableOutput("musica_contents"))
)