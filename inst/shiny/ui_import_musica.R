shinyPanelResult <- fluidPage(
  fluidRow(
    useShinyjs(),
    #wellPanel(id = "musica_data",
              #h1("Upload your Music Result Object"),
              radioButtons("musica_button","Upload Musica Result or Object",c("Musica Result" = "result",
                                                                           "Musica Object" = "object")),
              
              fileInput("musica_file", "Upload Musica Result or Object:",
                        multiple = TRUE,
                        accept = c(".rda","rds")),
              hr(),
              textInput("MusicaResultName", h3("Name your musica result object:")),
              actionButton(inputId = "upload_musica", label = "Upload"),
              hr(),
              downloadButton("download_musica_result", "Download Musica Variants"),
              actionButton(inputId = "reset_musica", label = "Clear Musica Summary"), 
              div(textOutput("musica_upload_summary")),
              hr(),
              div(dataTableOutput("musica_upload")),
              #bsTooltip("musica_data", "Upload your own Musica RResult Object in .rda or .rds format.", placement = "bottom", trigger = "hover",
              #options = NULL),
              bsTooltip("upload_musica", "Press button to upload your own Musica Result Object.", placement = "bottom", trigger = "hover",
              options = NULL),
              bsTooltip("download_musica_result", "Press button to download the musica variants.", placement = "bottom", trigger = "hover",
              options = NULL)
    
    #actionButton("get_musica", "Get Musica Data"),
    #actionButton("MusicaResults","Get Musica Result Object")
    
    
  ))
  
  
