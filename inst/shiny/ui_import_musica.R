shinyPanelResult <- fluidPage(
  #Adding box formatting
  fluidRow(box(width = 6,
    useShinyjs(),
    shinyalert::useShinyalert(),
    shinybusy::add_busy_spinner(spin = "fading-circle"),#adding spinner ui
    #wellPanel(id = "musica_data",
              #h1("Upload your Music Result Object"),
              #Creating radio buttons 
              radioButtons("musica_button","Upload Musica Result or Object",c("Musica Result" = "result",
                                                                           "Musica Object" = "object")),
              
              #Adding file upload feature
              fileInput("musica_file", "Upload Musica Result or Object:",
                        multiple = TRUE,
                        accept = c(".rda","rds")),
              hr(),
              uiOutput("MusicaResultName"),
              #textInput("MusicaResultName", value = "", h3("Name your musica result object:")),
              actionButton(inputId = "upload_musica", label = "Upload"),
              hr(),
              #adding download feature
              downloadButton("download_musica_result", "Download Variants"),
              #actionButton(inputId = "reset_musica", label = "Clear Musica Summary"), 
              div(textOutput("musica_upload_summary")),
              hr(),
              div(dataTableOutput("musica_upload"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
              #bsTooltip("musica_data", "Upload your own Musica Result Object in .rda or .rds format.", placement = "bottom", trigger = "hover",
              #options = NULL),
              #Adding help tooltips
              bsTooltip("upload_musica", "Press button to upload your own Musica Result Object.", placement = "bottom", trigger = "hover",
              options = NULL),
              bsTooltip("download_musica_result", "Press button to download the musica variants.", placement = "bottom", trigger = "hover",
              options = NULL)
    
    #actionButton("get_musica", "Get Musica Data"),
    #actionButton("MusicaResults","Get Musica Result Object")
    
    
  )))
  
  
