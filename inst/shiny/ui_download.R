shiny_panel_download <- fluidPage(
  hr(),
  #Adding box formatting
  fluidRow(box(uiOutput(outputId = "select_mus_obj_download"),
               #Adding download buttons
               downloadButton("download_mus_obj", "Download Musica Object"),
               hr(),
               #uiOutput(outputId = "select_res_download"),
               #downloadButton("download_res", "Download Musica Result")
)))
