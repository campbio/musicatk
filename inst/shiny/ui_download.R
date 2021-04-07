shinyPanelDownload <- fluidPage(
  hr(),
  fluidRow(box(uiOutput(outputId = "select_mus_obj_download"),
               downloadButton("download_mus_obj", "Download Musica Object"),
               hr(),
               uiOutput(outputId = "select_res_download"),
               downloadButton("download_res", "Download Musica Reesult")
  
)))