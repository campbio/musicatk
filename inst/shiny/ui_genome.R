shinyPanelGenome <- fluidPage(
  fluidRow(
    uiOutput("genome_list"),
  textOutput("genome_select"),
  #actionButton("custom_genome", "Upload custom Genome")
  
))
