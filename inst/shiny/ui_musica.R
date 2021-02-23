shinyPanelMusica <- fluidPage(
  fluidRow(
    uiOutput("genome_list"),
    textOutput("genome_select"),
    column(width = 12,
           h3("Step 3: Create Musica Object"),
           checkboxInput("ref_chr", "Check Reference Chromosomes", TRUE),
           checkboxInput("ref_bases", "Check Reference Bases", TRUE),
           checkboxInput("convert_dbs", "Convert DBS", TRUE),
           checkboxInput("stand_indels", "Standardize Indels", TRUE),
           actionButton("get_musica_object", "Get Musica Object"),
           helpText("Please find your uploaded object summary on the Import tab.")
    )),
    
    
)