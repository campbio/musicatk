shinyPanelImport <- fluidPage(
  
      
        box(
          fluidRow(
            column(width = 12,
                    h3("Step 1: Select File")
            )
          ),
          
          fileInput("file", "File:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          actionButton("get_musica", "Get Musica Data"),
          actionButton("MusicaResults","Get Musica Result Object"),
          
          )
        
)