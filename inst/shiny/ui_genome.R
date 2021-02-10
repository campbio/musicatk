shinyPanelGenome <- fluidPage(
  fluidRow(
    selectInput("GenomeSelect", "Step 2: Choose genome:",
                           list( "Common genomes" = list("hg38","hg17","hg18","hg19","mm8","mm9","mm10")
                                 #, "Genomes" = gg
                                 ), width ='100%'
                          
                           
  )
  ),
  textOutput("genome_select"),
  #actionButton("custom_genome", "Upload custom Genome")
  
)