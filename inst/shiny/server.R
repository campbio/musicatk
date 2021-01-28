library(musicatk)

source("Tables.R", local = T)

server <- function(input, output) {
  observeEvent(input$get_musica, {
    maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  })
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(input$file_vcf)
    req(input$file_maf)
    if(!is.NULL(input$file_vcf)){
      df <- read.vcfR(input$file_vcf$datapath)
    }
    else if(!is.NULL(input$file_maf)){
      df <- read.maf(input$file_maf$datapath)
    }
    return(df)

    
  })
  
###################### Nathan's Code ##########################################
  observeEvent(input$AddTable, {
    add_tables(input)
  })
  observeEvent(input$SelectTable, {
    if (input$SelectTable == "SBS192") {
      show(id = "StrandType")
    } else {
      hide(id = "StrandType")
    }
  })
  # Test when musica code has been generated
  # observeEvent(input$MusicaResults, {
  #   musica_result <- discover_signatures(
  #     output$musica, table_name = input$SelectTable, 
  #     num_signatures = input$NumberOfSignatures,
  #     method = input$Methods,
  #     seed = input$Seed,
  #     nstart = input$nStart)
  # })
###############################################################################
  
}
