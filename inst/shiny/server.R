library(musicatk)

server <- function(input, output) {
  observeEvent(input$get_musica, {
    maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  })
  observeEvent(input$get_musica_result,{
    musica <- data(res_annot)
    
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
  observeEvent(input$overwriteTable, { 
    overwrite <- T
  })
  observeEvent(input$keepTable, { 
    overwrite <- F 
  })
  observeEvent(input$AddTable, {
    # Currenly only using seelct genomes.
    g <- select_genome(input$Genome)
    overwrite = F
    if (input$SelectTable != "Custom") {
      # Check if table already exists
      if(input$SelectTable %in% names(extract_count_tables(musica))) {
        print(input$SelectTable)
        print(input$Genome)
        req(shinyalert::shinyalert(text = "Table already exists. 
                                   Do you want to overwrite the existing table?", 
                                   showCancelButton = T,
        callbackR = function(x) { if (x == T) {
          # build_standard_table(musica, g = g, 
          #                      table_name =input$SelectTable,
          #                      overwrite = T)}
          overwrite = T}}))
      }
      print(overwrite)
      build_standard_table(musica, g = g, table_name = input$SelectTable,
                           overwrite = overwrite)
      # print(names(extract_count_tables(musica)))
      # View(extract_count_tables(musica))
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
