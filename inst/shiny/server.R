library(musicatk)


options(shiny.maxRequestSize = 100*1024^2)
source("server_tables.R", local = T)
source("server_discover.R", local = T)

###################### Zainab's Code ##########################################
server <- function(input, output) {
  observeEvent(input$get_musica, {
    maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  })
  
  observeEvent(input$MusicaResults,{
    data(res_annot)
  })
  
  variants <- reactive({
    req(input$file)
    file_name <- input$file$datapath
    var <- extract_variants_from_maf_file(maf_file = file_name)
    return(var)
  })
  genome <- reactive({
    gen <- input$GenomeSelect
    gen <- select_genome(gen)
    return(gen)
  })
  output$genome_select <- renderText({
    paste("Genome selected:", input$GenomeSelect)
  })
  musica_contents <- eventReactive(input$get_musica_object,{
    musica <- create_musica(x = variants(), genome = genome())
    return(musica)
  })
  
  
  output$musica_contents <- renderTable({
    return(head(musica_contents()@variants))
    shinyjs::show(id="musica_contents")
    js$enableTabs();
  })
  output$download_musica <- downloadHandler(
    filename = function() {
      paste("musica_variants", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(musica_contents()@variants, file, row.names = FALSE)
    }
  )
  
    

  
  
  
        
  
    
###############################################################################  
    
    
###################### Nathan's Code ##########################################
  # Create dynamic table
  vals <- reactiveValues(
    musica = musica,
    musica_result = NULL
  ) 
  
  output$DiscoverTable <- renderUI({
    tagList(
      selectInput("SelectDiscoverTable", h3("Select Count Table"),
                  choices = names(extract_count_tables(vals$musica)))
    )
  })
  observeEvent(input$AddTable, {
    if(input$SelectTable %in% names(extract_count_tables(vals$musica))) {
      showModal(modalDialog(
        title = "Existing Table.",
        "Do you want to overwrite the existing table?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirmOverwrite", "OK"),
          modalButton("Cancel"))
        ))
    } else{
      add_tables(input, vals)
    }})
  
  sigFun <- renderPrint({
    cosmic_v2_subtype_map(input$CosmicV2Signatures)
  })
  # observeEvent(input$CosmicV2Signatures, {
    output$CosmicSignatures <- renderUI({
      tagList(
      textInput("CosmicSignatures", h3("Select Cosmic signatures to predict."),
                  value = sigFun)
    )
    })
    # })
  
  observeEvent(input$confirmOverwrite, {
    removeModal()
    add_tables(input, vals)
  })
  
  
  # Initially hidden additional required option for SBS192 and Custom
  # table creation
  observeEvent(input$SelectTable, {
    if (input$SelectTable == "SBS192") {
      show(id = "StrandType")
      hide(id = "GetTableName")
    } else if (input$SelectTable == 5) {
      hide(id = "StrandType")
      hide(id = "GRangeFile")
      show(id = "GetTableName")
    } else {
      hide(id = "GetTableName")
      hide(id = "StrandType")
      hide(id = "GRangeFile")
    }
  })

  observeEvent(input$StrandType, {
    if (input$StrandType == "Replication_Strand") {
      show(id = "GRangeFile")
    } else {
      hide(id = "GRangeFile")
    }
  })
  

  # Test when musica code has been generated
  observeEvent(input$MusicaResults, {
    vals$musica_result <- discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      method = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$nStart))
  })
###############################################################################
  
}
