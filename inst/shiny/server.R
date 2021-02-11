library(musicatk)


options(shiny.maxRequestSize = 100*1024^2)
source("server_tables.R", local = T)

###################### Zainab's Code ##########################################
server <- function(input, output) {
  observeEvent(input$get_musica, {
    maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  })
  
  # observeEvent(input$MusicaResults,{
  #   data(res_annot)
  # })
  
  variants <- reactive({
    req(input$file)
    file_name <- input$file$datapath
    var <- extract_variants_from_maf_file(maf_file = file_name)
    return(var)
  })
  output$genome_list <- renderUI({
    g <- BSgenome::available.genomes()
    g <-strsplit(g,",")
    gg <- gsub("^.*?\\.","", g)
    selectInput("GenomeSelect", "Step 2: Choose genome:",
                list( "Common genomes" = gg, 
                      width ='100%')
    )
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
    musica_result = NULL,
    pred_result = NULL
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
  
  cosmic <- data("cosmic_v2_sigs")
  cosmic_dbs <- data("cosmic_v3_dbs_sigs")
  cosmic_indel <- data("cosmic_v3_indel_sigs")
  
  observeEvent(input$DiscoverSignatures, {
    vals$musica_result <- discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      method = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$nStart))
  })
  
  output$PredictTable <- renderUI({
    tagList(
      selectInput("SelectPredTable", "Select Counts Table",
                  choices = names(tables(vals$musica)))
    )
  })
  
  observeEvent(input$CosmicCountTable, {
    if (input$CosmicCountTable == "SBS") {
      show(id = "CosmicSBSSigs")
      hide(id = "CosmicDBSSigs")
      hide(id = "CosmicINDELSigs")
    } else if (input$CosmicCountTable == "DBS") {
      hide(id = "CosmicSBSSigs")
      show(id = "CosmicDBSSigs")
      hide(id = "CosmicINDELSigs")
    } else {
      hide(id = "CosmicSBSSigs")
      hide(id = "CosmicDBSSigs")
      show(id = "CosmicINDELSigs")
    }
  })
  
  observeEvent(input$PredictCosmic, {
    if (input$CosmicCountTable == "SBS") {
      sigs <- input$CosmicSBSSigs
      res <- cosmic_v2_sigs
    } else if (input$CosmicCountTable == "DBS") {
      sigs <- input$CosmicDBSSigs
      res <- cosmic_v3_dbs_sigs
    } else {
      sigs <- input$CosmicINDELSigs
      res <- cosmic_v3_indel_sigs
    }
    vals$pred_res <- predict_exposure(vals$musica, g = genome, 
                     table_name = input$SelectPredTable,
                     signature_res = res,
                     algorithm = input$PredictAlgorithm,
                     signatures_to_use = as.numeric(sigs))
    # if (!is.null(vals$pred_res)) {
    #   show(id = "Threshold")
    #   show(id = "Compare")
    # } else {
    #   hide(id = "Threshold")
    #   hide(id = "Compare")
    # }
  })
  
  observeEvent(input$Compare, {
    tryCatch( {
     output$ComparePlot <- renderPlot({comparisons <- compare_results(vals$musica_result,
                    vals$pred_res,
                    threshold = input$Threshold)})
     }, error = function(cond) {
       shinyalert::shinyalert(title = "Error", text = cond$message)
       return()
     }
    )
  })
  # output$CosmicSignatures <- renderUI({
  #   sigsDBS <- toString(ncol(signatures(cosmic_v3_dbs_sigs)))
  #   #sigsDBS <- 1:ncol(signatures(cosmic_v3_dbs_sigs))
  #   browser()
  #   tagList(
  #     # textInput("CosmicSignatures", h3("Select Cosmic signatures to predict."),
  #     #             value = sigs)
  #     checkboxInput("CosmicDBSSigs", "Cosmic V3 DBS Signatures",
  #                   value = list(sigsDBS)),
  #     
  #     )
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


###############################################################################
  
}
