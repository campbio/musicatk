library(musicatk)

source("server_tables.R", local = T)
source("server_discover.R", local = T)

server <- function(input, output, session) {
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
