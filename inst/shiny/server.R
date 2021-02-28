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
  musica_contents <- eventReactive(input$get_result_object,{
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
    result_objects = list(),
    annotations = NULL,
    comparison = NULL
  ) 
  
  # rep_range needed for SBS192 replication strand
  data(rep_range)
  # needed to predict cosmic sigs
  data("cosmic_v2_sigs")
  data("cosmic_v3_dbs_sigs")
  data("cosmic_v3_indel_sigs")
  cosmic_objects <- list("cosmic_v2_sigs" = cosmic_v2_sigs, 
                      "cosmic_v3_dbs_sigs" = cosmic_v3_dbs_sigs,
                      "cosmic_v3_indel_sigs" = cosmic_v3_indel_sigs)
  
  output$DiscoverTable <- renderUI({
    tagList(
      selectInput("SelectDiscoverTable", h3("Select Count Table"),
                  choices = names(
                    extract_count_tables(vals$musica)))
    )
  })
  
  output$AllowTable <- renderUI({
    if (!is.null(vals$musica)) {
      actionButton("AddTable", h3("Create Table"))
    } else {
      helpText("You must first create or upload a musica object to generate
               count tables.")
    }
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
      showNotification("Table created.")
    }
  })
  
  observeEvent(input$DiscoverSignatures, {
    if (input$MusicaResultName == "" |
        input$NumberOfSignatures == "" |
        input$nStart == "") {
      output$DiscoverWarning <- renderText({
        validate(
          need(input$MusicaResultName != "",
               'You must provide a name for the new result object.'),
          need(input$NumberOfSignatures != "",
               'You must specify the number of expected signatures.'),
          need(input$nStart != "",
               "Please specify the number of random starts.")
        )
      })
      return ()
    }
    if(input$MusicaResultName %in% names(vals$result_objects)) {
      showModal(modalDialog(
        title = "Existing Result Object.",
        "Do you want to overwrite the existing result object?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirmResultOverwrite", "OK"),
          modalButton("Cancel"))
      ))
    } else {
      getResult(input, vals)
      showNotification(paste0("Musica Result object, ", input$MusicaResultName,
                            ", was successfully generated"))
    }
  })
  
  getResult <- function(input, vals) {
    vals$result_objects[[input$MusicaResultName]] <- discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      method = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$nStart))
  }
  
  observeEvent(input$confirmResultOverwrite, {
    removeModal()
    getResult(input, vals)
    showNotification("Existing result object overwritten.")
  })
  
  output$PredictTable <- renderUI({
    tagList(
      selectInput("SelectPredTable", "Select Counts Table",
                  choices = names(tables(vals$musica)))
    )
  })
  
  output$AnnotationMusicaList <- renderUI({
    tagList(
      selectInput("AnnotationMusicaList", h3("Select Musica Object"),
                  choices = c("musica", names(vals$result_objects),
                  names(vals$result_objects)))
    )
  })
  
  output$DiscoverMusicaList <- renderUI({
    tagList(
      selectInput("DiscoverMusicaList", h3("Select Musica Object"),
                  choices = names(vals$result_objects))
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
    if (input$PredictResultName == "" | is.null(sigs)) {
      output$PredictWarning <- renderText({
        validate(
          need(input$PredictResultName != "",
               'You must provide a name for the new result object.'),
          need(!is.null(sigs), 
               'Please select signatures to predict.')
        )
      })
      return ()
    }
    browser()
    vals$result_objects[[input$PredictResultName]] <-
      predict_exposure(vals$musica, g = genome, 
                     table_name = input$SelectPredTable,
                     signature_res = res,
                     algorithm = input$PredictAlgorithm,
                     signatures_to_use = as.numeric(sigs))
    showNotification(paste0("New result object, ", input$PredictResultName,
                            "was created."))
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
    showNotification("Existing table overwritten.")
  })

  output$annotations <- renderTable({
    file <- input$AnnotationsFile
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext %in% c("txt", "csv"), "Please upload a txt file"))
    vals$annotations <- read.csv(file$datapath, header = input$AnnotationHeader)
    vals$annotations
  })
  
  #Add Annotations to Musica object
  observeEvent(input$AddAnnotation, {
    if (!is.null(vals$musica)) {
      tryCatch( {
      sapply(names(vals$annotations), FUN = function(a) {
          samp_annot(vals$musica, a) <- vals$annotations[[a]]
        })}, error = function(cond) {
          shinyalert::shinyalert(title = "Error", text = cond$message)
          return()
        })
    } else if (!is.null(vals$result_objects[[input$AnnotationMusicaList]])) {
      tryCatch( {
        sapply(names(vals$annotations), FUN = function(a) {
          samp_annot(vals$result_objects[[input$AnnotationMusicaList]], a) <- 
          vals$annotations[[a]]
        })}, error = function(cond) {
          shinyalert::shinyalert(title = "Error", text = cond$message)
          return()
        })
    } else {
      print("Error: selected object does not exist")
    }
  })

  output$CompareResultA <- renderUI({
    tagList(
      selectInput("SelectResultA", h3("Select result object"),
                  choices = c(names(vals$result_objects)))
    )
  })
  
  output$CompareResultB <- renderUI({
    tagList(
      selectInput("SelectResultB", h3("Select comparison result object"),
                  choices = c("cosmic_v2_sigs", "cosmic_v3_dbs_sigs", 
                              "cosmic_v3_indel_sigs", 
                              names(vals$result_objects)))
    )
  })

  observeEvent(input$CompareResults, {
    if (is.null(input$SelectResultA) | input$SelectResultA == "" | 
        input$Threshold == "") {
      browser()
      output$CompareValidate <- renderText({
        validate(
          need(input$SelectResultA != "",
               'Please select a result object to compare.'),
          need(input$Threshold == "",
               'Please provide a similarity threshold from 0 to 1.')
        )
      })
      return()
    }
    if(input$SelectResultB %in% names(cosmic_objects)) {
      other <- cosmic_objects[[input$SelectResultB]]
    } else {
      other <- isolate(vals$result_objects[[input$SelectResultB]])
    }
    tryCatch( {
    output$ComparePlot <- renderPlot({
    #vals$comparisonTable <-
      browser()
      compare_results(isolate(vals$result_objects[[input$SelectResultA]]),
                                            other, threshold = input$Threshold,
                                            result_name = paste(input$CompareResultA, "Signatures"),
                                            other_result_name =
                                              paste(input$CompareResultB, "Signatures"))
  })
    }, error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
    })
  
    # output$CompareTable <- renderTable({
    #   vals$comparisonTable <- compare_results(vals$result_objects[[input$SelectResultA]],
    #                   other, threshold = input$Threshold,
    #                   result_name = paste(input$CompareResultA, "Signatures"),
    #                   other_result_name =
    #                     paste(input$CompareResultB, "Signatures"))
    # })
  })

###############################################################################
  
}
