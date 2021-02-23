library(musicatk)


options(shiny.maxRequestSize = 100*1024^2)
source("server_tables.R", local = T)

###################### Zainab's Code ##########################################
server <- function(input, output) {
  #observeEvent(input$get_musica, {
   # maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  #})
  
  #observeEvent(input$MusicaResults,{
   # data(res_annot)
 # })
  # variants <- eventReactive(input$import,{
  #   req(input$file)
  #   file_name <- input$file$datapath
  #   var <- extract_variants(c(file_name))
  #   removeUI(selector = "div#file_id")
  #   return(var)
  # })
  observeEvent(input$import,{
    req(input$file)
    file_name <- input$file$datapath
    vals$var <- extract_variants(c(file_name))
    removeUI(selector = "div#file_id")
  })
  
  
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
                list( "Common genomes" = list("hg18","hg19","hg38","mm9","mm10"),
                "Genomes" = gg), 
                      width ='100%')
    })
  genome <- reactive({
    gen <- input$GenomeSelect
    gen <- select_genome(gen)
    return(gen)
  })
  output$genome_select <- renderText({
    paste("Genome selected:", input$GenomeSelect)
  })
  
  check_chr <- reactive({ 
    chr <- input$ref_chr
    return(chr)
  })
  check_bases <- reactive({
    bases <- input$ref_bases
    return(bases)
  })
   convert_dbs <- reactive({
    conv_dbs <- input$convert_dbs
    return(conv_dbs)
  })
  stand_indels <- reactive({
    stand_indels <- input$stand_indels
    return(stand_indels)
  })

  })
  
  observeEvent(input$get_musica_object,{
    vals$musica_contents <- create_musica(x = vals$var, genome = genome(),check_ref_chromosomes = check_chr(),check_ref_bases = check_bases(),
                            convert_dbs = convert_dbs(),standardize_indels = stand_indels())
    
  })
  
  
  
  output$musica_contents <- renderDataTable({
    req(vals$musica_contents)
    return(head(vals$musica_contents@variants))
    shinyjs::show(id="musica_contents")
    js$enableTabs();
  })
  
  output$musica_contents_summary <- renderText({
    req(vals$musica_contents)
    vt <- unique(vals$musica_contents@variants$Variant_Type) #variant types
    nvt<- table(vals$musica_contents@variants$Variant_Type)
    ns <- length(vals$musica_contents@variants$sample) #sample length
    mylist <- c("No. of Samples:\n",ns,"\n","Variant types",vt,"\n",nvt)
    return(mylist)
    shinyjs::show(id="musica_contents_summary")
    js$enableTabs();
  })
  
  observeEvent(input$reset, {
    removeUI("#musica_contents")
    removeUI("#musica_contents_summary")
  })
  
  observeEvent(input$musica_file,{
    req(input$musica_file)
    vals$musica_upload <- load(input$musica_file$datapath,.GlobalEnv)
    print(input$musica_file$datapath)
    
  })
  
  output$musica_upload <- renderTable({
    req(vals$musica_upload)
    return(head(vals$musica_upload@musica@variants))
    shinyjs::show(id="musica_upload")
    js$enableTabs();
  },striped = TRUE)
  
  output$download_musica <- downloadHandler(
    filename = function() {
      paste("musica_variants", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$musica_contents@variants, file, row.names = FALSE)
    }
  )
  
    

  
  
  
        
  
    
###############################################################################  
    
    
###################### Nathan's Code ##########################################
  # Create dynamic table
  vals <- reactiveValues(
    musica = NULL,
    result_objects = list(),
    vals_annotatios = NULL
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
                    extract_count_tables(
                      vals$musica)))
    )
  })
  
  observeEvent(input$AddTable, {
    table_name = input$SelectTable
    if (table_name == "SBS192") {
      table_name = input$StrandType
    }
    if(table_name %in% names(extract_count_tables(vals$musica))) {
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
    }
  })
  
  observeEvent(input$DiscoverSignatures, {
    vals$result_objects[[input$MusicaResultName]] <- discover_signatures(
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
    
    vals$result_objects[[input$PredictResultName]] <-
      predict_exposure(vals$musica, g = genome, 
                     table_name = input$SelectPredTable,
                     signature_res = res,
                     algorithm = input$PredictAlgorithm,
                     signatures_to_use = as.numeric(sigs))
  })
  
  observeEvent(input$Compare, {
    tryCatch( {
     output$ComparePlot <- renderPlot({comparisons <- 
       compare_results(vals$result_objects,
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
      show(id = "GetTableName")
    } else {
      hide(id = "GetTableName")
      hide(id = "StrandType")
    }
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
    validate(
      need(input$SelectResultA != "", 
           'Please select a result object to compare')
    )
    if(input$SelectResultB %in% names(cosmic_objects)) {
      other <- cosmic_objects[[input$SelectResultB]]
    } else {
      other <- vals$result_objects[[input$SelectResultB]]
    }
    compare_results(vals$result_objects[[input$SelectResultA]],
                    other,
                    threshold = input$Threshold,
                    metric = input$CompareMetric,
                    result_name = paste(input$CompareResultA, "Signatures"),
                    other_result_name = 
                      paste(input$CompareResultB, "Signatures"))
  })
###############################################################################
  
}
