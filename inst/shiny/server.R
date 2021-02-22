library(musicatk)
library(plotly)
library(sortable)

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
 
##################Visualization#################   
  v <- reactiveValues(
    data = NULL,
    entry_id = NULL,
    j = 0,
    del_ind = 0,
    sig_list = NULL,
    point_ind = 0,
    annot = NULL
  )
  
  observeEvent(input$get_res,{
    v$data <- res_annot
  })
  
  observeEvent(input$rename,{
    n <- ncol(v$data@signatures)
    for(i in 1:n){
      id <- paste0("sig", i)
      if(input$rename){
        insertUI(
          selector = "#signame",
          ui = textInput(inputId = id, paste0("Signature", i))
        )
      }
      else{
        removeUI(selector = paste0("div:has(> #", id, ")"))
      }
    }
  }, ignoreInit = TRUE)
  
  get_sig_option <- function(input){
    n <- ncol(v$data@signatures)
    if(input$rename){
      ids <- vector()
      for(i in 1:n){
        ids <- c(ids, input[[paste0("sig", i)]])
      }
      name_signatures(result = v$data, ids)
    }
    legend <- input$legend1
    text_size <- input$textsize1
    facet_size <- input$facetsize
    show_x_labels <- input$xlab1
    same_scale <- input$scale1
    plotly = input$plotly1
    options <- list(legend, text_size, facet_size, show_x_labels, same_scale, plotly)
    return(options)
  }
  
  observeEvent(input$get_plot1,{
    options <- get_sig_option(input)
    if(options[[6]]){
      removeUI(selector = "#sigplot_plot")
      insertUI(
        selector = "#plotdiv1",
        ui = plotlyOutput(outputId = "sigplot_plotly")
      )
      output$sigplot_plotly <- renderPlotly(
        plot_signatures(
          result = v$data, 
          legend = options[[1]],
          plotly = options[[6]],
          text_size = options[[2]],
          facet_size = options[[3]],
          show_x_labels = options[[4]],
          same_scale = options[[5]]
        )
      )
    }
    else{
      removeUI(selector = "#sigplot_plotly")
      insertUI(
        selector = "#plotdiv1",
        ui = plotOutput(outputId = "sigplot_plot")
      )
      output$sigplot_plot <- renderPlot(
        plot_signatures(
          result = v$data, 
          legend = options[[1]],
          plotly = options[[6]],
          text_size = options[[2]],
          facet_size = options[[3]],
          show_x_labels = options[[4]],
          same_scale = options[[5]]
        )
      )
    }
  })
  
  observeEvent(input$plottype, {
    if(input$plottype %in% c("box", "violin") & v$point_ind == 0){
      insertUI(
        selector = "#points",
        ui = tags$div(
          id = "point_opt",
          checkboxInput(inputId = "addpoint", label = "Add Points", value = TRUE),
          numericInput(inputId = "pointsize", label = "Point Size", value = 2) 
        )
      )
      removeUI(selector = "#group1")
      insertUI(
        selector = "#insert_group",
        ui = radioButtons(
          inputId = "group1",
          label = "Group By",
          choices = list("Signature" = "signature",
                         "Annotation" = "annotation"),
          inline = TRUE,
          selected = "signature"
        )
      )
      v$point_ind <- 1
    }
    else if(input$plottype == "bar"){
      removeUI(selector = "#point_opt")
      removeUI(selector = "#group1")
      insertUI(
        selector = "#insert_group",
        ui = radioButtons(
          inputId = "group1",
          label = "Group By",
          choices = list("None" = "none","Signature" = "signature",
                         "Annotation" = "annotation"),
          inline = TRUE,
          selected = "none"
        )
      )
      v$point_ind <- 0
    }
    else{
      return(NULL)
    }
  })
  
  observeEvent(input$group1,{
    if(input$group1 == "annotation" & input$color != "annotation"){
      v$annot <- as.list(colnames(samp_annot(v$data))[-1])
      names(v$annot) <- colnames(samp_annot(v$data))[-1]
      insertUI(
        selector = "#insertannot",
        ui = selectInput(
          inputId = "annotation",
          label = "Annotation",
          choices = v$annot
        )
      )
    }
    else{
      if(input$color != "annotation"){
        removeUI(selector = "div:has(>> #annotation)")
      }
    }
  })
  
  observeEvent(input$color,{
    if(input$color == "annotation" & input$group1 != "annotation"){
      insertUI(
        selector = "#insertannot",
        ui = selectInput(
          inputId = "annotation",
          label = "Annotation",
          choices = v$annot
        )
      )
    }
    else{
      if(input$group1 != "annotation"){
        removeUI(selector = "div:has(>> #annotation)")
      }
    }
  })
  
  observeEvent(input$sort, {
    if(input$sort == "signature"){
      insertUI(
        selector = "#sortbysig",
        ui = tags$div(
          id = "insertsig",
          bucket_list(
            header = "Select signatures to sort",
            group_name = "bucket",
            orientation = "horizontal",
            add_rank_list(
              text = "From:",
              labels = as.list(colnames(v$data@signatures)),
              input_id = "sig_from"
            ),
            add_rank_list(
              text = "To:",
              labels = NULL,
              input_id = "sig_to"
            )
          )
        )
      )
    }
    else{
      removeUI(selector = "#insertsig")
    }
  })
  
  get_exp_option <- function(input){
    plot_type <- input$plottype
    proportional <- input$proportional
    group_by <- input$group1
    color_by <- input$color
    if(input$group1 == "annotation" | input$color == "annotation"){
      annot <- input$annotation
    }
    else{
      annot <- NULL
    }
    if(!is.numeric(input$numsamp)){
      num_samples <- NULL
    }
    else{
      num_samples <- input$numsamp
    }
    sort_by <- input$sort
    if(sort_by == "signature"){
      sort_samples <- input$sig_to
    }
    else{
      sort_samples <- sort_by
    }
    if(!is.numeric(input$theta)){
      threshold <- NULL
    }
    else{
      threshold <- input$theta
    }
    same_scale <- input$scale2
    label_x_axis <- input$xlab2
    legend <- input$legend2
    if(length(input$addpoint) == 0){
      add_points <- FALSE
    }
    else{
      add_points <- input$addpoint
    }
    point_size <- input$pointsize
    plotly <- input$plotly2
    options <- list(plot_type, proportional, group_by, color_by,
                    annot, num_samples, sort_samples, threshold,
                    same_scale, label_x_axis, legend, add_points,
                    point_size, plotly)
    return(options)
  }
  
  observeEvent(input$get_plot2,{
    options <- get_exp_option(input)
    if(options[[14]]){
      removeUI(selector = "#expplot_plot")
      insertUI(
        selector = "#plotdiv2",
        ui = plotlyOutput(outputId = "expplot_plotly")
      )
      output$expplot_plotly <- renderPlotly(
        plot_exposures(
          result = v$data, 
          plot_type = options[[1]], 
          proportional = options[[2]],
          group_by = options[[3]],
          color_by = options[[4]],
          annotation = options[[5]],
          num_samples = options[[6]],
          sort_samples = options[[7]],
          threshold = options[[8]],
          same_scale = options[[9]],
          label_x_axis = options[[10]],
          legend = options[[11]],
          add_points = options[[12]],
          point_size = options[[13]],
          plotly = options[[14]]
        )
      )
    }
    else{
      removeUI(selector = "#expplot_plotly")
      insertUI(
        selector = "#plotdiv2",
        ui = plotOutput(outputId = "expplot_plot")
      )
      output$expplot_plot <- renderPlot(
        plot_exposures(
          result = v$data, 
          plot_type = options[[1]], 
          proportional = options[[2]],
          group_by = options[[3]],
          color_by = options[[4]],
          annotation = options[[5]],
          num_samples = options[[6]],
          sort_samples = options[[7]],
          threshold = options[[8]],
          same_scale = options[[9]],
          label_x_axis = options[[10]],
          legend = options[[11]],
          add_points = options[[12]],
          point_size = options[[13]],
          plotly = options[[14]]
        )
      )
    }
  })

################################################
}
