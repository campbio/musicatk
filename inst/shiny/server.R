library(musicatk)
library(plotly)
library(sortable)

options(shiny.maxRequestSize = 100*1024^2)
source("server_tables.R", local = T)
vals <- reactiveValues(
  musica = NULL,
  files = list(),
  result_objects = list(),
  vals_annotatios = NULL,
  df = NULL,
  musica_contents = NULL,
  musica_upload = NULL)
rv <- reactiveValues(
  data = NULL,
  deletedRows = NULL,
  deletedRowIndices = list()
)
###################### Zainab's Code ##########################################
server <- function(input, output) {
  # Create dynamic table
  vals <- reactiveValues(
    var = NULL,
    genome = NULL,
    musica = NULL,
    musica_upload = NULL,
    result_objects = list(),
    vals_annotatios = NULL,
    data = NULL,
    point_ind = 0,
    annot = NULL
  )
  
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
    showNotification("Import successfully completed!")
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

  observeEvent(input$get_musica_object,{
    vals$genome <- genome()
    vals$musica <- create_musica(x = vals$var, genome = genome(),check_ref_chromosomes = check_chr(),check_ref_bases = check_bases(),
                            convert_dbs = convert_dbs(),standardize_indels = stand_indels())
    showNotification("Musica Object successfully created! ")
    
  })
  output$musica_contents <- renderDataTable({
    req(vals$musica)
    return(head(vals$musica@variants))
    shinyjs::show(id="musica_contents")
    js$enableTabs()
  })
  
  output$musica_contents_summary <- renderText({
    req(vals$musica)
    vt <- unique(vals$musica@variants$Variant_Type) #variant types
    nvt<- table(vals$musica@variants$Variant_Type)
    ns <- length(vals$musica@variants$sample) #sample length
    mylist <- c("No. of Samples:\n",ns,"\n","Variant types",vt,"\n",nvt)
    return(mylist)
    shinyjs::show(id="musica_contents_summary")
    js$enableTabs();
  })
  
  observeEvent(input$reset, {
    removeUI("#musica_contents")
    removeUI("#musica_contents_summary")
    #removeUI("#musica_upload")
    showNotification("Tables cleared!")
  })
  
  observeEvent(input$musica_file,{
    req(input$musica_file)
    vals$musica_upload <- load(input$musica_file$datapath)
    vals$musica_upload <- get(vals$musica_upload)
    vals$musica <- vals$musica_upload
    showNotification("Musica Result Object successfully imported!")
  })
  
  output$musica_upload <- renderDataTable({
    req(vals$musica_upload)
    return(head(vals$musica_upload@musica@variants))
    shinyjs::show(id="musica_upload")
    js$enableTabs();
  })
  output$musica_upload_summary <- renderText({
    req(vals$musica_upload)
    vt <- unique(vals$musica_upload@musica@variants$Variant_Type) #variant types
    nvt<- table(vals$musica_upload@musica@variants$Variant_Type)
    ns <- length(vals$musica_upload@musica@variants$sample) #sample length
    mylist <- c("No. of Samples:\n",ns,"\n","Variant types",vt,"\n",nvt)
    return(mylist)
    shinyjs::show(id="musica_upload_summary")
    js$enableTabs();
  })
  
  observeEvent(input$reset_musica, {
    removeUI("#musica_upload")
    removeUI("#musica_upload_summary")
    showNotification("Tables cleared!")
  })
  
  output$download_musica <- downloadHandler(
    filename = function() {
      paste("musica_variants", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$musica@variants, file, row.names = FALSE)
    }
  )
  
  output$download_musica_result <- downloadHandler(
    filename = function() {
      paste("musica_variants", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$musica_upload@musica@variants, file, row.names = FALSE)
    }
  )
  
  
  
  observeEvent(input$upload, {
    rv$data <- data.frame(input$file[['name']])
    # Clear the previous deletions
    rv$deletedRows <- NULL
    rv$deletedRowIndices = list()
    })
    
    observeEvent(input$deletePressed, {
    rowNum <- parseDeleteEvent(input$deletePressed)
    dataRow <- rv$data[rowNum,]
    
    # Put the deleted row into a data frame so we can undo
    # Last item deleted is in position 1
    rv$deletedRows <- rbind(dataRow, rv$deletedRows)
    rv$deletedRowIndices <- append(rv$deletedRowIndices, rowNum, after = 0)
    
    # Delete the row from the data frame
    rv$data <- rv$data[-rowNum,]
  })
  
  observeEvent(input$undo, {
    if(nrow(rv$deletedRows) > 0) {
      row <- rv$deletedRows[1, ]
      rv$data <- addRowAt(rv$data, row, rv$deletedRowIndices[[1]])
      
      # Remove row
      rv$deletedRows <- rv$deletedRows[-1,]
      # Remove index
      rv$deletedRowIndices <- rv$deletedRowIndices[-1]
    }
  })
  
  # Disable the undo button if we have not deleted anything
  output$undoUI <- renderUI({
    if(!is.null(rv$deletedRows) && nrow(rv$deletedRows) > 0) {
      actionButton('undo', label = 'Undo delete', icon('undo'))
    } else {
      actionButton('undo', label = 'Undo delete', icon('undo'), disabled = TRUE)
    }
  })
  
  output$dtable <- DT::renderDataTable({
    # Add the delete button column
    deleteButtonColumn(rv$data, 'delete_button')
  })


#' Adds a row at a specified index
#'
#' @param df a data frame
#' @param row a row with the same columns as \code{df}
#' @param i the index we want to add row at.
#' @return the data frame with \code{row} added to \code{df} at index \code{i}
addRowAt <- function(df, row, i) {
  # Slow but easy to understand
  if (i > 1) {
    rbind(df[1:(i - 1), ], row, df[-(1:(i - 1)), ])
  } else {
    rbind(row, df)
  }
  
}

#' A column of delete buttons for each row in the data frame for the first column
#'
#' @param df data frame
#' @param id id prefix to add to each actionButton. The buttons will be id'd as id_INDEX.
#' @return A DT::datatable with escaping turned off that has the delete buttons in the first column and \code{df} in the other
deleteButtonColumn <- function(df, id, ...) {
  # function to create one action button as string
  f <- function(i) {
    # https://shiny.rstudio.com/articles/communicating-with-js.html
    as.character(actionButton(paste(id, i, sep="_"), label = NULL, icon = icon('trash'),
                              onclick = 'Shiny.setInputValue(\"deletePressed\",  this.id, {priority: "event"})'))
  }
  
  deleteCol <- unlist(lapply(seq_len(nrow(df)), f))
  
  # Return a data table
  DT::datatable(cbind(delete = deleteCol, df),
                # Need to disable escaping for html as string to work
                escape = FALSE,
                options = list(
                  # Disable sorting for the delete column
                  columnDefs = list(list(targets = 1, sortable = FALSE))
                ))
}

#' Extracts the row id number from the id string
#' @param idstr the id string formated as id_INDEX
#' @return INDEX from the id string id_INDEX
parseDeleteEvent <- function(idstr) {
  res <- as.integer(sub(".*_([0-9]+)", "\\1", idstr))
  if (! is.na(res)) res
}

  
  
  
        
  
    
###############################################################################  
    
    
###################### Nathan's Code ##########################################
  
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
  
  # observeEvent(input$DiscoverSignatures, {
  #   vals$result_objects[[input$MusicaResultName]] <- discover_signatures(
  #     vals$musica, table_name = input$SelectDiscoverTable,
  #     num_signatures = as.numeric(input$NumberOfSignatures),
  #     method = input$Method,
  #     #seed = input$Seed,
  #     nstart = as.numeric(input$nStart))
  # })
  
  #The latest argument for method in discover_signatures() is "algorithm"
  observeEvent(input$DiscoverSignatures, {
    vals$result_objects[[input$MusicaResultName]] <- discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      algorithm = input$Method,
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
      shinyjs::show(id = "CosmicSBSSigs")
      shinyjs::hide(id = "CosmicDBSSigs")
      shinyjs::hide(id = "CosmicINDELSigs")
    } else if (input$CosmicCountTable == "DBS") {
      shinyjs::hide(id = "CosmicSBSSigs")
      shinyjs::show(id = "CosmicDBSSigs")
      shinyjs::hide(id = "CosmicINDELSigs")
    } else {
      shinyjs::hide(id = "CosmicSBSSigs")
      shinyjs::hide(id = "CosmicDBSSigs")
      shinyjs::show(id = "CosmicINDELSigs")
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
  
  # observeEvent(input$Compare, {
  #   tryCatch( {
  #    output$ComparePlot <- renderPlot({comparisons <- 
  #      compare_results(vals$result_objects,
  #                   vals$pred_res,
  #                   threshold = input$Threshold)})
  #    }, error = function(cond) {
  #      shinyalert::shinyalert(title = "Error", text = cond$message)
  #      return()
  #    }
  #   )
  # })
  
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
      shinyjs::show(id = "StrandType")
      shinyjs::hide(id = "GetTableName")
    } else if (input$SelectTable == 5) {
      shinyjs::hide(id = "StrandType")
      shinyjs::show(id = "GetTableName")
    } else {
      shinyjs::hide(id = "GetTableName")
      shinyjs::hide(id = "StrandType")
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
  output$select_res1 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res1",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      )
    )
  })
  
  output$select_res2 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res2",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      )
    )
  })
  
  # observeEvent(input$select_button, {
  #   vals$data <- vals$result_objects[[input$selected_res]]
  # })
  
  observeEvent(input$get_res,{
    vals$data <- res_annot
  })
  
  observeEvent(input$rename,{
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
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
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
    if(input$rename){
      ids <- vector()
      for(i in 1:n){
        ids <- c(ids, input[[paste0("sig", i)]])
      }
      name_signatures(result = vals$result_objects[[input$selected_res1]], ids)
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
          result = vals$result_objects[[input$selected_res1]], 
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
          result = vals$result_objects[[input$selected_res1]], 
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
    if(input$plottype %in% c("box", "violin") & vals$point_ind == 0){
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
      vals$point_ind <- 1
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
      vals$point_ind <- 0
    }
    else{
      return(NULL)
    }
  })
  
  observeEvent(input$group1,{
    if(input$group1 == "annotation" & input$color != "annotation"){
      vals$annot <- as.list(colnames(samp_annot(vals$result_objects[[input$selected_res2]]))[-1])
      names(vals$annot) <- colnames(samp_annot(vals$result_objects[[input$selected_res2]]))[-1]
      insertUI(
        selector = "#insertannot",
        ui = selectInput(
          inputId = "annotation",
          label = "Annotation",
          choices = vals$annot
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
          choices = vals$annot
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
              text = "Available Signatures:",
              labels = as.list(colnames(vals$result_objects[[input$selected_res2]]@signatures)),
              input_id = "sig_from"
            ),
            add_rank_list(
              text = "Selected Signatures:",
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
          result = vals$result_objects[[input$selected_res2]], 
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
          result = vals$result_objects[[input$selected_res2]], 
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

####################Heatmap##############  
  output$select_res_heatmap <- renderUI({
    tagList(
      selectInput(
        inputId = "select_res_heatmap",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      )
    )
  })
  propor <- reactive({
    props <- input$prop
    return(props)
  })
  sel_col_names <- reactive({
    cc <- input$col_names
    return(cc)
  })
  sel_row_names <- reactive({
    rr <- input$row_names
    return(rr)
  })
  zscale <- reactive({
    zscale <- input$scale
    return(zscale)
  })
  
  observeEvent(input$subset, {
    if(input$subset == "signature"){
      insertUI(
        selector = "#sortbysigs",
        ui = tags$div(
          id = "insertsig",
          bucket_list(
            header = "Select signatures to sort",
            group_name = "bucket",
            orientation = "horizontal",
            add_rank_list(
              text = "Available Signatures:",
              labels = as.list(colnames(vals$result_objects[[input$select_res_heatmap]]@signatures)),
              input_id = "sig_from"
            ),
            add_rank_list(
              text = "Selected Signatures:",
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
  observeEvent(input$subset_tum, {
    if(input$subset_tum == "tumors"){
      insertUI(
        selector = "#sortbytum",
        ui = tags$div(
          id = "inserttum",
          selectInput("tum_val","",choices = as.list(unique(vals$result_objects[[input$select_res_heatmap]]@musica@sample_annotations$Tumor_Subtypes)))
          )
        )
      
    }
    else{
      removeUI(selector = "#inserttum")
    }
  })
  
  observeEvent(input$subset_annot, {
    if(input$subset_annot == "annotation"){
      insertUI(
        selector = "#sortbyannot",
        ui = tags$div(
          id = "#insertannot",
          selectInput("annot_val","",choices = as.list(colnames(samp_annot(vals$result_objects[[input$select_res_heatmap]])))
        )
      ))
    }
    else{
      removeUI(selector = "#insertannot")
    }
  })
  observeEvent(input$get_heatmap,{
    output$heatmap <- renderPlot({
    #paste0("Heatmap")
    plot_heatmap(res_annot = vals$result_objects[[input$select_res_heatmap]],proportional = propor(),show_row_names = sel_row_names(),show_column_names = sel_col_names(),scale = zscale(),subset_signatures = c(input$sig_to),subset_tumor = input$tum_val,annotation = input$annot_val)
  })
  }) 
  ########################################

  ##############Clustering################
  output$select_res3 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res3",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      ),
      bsTooltip(id = "selected_res3", title = "Select one musica_result object to visualize signatures.",
                placement = "right", options = list(container = "body"))
    )
  })
  
  observeEvent(input$selected_res3, {
    output$no_cluster1 <- renderUI(
      tagList(
        numericInput(inputId = "numclust1", label = "Max Number of Clusters", 
                     value = 10,
                     min = 2,
                     max = dim(vals$result_objects[[input$selected_res3]]@exposures)[2])
      )
    )
    output$no_cluster2 <- renderUI(
      tagList(
        numericInput(inputId = "numclust2", label = "Max Number of Clusters", 
                     value = 2,
                     min = 2,
                     max = dim(vals$result_objects[[input$selected_res3]]@exposures)[2])
      )
    )
  })
  
  observeEvent(input$explore,{
    method <- input$metric
    clust.method <- input$algorithm1
    n <- input$numclust1
    proportional <- input$proportional2
    output$explore_plot <- renderPlot(
      k_select(
        result = vals$result_objects[[input$selected_res3]],
        method = method,
        clust.method = clust.method,
        n = n,
        proportional = proportional
      )
    )
  })
  
  observeEvent(input$algorithm2, {
    choices <- list(hkmeans = c("Euclidean" = "euclidean", "Manhattan" = "manhattan", "Canberra" = "canberra"),
                    clara = c("Euclidean" = "euclidean", "Manhattan" = "manhattan", "Jaccard" = "jaccard"),
                    kmeans = c("Euclidean" = "euclidean", "Manhattan" = "manhattan", "Jaccard" = "jaccard", 
                               "Cosine" = "cosine", "Canberra" = "canberra"),
                    hclust = c("Euclidean" = "euclidean", "Manhattan" = "manhattan", "Jaccard" = "jaccard", 
                               "Cosine" = "cosine", "Canberra" = "canberra"),
                    pam = c("Euclidean" = "euclidean", "Manhattan" = "manhattan", "Jaccard" = "jaccard", 
                            "Cosine" = "cosine", "Canberra" = "canberra"))
    output$diss <- renderUI(
      selectInput(
        inputId = "diss_method",
        label = "Method for Dissimilarity Matrix",
        choices = choices[[input$algorithm2]]
      )
    )
    if(input$algorithm2 == "hclust"){
      insertUI(
        selector = "#hclust",
        ui = selectInput(
          inputId = "hclust_method",
          label = "Hierarchical Clustering Method",
          choices = c("ward.D" = "ward.D", "ward.D2" = "ward.D2", "single" = "single", "complete" = "complete", 
                      "average" = "average", "mcquitty" = "mcquitty", "median" = "median", 
                      "centroid" = "centroid")
        )
      )
    }
    else{
      removeUI(selector = "div:has(>> #hclust_method)")
    }
    if(input$algorithm2 == "clara"){
      insertUI(
        selector = "#clara",
        ui = numericInput(inputId = "clara_num",
                          label = "No. of Samples for CLARA",
                          value = 5,
                          min = 1,
                          max = dim(vals$result_objects[[input$selected_res3]]@exposures)[2])
      )
    }
    else{
      removeUI(selector = "div:has(>> #clara_num)")
    }
    if(input$algorithm2 %in% c("kmeans", "hkmeans")){
      insertUI(
        selector = "iter",
        ui = numericInput(inputId = "max_iter",
                          label = "Max No. of Iterations",
                          value = 10,
                          min = 1)
      )
    }
    else{
      removeUI(selector = "div:has(>> #max_iter)")
    }
  })
  
  observeEvent(input$group2, {
    if(input$group2 == "annotation"){
      vals$annot <- as.list(colnames(samp_annot(vals$result_objects[[input$selected_res3]]))[-1])
      names(vals$annot) <- colnames(samp_annot(vals$result_objects[[input$selected_res3]]))[-1]
      insertUI(
        selector = "#insertannot2",
        ui = tagList(
          selectInput(
            inputId = "annotation2",
            label = "Annotation",
            choices = vals$annot
          )
        )
      )
    }
    else{
      removeUI(selector = "div:has(>> #annotation2)")
    }
  })
  
  observeEvent(input$cluster_calc, {
    result <- vals$result_objects[[input$selected_res3]]
    nclust <- input$numclust2
    proportional <- input$proportional3
    method <- input$algorithm2
    dis.method <- input$diss_method
    if(!is.null(input$hclust_method)){
      hc.method <- input$hclust_method
    }
    else{
      hc.method <- "ward.D"
    }
    if(!is.numeric(input$clara_num)){
      clara.samples <- 5
    }
    else{
      clara.samples <- input$clara_num
    }
    if(!is.numeric(input$max_iter)){
      iter.max <- 10
    }
    else{
      iter.max <- input$max_iter
    }
    vals$cluster <- cluster_exposure(result = result,
                                     nclust = nclust,
                                     proportional = proportional,
                                     method = method,
                                     dis.method = dis.method,
                                     hc.method = hc.method,
                                     clara.samples = clara.samples,
                                     iter.max = iter.max)
    insertUI(
      selector = "#insert_cluster_table",
      ui = DT::dataTableOutput("cluster_table")
    )
    annot <- samp_annot(vals$result_objects[[input$selected_res3]])
    row.names(annot) <- annot$Samples
    dat <- cbind(annot, vals$cluster)
    output$cluster_table <- DT::renderDataTable(
      DT::datatable(dat[,-1])
    )
  })
  
  observeEvent(input$cluster_vis, {
    if(length(umap(vals$result_objects[[input$selected_res3]])) == 0){
      create_umap(vals$result_objects[[input$selected_res3]])
      result <- vals$result_objects[[input$selected_res3]]
    }
    else{
      result <- vals$result_objects[[input$selected_res3]]
    }
    clusters <- vals$cluster
    group <- input$group2
    if(group == "annotation"){
      annotation <- input$annotation2
    }
    else{
      annotation <- NULL
    }
    output$cluster_plot <- renderPlot(
      plot_cluster(result = result,
                   clusters = clusters,
                   group = group,
                   annotation = annotation)
    )
  })
  ########################################
}
