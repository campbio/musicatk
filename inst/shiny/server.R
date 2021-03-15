library(musicatk)
library(plotly)
library(sortable)

options(shiny.maxRequestSize = 100*1024^2)
source("server_tables.R", local = T)

server <- function(input, output, session) {
  
#################### GENERAL ##################################################  
  vals <- reactiveValues(
    genome = NULL,
    musica = NULL,
    files = NULL,
    result_objects = list(),
    cSigs = NULL, #cosmic sig
    cRes = NULL, #cosmic result object
    annotations = NULL,
    df = NULL,
    musica_contents = NULL,
    musica_upload = NULL,
    genome = genome,
    data = NULL,
    point_ind = 0,
    annot = NULL,
    deletedRows = NULL,
    deletedRowIndices = list()
  )
  

###################### Zainab's Code ##########################################
  
  # Create dynamic table

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
    withProgress(message = "Importing",{
      for (i in 1:15) {
        incProgress(1/15)
        Sys.sleep(0.25)
      }
    })
    
    file_name <- vals$data$datapath
    withCallingHandlers(
      vals$var <- extract_variants(c(file_name)),
      message = function(m)output$text <- renderPrint(m$message))
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

tryCatch({  observeEvent(input$get_musica_object,{
    vals$musica <- create_musica(x = vals$var, genome = genome(),check_ref_chromosomes = check_chr(),check_ref_bases = check_bases(),
                            convert_dbs = convert_dbs(),standardize_indels = stand_indels())
    showNotification("Musica Object successfully created! ")
    
  })},
  error = function(cond){
    shinyalert::shinyalert(title = "Error", text = cond$message)
  },
  warning = function(cond) {
    print(cond$message)
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
    ns <- length(vals$musica@variants$sample) #sample length
    #mylist <- c("No. of Samples:\n",ns,"\n","Variant types",vt)
    mylist <- c("No. of Samples:\n",ns)
    return(mylist)
    shinyjs::show(id="musica_contents_summary")
    js$enableTabs();
  })
  output$musica_contents_table <- renderDataTable({
    req(vals$musica)
    nvt<- as.data.frame(table(vals$musica@variants$Variant_Type))
    return(nvt)
    shinyjs::show(id="musica_contents_table")
    js$enableTabs();
    })
  
  observeEvent(input$reset, {
    removeUI("#musica_contents")
    removeUI("#musica_contents_summary")
    removeUI("#musica_contents_table")
    #removeUI("#musica_upload")
    showNotification("Tables cleared!")
  })
  
  observeEvent(input$musica_button,{
    if(input$musica_button == "result"){
      shinyjs::show(id = "MusicaResultName")
    }
    else if(input$musica_button == "object"){
      shinyjs::hide(id = "MusicaResultName")}
  })
  observeEvent(input$upload_musica,{
    req(input$musica_file)
    if(input$musica_button == "result"){
      vals$musica_upload <- load(input$musica_file$datapath)
      vals$musica_upload <- get(vals$musica_upload)
      vals$result_objects[[input$MusicaResultName]] <- vals$musica_upload
      showNotification("Musica Result Object successfully imported!")
    }
    else if(input$musica_button == "object"){
      vals$musica_upload <- load(input$musica_file$datapath)
      vals$musica_upload <- get(vals$musica_upload)
      vals$musica <- vals$musica_upload 
      showNotification("Musica Object successfully imported!")
    }
    
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
    # Clear the previous deletions
    vals$files <- list(input$file[['name']])
    vals$files <- unlist(vals$files)
    vals$files <- c(vals$files)
    dt <- list(input$file[['datapath']])
    dt <- unlist(dt)
    dt <- c(dt)
    
    #vals$files <- list(a = list(input$file[['name']]))
    #vals$files <- rbindlist(vals$files)
    #vals$files <- setDF(vals$files)
    #vals$files <- split(vals$files,seq(nrow(vals$files)))
    vals$files <- data.frame(files = vals$files,datapath = dt,stringsAsFactors = FALSE)
    vals$data <- vals$files
    vals$deletedRows <- NULL
    vals$deletedRowIndices = list()
    })
    
    observeEvent(input$deletePressed, {
    rowNum <- parseDeleteEvent(input$deletePressed)
    dataRow <- vals$data[rowNum,]
    
    # Put the deleted row into a data frame so we can undo
    # Last item deleted is in position 1
    vals$deletedRows <- rbind(dataRow, vals$deletedRows)
    vals$deletedRowIndices <- append(vals$deletedRowIndices, rowNum, after = 0)
    
    # Delete the row from the data frame
    vals$data <- vals$data[-rowNum,]
  })
  
  observeEvent(input$undo, {
    if(nrow(vals$deletedRows) > 0) {
      row <- vals$deletedRows[1, ]
      vals$data <- addRowAt(vals$data, row, vals$deletedRowIndices[[1]])

      # Remove row
      vals$deletedRows <- vals$deletedRows[-1,]
      # Remove index
      vals$deletedRowIndices <- vals$deletedRowIndices[-1]
    }
  })
  
  #Disable the undo button if we have not deleted anything
  output$undoUI <- renderUI({
    if(!is.null(vals$deletedRows) && nrow(vals$deletedRows) > 0) {
      actionButton('undo', label = 'Undo delete', icon('undo'))
    } else {
      actionButton('undo', label = 'Undo delete', icon('undo'), disabled = TRUE)
    }
  })
  
  output$dtable <- DT::renderDataTable({
    # Add the delete button column
    req(vals$data)
    deleteButtonColumn(vals$data, 'delete_button')
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
  
  deleteCol <- unlist(lapply(seq(nrow(df)), f))
  
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
      selectInput("SelectDiscoverTable", "Select Count Table",
                  choices = names(
                    extract_count_tables(vals$musica))),
      bsTooltip("SelectDiscoverTable",
                "Name of the table to use for signature discovery.", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  output$CombineTable <- renderUI({
    if (length(names(tables(vals$musica))) > 1) {
      tagList(
        box(width = 12,
            helpText("Combine any 2 or more tables contained in your musica object. This is optional."),
        checkboxGroupInput("CombineTables", "Tables to Combine",
                    choices = names(tables(vals$musica))),
        textInput("CombinedTableName", "Name of combined table"),
        uiOutput("CombineWarning"),
        actionButton("Combine", "Build Combined Table"),
        bsTooltip("CombinedTableName",
                  "Combine tables into a single table that can be used for
                  discovery/prediction.", 
                  placement = "bottom", trigger = "hover", options = NULL),
        bsTooltip("Combine",
                  "Combines tables into a single table that can be used for discovery/prediction.", 
                  placement = "bottom", trigger = "hover", options = NULL),
        bsTooltip("CombineTables",
                  "Tables to combine.", 
                  placement = "left", trigger = "hover", options = NULL),
        bsTooltip("CombinedTableName",
                  "Name for the combined table.", 
                  placement = "bottom", trigger = "hover", options = NULL)
        )
      )
    }

  })
  
  observeEvent(input$Combine, {
    if (input$CombinedTableName == "" | length(input$CombineTables) < 2) {
      output$CombineWarning <- renderText({
        validate(
          need(input$CombinedTableName != "",
               'You must provide a name for the new result object.'),
          need(length(input$CombineTables) < 2, 
               'You must select two or more tables to combine.')
        )
      })
      return ()
    }
    shinybusy::show_spinner()
    tryCatch( {
      combine_count_tables(vals$musica, input$CombineTables, 
                         input$CombinedTableName)
    }, error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
      shinybuys::hide_spinner()
      return()
    })
    shinybusy::hide_spinner()
    showNotification("Table created.")
    
  })
  
  output$AllowTable <- renderUI({
    if (!is.null(vals$musica)) {
      tagList(
      actionButton("AddTable", "Create Table"),
      bsTooltip("AddTable", "Create a table containig the mutationl count information of each sample.", placement = "bottom", trigger = "hover",
                options = NULL)
      )

    } else {
      helpText("You must first create or upload a musica object to generate
               count tables.")
    }
  })
  
  observeEvent(input$AddTable, {
    tableName <- input$SelectTable
    if(tableName == "SBS192 - Replication_Strand") {
      tableName <- "SBS192_Rep"
    }
    if(tableName == "SBS192 - Transcript_Strand") {
      tableName <- "SBS192_Trans"
    }
    if(tableName %in% names(tables(vals$musica))) {
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

  #The latest argument for method in discover_signatures() is "algorithm"
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
    shinybusy::show_spinner()
    vals$result_objects[[input$MusicaResultName]] <- discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      algorithm = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$nStart))
    shinybusy::hide_spinner()
  }
  
  observeEvent(input$confirmResultOverwrite, {
    removeModal()
    getResult(input, vals)
    showNotification("Existing result object overwritten.")
  })
  
  output$PredictTable <- renderUI({
    tagList(
      selectInput("SelectPredTable", "Select Counts Table",
                  choices = names(tables(vals$musica))),
      bsTooltip("SelectPredTable",
                "Name of the table used for posterior prediction", 
                placement = "bottom", trigger = "hover", options = NULL)
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
     shinyjs:: show(id = "CosmicSBSSigs")
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
      vals$cSigs <- "CosmicSBSSigs"
      vals$cRes <- cosmic_v2_sigs
    } else if (input$CosmicCountTable == "DBS") {
      vals$cSigs <- "CosmicDBSSigs"
      vals$cRes <- cosmic_v3_dbs_sigs
      } else {
      vals$cSigs <- "CosmicINDELSigs"
      vals$cRes <- cosmic_v3_indel_sigs
    }
    if (input$PredictResultName == "" | length(input[[vals$cSigs]]) < 2) {
      output$PredictWarning <- renderText({
        validate(
          need(input$PredictResultName != "",
               'You must provide a name for the new result object.'),
          need(length(c(input[[vals$cSigs]])) >= 2,
               'You must select two or more signatures to predict.')
        )
      })
      return ()
    }
    if(input$PredictResultName %in% names(vals$result_objects)) {
      showModal(modalDialog(
        title = "Existing Result Object.",
        "Do you want to overwrite the existing result object?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirmPredictOverwrite", "OK"),
          modalButton("Cancel"))
      ))
    } else {
      getPredict(input, vals)
      showNotification(paste0("Musica result object, ", input$PredictResultName,
                              ", was created"))
    }

  })

  getPredict <- function(inputs, vals) {
    shinybusy::show_spinner()
    vals$result_objects[[input$PredictResultName]] <-
      predict_exposure(vals$musica, g = vals$genome, 
                       table_name = input$SelectPredTable,
                       signature_res = vals$cRes,
                       algorithm = input$PredictAlgorithm,
                       signatures_to_use = c(as.numeric(input[[vals$cSigs]])))
    shinybusy::hide_spinner()
  }
  
  observeEvent(input$confirmPredictOverwrite, {
    removeModal()
    getPredict(inputs, vals)
    showNotification("Existing result overwritten.")
  })
  
  observeEvent(input$confirmOverwrite, {
    removeModal()
    add_tables(input, vals)
    showNotification("Existing table overwritten.")
  })

  output$annotations <- renderDataTable({
    file <- input$AnnotationsFile
    ext <- tools::file_ext(file$datapath)
    req(file)
    vals$annotations <- read.delim(file$datapath, 
                                   header = input$AnnotationHeader,
                                   sep = input$AnnotationDelimiter)
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
      selectInput("SelectResultA", "Select result object",
                  choices = c(names(vals$result_objects))),
      bsTooltip("SelectResultA",
                "A musica result object", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  output$CompareResultB <- renderUI({
    tagList(
      selectInput("SelectResultB", "Select comparison result object",
                  choices = c("cosmic_v2_sigs", "cosmic_v3_dbs_sigs", 
                              "cosmic_v3_indel_sigs", 
                              names(vals$result_objects))),
      bsTooltip("SelectResultB",
                "A second musica result object", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })

  observeEvent(input$CompareResults, {
    if (is.null(input$SelectResultA) | input$SelectResultA == "" | 
        input$Threshold == "") {
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
    tryCatch({
      print("comparing")
      output$CompareTable <- renderDataTable({
        shinybusy::show_spinner()
        
      isolate(vals$comparison <-
        compare_results(isolate(vals$result_objects[[input$SelectResultA]]),
                                            other, threshold = as.numeric(input$Threshold)))#,
                                            # result_name = paste(input$CompareResultA, "Signatures"),
                                            # other_result_name =
                                            #   paste(input$CompareResultB, "Signatures"))
        shinybusy::hide_spinner()
        return(isolate(vals$comparison))
      })
      output$DownloadComparison <- renderUI({
        tagList(
        downloadButton("DownloadCompare", "Download"),
        bsTooltip("DownloadCompare",
                  "Download the comparison table", 
                  placement = "bottom", trigger = "hover", options = NULL)
        )
      })
    }, error = function(cond) {
    shinybusy::hide_spinner()
    shinyalert::shinyalert(title = "Error", text = cond$message)
    })

  })
  
  output$DownloadCompare <- downloadHandler(
    filename = function() { paste0("Sig-Compare", Sys.Date(), ".csv")
      },
    content = function(file) {
      write.csv(vals$comparison, file)
    }
  )


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
}
