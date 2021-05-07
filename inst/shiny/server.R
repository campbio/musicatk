library(musicatk)
library(plotly)
library(sortable)
library(shinyBS)
library(shinyalert)
library(TCGAbiolinks)

options(shiny.maxRequestSize = 1000*1024^2)
source("server_tables.R", local = T)

server <- function(input, output, session) {
#################### GENERAL ##################################################  
  addCssClass(selector = "a[data-value='musica']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='annotations']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='tables']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='discover']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='predict']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='visualization']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='compare']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='differentialanalysis']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='cluster']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='heatmap']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='download']", class = "inactiveLink")
  
  vals <- reactiveValues(
    genome = NULL,
    musica = NULL,
    files = NULL,
    result_objects = list(),
    pSigs = NULL, #cosmic sig
    pRes = NULL, #cosmic result object
    annotations = NULL,
    diff = NULL,
    df = NULL,
    musica_upload = NULL,
    data = NULL,
    point_ind = 0,
    annot = NULL,
    deletedRows = NULL,
    deletedRowIndices = list(),
    cluster = NULL,
    var = NULL,
    musica_name_user = NULL,
    musica_message = NULL,
    sort_sigs = NULL
  )
  
  observeEvent(input$menu, {
    if(input$menu == "musica"){
      if(is.null(vals$var)){
        shinyalert::shinyalert("Error", "No data was uploaded. Please go to \"Import\" and upload your data.", "error")
        updateTabItems(session, "menu", "import")
      }
      else{
        removeCssClass(selector = "a[data-value='musica']", class = "inactiveLink")
      }
    }
    else if(input$menu == "annotations") {
      if(is.null(vals$musica) && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error", "No musica object was found. Please go to \"Import\" or \"Create Musica Object\" to upload or create an object.", "error")
        updateTabItems(session, "menu", "import")
      }
    }
    else if(input$menu == "download") {
      if(is.null(vals$musica) && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error", "No musica object was found. Please go to \"Import\" or \"Create Musica Object\" to upload or create an object.", "error")
        updateTabItems(session, "menu", "import")
      }
    }
    else if(input$menu == "tables") {
      if(is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload)){
        shinyalert::shinyalert("Error", "No data was uploaded. Please go to \"Import\" and upload your data.", "error")
        updateTabItems(session, "menu", "import")
      }
      else if(!is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload)){
        shinyalert::shinyalert("Error", "No musica object was created. Please go to \"Create Musica Object\" 
                               to create a musica object or go to \"Import\" -> \"Import Musica Result Object\" to upload a
                               muscia object.", "error")
        updateTabItems(session, "menu", "musica")
      }
      else{
        removeCssClass(selector = "a[data-value='tables']", class = "inactiveLink")
      }
    }
    else if(input$menu %in% c("discover", "predict")){
      if(is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload)){
        shinyalert::shinyalert("Error", "No data was uploaded. Please go to \"Import\" and upload your data.", "error")
        updateTabItems(session, "menu", "import")
      }
      else if(!is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload)){
        shinyalert::shinyalert("Error", "No musica object was created. Please go to \"Create Musica Object\" 
                               to create a musica object or go to \"Import\" -> \"Import Musica Result Object\" to upload a
                               muscia object.", "error")
        updateTabItems(session, "menu", "musica")
      }
      else if(!is.null(vals$var) && !is.null(vals$musica) && length(tables(vals$musica)) == 0 && is.null(vals$musica_upload)){
        shinyalert::shinyalert("Error", "No mutation count table was created. Please go to \"Build Tables\" to create count table.",
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if(!is.null(vals$var) && is.null(vals$musica) && !is.null(vals$musica_upload) && length(tables(vals$musica_upload)) == 0){
        shinyalert::shinyalert("Error", "No mutation count table was created. Please go to \"Build Tables\" to create count table.",
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else{
        removeCssClass(selector = "a[data-value='discover']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='predict']", class = "inactiveLink")
      }
    }
    else if(input$menu %in% c("visualization", "compare", "differentialanalysis", "cluster", "heatmap")){
      if(is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload) && length(vals$result_objects) == 0){
        shinyalert::shinyalert("Error", "No data was uploaded. Please go to \"Import\" and upload your data.", "error")
        updateTabItems(session, "menu", "import")
      }
      else if(!is.null(vals$var) && is.null(vals$musica) && is.null(vals$musica_upload) && length(vals$result_objects) == 0){
        shinyalert::shinyalert("Error", "No musica object was created. Please go to \"Create Musica Object\" 
                               to create a musica object or go to \"Import\" -> \"Import Musica Result Object\" to upload a
                               muscia object.", "error")
        updateTabItems(session, "menu", "musica")
      }
      else if(!is.null(vals$var) && !is.null(vals$musica) && length(tables(vals$musica)) == 0 && is.null(vals$musica_upload) && length(vals$result_objects) == 0){
        shinyalert::shinyalert("Error", "No mutation count table was created. Please go to \"Build Tables\" to create count table.",
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if(!is.null(vals$var) && is.null(vals$musica) && !is.null(vals$musica_upload) && length(tables(vals$musica_upload)) == 0 && length(vals$result_objects) == 0){
        shinyalert::shinyalert("Error", "No mutation count table was created. Please go to \"Build Tables\" to create count table.",
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if(!is.null(vals$var) && (!is.null(vals$musica) || !is.null(vals$musica_upload)) && (length(tables(vals$musica)) != 0 || length(tables(vals$musica_upload)) != 0) && length(vals$result_objects) == 0){
        shinyalert::shinyalert("Error", "No musica_result object was created. Please go to \"Signatures and Exposures\" -> 
                               \"Discover Signatures and Exposures\" to create musica_result object.", "error")
        updateTabItems(session, "menu", "discover")
      }
      else{
        removeCssClass(selector = "a[data-value='visualization']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='compare']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='differentialanalysis']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='cluster']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='heatmap']", class = "inactiveLink")
      }
    }
    else{
      return()
    }
  })
  
  observeEvent(vals$var, {
    if(!is.null(vals$var)){
      removeCssClass(selector = "a[data-value='musica']", class = "inactiveLink")
    }
  })
  
  observeEvent(vals$musica, {
    if(!is.null(vals$musica)){
      removeCssClass(selector = "a[data-value='tables']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='annotations']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='download']", class = "inactiveLink")
    }
    if(length(tables(vals$musica)) != 0){
      removeCssClass(selector = "a[data-value='discover']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='predict']", class = "inactiveLink")
    }
  })
  
  observeEvent(vals$musica_upload, {
    if(!is.null(vals$musica_upload)){
      removeCssClass(selector = "a[data-value='tables']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='annotations']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='download']", class = "inactiveLink")
    }
    if(length(tables(vals$musica_upload)) != 0){
      removeCssClass(selector = "a[data-value='discover']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='predict']", class = "inactiveLink")
    }
  })
  
  observeEvent(vals$result_objects, {
    if(length(vals$result_objects) > 0){
      removeCssClass(selector = "a[data-value='annotations']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='visualization']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='compare']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='differentialanalysis']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='cluster']", class = "inactiveLink")
      removeCssClass(selector = "a[data-value='heatmap']", class = "inactiveLink")
    }
  })
  
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
  output$tcga_tumor <- renderUI({
    hr()
    p <- TCGAbiolinks:::getGDCprojects()$project_id
    t <- grep("TCGA",p)
    p <- p[t]
    p <- gsub(".*-","",p)
    names(p) <- c("UCEC - Uterine Corpus Endometrial Carcinoma","LGG - Brain Lower Grade Glioma",
                  "SARC - Sarcoma","PAAD - Pancreatic adenocarcinoma","ESCA - Esophageal carcinoma",
                  "PRAD - Prostate adenocarcinoma","LAML - Acute Myeloid Leukemia",
                  "KIRC - Kidney renal clear cell carcinoma","PCPG - Pheochromocytoma and Paraganglioma",
                  "HNSC - Head and Neck squamous cell carcinoma","OV - Ovarian serous cystadenocarcinoma",
                  "GBM - Glioblastoma multiforme","UCS - Uterine Carcinosarcoma","MESO - Mesothelioma","TGCT - Testicular Germ Cell Tumors",
                  "KICH - Kidney Chromophobe","READ - Rectum adenocarcinoma","UVM - Uveal Melanoma",
                  "THCA - Thyroid carcinoma","LIHC - Liver hepatocellular carcinoma","THYM - Thymoma",
                  "CHOL - Cholangiocarcinoma","DLBC - Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                  "KIRP - Kidney renal papillary cell carcinoma","BLCA - Bladder Urothelial Carcinoma","BRCA - Breast invasive carcinoma",
                  "COAD - Colon adenocarcinoma","CESC - Cervical squamous cell carcinoma and endocervical adenocarcinoma",
                  "LUSC - Lung squamous cell carcinoma","STAD - Stomach adenocarcinoma",
                  "SKCM - Skin Cutaneous Melanoma","LUAD - Lung adenocarcinoma","ACC - Adrenocortical carcinoma")
    textInput("tcga_tumor","Enter TCGA tumor type")
    tags$style("#tcga_tumor {
                    font-size:8px;
                    height:10px;
           }")
    checkboxGroupInput("tcga_tumor"," ",choices = as.list(p))
    #selectInput("tcga_tumor","Select Tumor",choices = as.list(p),multiple = TRUE)
  })
  
  observeEvent(input$import_tcga,{
    req(input$import_tcga)
    # validate(
    #   need(input$import_tcga, 'Check at least one letter!')
    # )
    if(!is.null(input$tcga_tumor)){
      shinybusy::show_spinner()
      maf <- TCGAbiolinks::GDCquery_Maf(input$tcga_tumor,pipelines = "mutect")
      vals$var <- extract_variants_from_maf_file(maf)
      showNotification("TCGA dataset successfully imported!")
      shinybusy::hide_spinner()
    }
    else if(is.null(input$tcga_tumor)){
      shinyalert("Error: No tumor found. Please select a tumor from the list!")
      }
  })
  
  output$tcga_contents <- renderDataTable({
    req(vals$var)
    req(input$tcga_tumor)
    return(head(vals$var))
    shinyjs::show(id="tcga_contents")
    js$enableTabs()
  })
  
  observeEvent(input$import,{
    req(input$file)

    # withProgress(message = "Importing",{
    #   for (i in 1:15) {
    #     incProgress(1/15)
    #     Sys.sleep(0.25)
    #   }
    # })
    
    file_name <- vals$data$datapath
    if(all(tools::file_ext(file_name) != c("maf","vcf"))){
      shinyalert::shinyalert("Error: File format not supported! Please upload .maf or .vcf files")
    }
    else{
     #withCallingHandlers(
        shinybusy::show_spinner()
        vals$var <- extract_variants(c(file_name))
        shinybusy::hide_spinner()
        showNotification("Import successfully completed!")
        #message = function(m)output$text <- renderPrint(m$message))
        #removeUI(selector = "div#file_id")
        
        }
  })

  output$genome_list <- renderUI({
    g <- BSgenome::available.genomes()
    g <-strsplit(g,",")
    gg <- gsub("^.*?\\.","", g)
    selectInput("GenomeSelect", "Choose genome:",
                list( "Common genomes" = list("hg38","hg19","hg18","mm9","mm10"),
                      "Genomes" = gg), 
                width ='100%')
  })
  
  # genome <- reactive({
  #   gen <- input$GenomeSelect
  #   gen <- select_genome(gen)
  #   return(gen)
  # })
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

tryCatch({observeEvent(input$get_musica_object,{
    shinybusy::show_spinner()
    vals$genome <- input$GenomeSelect
    #vals$genome <- select_genome(vals$genome)
    if(!is.null(vals$var)){
      vals$musica <- create_musica(x = vals$var, genome = select_genome(vals$genome),check_ref_chromosomes = check_chr(),check_ref_bases = check_bases(),
                            convert_dbs = convert_dbs(),standardize_indels = stand_indels())
     #vals$musica_message <- capture.output(data <- )
      
      if(req(input$get_musica_object)){
        shinyjs::show(id = "download_musica")}
      else {
        shinyjs::hide(id = "download_musica")
      }
      shinybusy::hide_spinner()
      showNotification("Musica Object successfully created! ")
    }
    else{
      shinyalert("Error: Please import files first in the Import section!")
    }
  })},
  error = function(cond){
    shinyalert::shinyalert(title = "Error", text = cond$message)
  },
  warning = function(cond) {
    shinyalert::shinyalert(title = "Error", text = cond$message)
  })
 
 output$musica_console <- renderPrint({
   return(print(vals$musica_message))
 })
 
  output$musica_contents <- renderDataTable({
    req(vals$var)
    return(head(vals$var))
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
  
  
  observe(
    if(!is.null(req(input$musica_file))){
    vals$musica_name_user <- tools::file_path_sans_ext(input$musica_file$name)
  })
  
  output$MusicaResultName <- renderUI(
    textInput("MusicaResultName", value = paste0(vals$musica_name_user), h3("Name your musica result object:"))
  )
  observeEvent(input$musica_button,{
    if(input$musica_button == "result"){
      shinyjs::show(id = "MusicaResultName")
    }
    else if(input$musica_button == "object"){
      shinyjs::hide(id = "MusicaResultName")}
  })
  observeEvent(input$upload_musica,{
    req(input$musica_file)
    if(all(tools::file_ext(input$musica_file$name) != c("rda","rds"))){
      shinyalert::shinyalert("Error: File format not supported! Please upload .rda or .rds files")
    }
    else{
      if(input$musica_button == "result"){
        if(all(tools::file_ext(input$musica_file$name) == "rda")){
          vals$musica_upload <- load(input$musica_file$datapath)
          vals$musica_upload <- get(vals$musica_upload)
          vals$result_objects[[input$MusicaResultName]] <- vals$musica_upload
        }
        else if(all(tools::file_ext(input$musica_file$name) == "rds")){
          vals$musica_upload <- readRDS(input$musica_file$datapath)
          vals$result_objects[[input$MusicaResultName]] <- vals$musica_upload
        }
        showNotification("Musica Result Object successfully imported!")
      }
      else if(input$musica_button == "object"){
        if(all(tools::file_ext(input$musica_file$name) == "rda")){
          vals$musica_upload <- load(input$musica_file$datapath)
          vals$musica_upload <- get(vals$musica_upload)
          vals$musica <- vals$musica_upload 
        }
        else if(all(tools::file_ext(input$musica_file$name) == "rds")){
          vals$musica_upload <- readRDS(input$musica_file$datapath)
          vals$musica <- vals$musica_upload 
        }
        showNotification("Musica Object successfully imported!")
      }}
    
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
  
  output$download_musica_object <- downloadHandler(
    filename = function() {
      paste("musica_object", ".rda", sep = "")
    },
    content = function(file) {
      save(as.data.frame(vals$musica), file = filename)
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
  data("cosmic_v3_sbs_sigs")
  data("cosmic_v3_dbs_sigs")
  data("cosmic_v3_indel_sigs")
  cosmic_objects <- list("cosmic_v3_sbs_sigs" = cosmic_v3_sbs_sigs,
                        # "cosmic_v3_sbs_sigs_exome" = cosmic_v3_sbs_sigs_exome,
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
    if (length(names(extract_count_tables(vals$musica))) > 1) {
      tagList(
        box(width = 12,
            helpText("Combine any 2 or more tables contained in your musica object. This is optional."),
        checkboxGroupInput("CombineTables", "Tables to Combine",
                    choices = names(extract_count_tables(vals$musica))),
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
      tagList(
      helpText("You must first create or upload a musica object to generate
               count tables.")
      )
    }
  })
  
  observeEvent(input$AddTable, {
    # if
    # output$TableGenomeWarning <- renderText({
    #   validate(
    #     need(!is.null(vals$genome),
    #          'Please select a reference genome in .')
    #   )
    # })
    tableName <- input$SelectTable
    if(tableName == "SBS192 - Replication_Strand") {
      tableName <- "SBS192_Rep"
    }
    if(tableName == "SBS192 - Transcript_Strand") {
      tableName <- "SBS192_Trans"
    }
    if(tableName %in% names(extract_count_tables(vals$musica))) {
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

  #The latest argument for method in discover_signatures() is "algorithm"
  observeEvent(input$DiscoverSignatures, {
    if (input$DiscoverResultName == "" |
        input$NumberOfSignatures == "" |
        input$nStart == "" |
        extract_count_tables(vals$musica)[["SBS96"]]@count_table < 2 |
        input$NumberOfSignatures < 2) {
      output$DiscoverWarning <- renderText({
        validate(
          need(input$DiscoverResultName != "",
               'You must provide a name for the new result object.'),
          need(input$NumberOfSignatures != "",
               'You must specify the number of expected signatures.'),
          need(input$nStart != "",
               "Please specify the number of random starts."),
          need(input$NumberOfSignatures >= 2,
               "Must specify 2 or more signatures."),
          need(extract_count_tables(vals$musica)[["SBS96"]]@count_table >= 2,
            "You must provide 2 or more samples")
        )
      })
      return ()
    }
    if(input$DiscoverResultName %in% names(vals$result_objects)) {
      showModal(modalDialog(
        title = "Existing Result Object.",
        "Do you want to overwrite the existing result object?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirmResultOverwrite", "OK"),
          modalButton("Cancel"))
      ))
    } else {
      discSigs(input, vals)
      showNotification(paste0("Musica Result object, ", input$DiscoverResultName,
                            ", was created"))
    }
  })
  
  discSigs <- function(input, vals) {
    #shinybusy::show_spinner()
    setResult(input$DiscoverResultName, discover_signatures(
      vals$musica, table_name = input$SelectDiscoverTable,
      num_signatures = as.numeric(input$NumberOfSignatures),
      algorithm = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$nStart)))
    #shinybusy::hide_spinner()
  }
  
  observeEvent(input$confirmResultOverwrite, {
    removeModal()
    discSigs(input, vals)
    showNotification("Existing result object overwritten.")
  })
  
  output$PredictTable <- renderUI({
    tagList(
      selectInput("SelectPredTable", "Select Counts Table",
                  choices = names(extract_count_tables(vals$musica))),
      bsTooltip("SelectPredTable",
                "Name of the table used for posterior prediction", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  output$AnnotationMusicaList <- renderUI({
    if(is.null(vals$musica)) {
      tagList(
        selectInput("AnnotationMusicaList", "Select object",
                    choices = list("result objects" = list(names(vals$result_objects))))        
      )
    } else {
      tagList(
        selectInput("AnnotationMusicaList", "Select object",
                    choices = list("musica object" = list("musica"),
                    "result objects" = list(names(vals$result_objects))))
      )
    }
  })

  # Discover Musica Result Object
  output$DiscoverResultName <- renderUI({
    name <- input$SelectDiscoverTable
    tagList(
      textInput("DiscoverResultName", "Name for musica result object", 
                value = paste0(name, "-Result")),
      bsTooltip("DiscoverResultName",
                "Name for the newly created musica result object.", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  # observeEvent(input$SelectDiscoverTable, {
  #   updateTextInput(session, "MusicaResultName",
  #   value = paste0(input$SelectDiscoverTable, "-Result"))
  # })
  # 
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
  
  output$PredictResultName <- renderUI({
    name <- names(extract_count_tables(vals$musica))[1]
    tagList(
      textInput("PredictResultName", "Name for musica result object",
        value = paste0(name, "-Predict-Result")),
      bsTooltip("PredictResultName",
                "Name for the newly created musica result object.", 
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  output$PredictedResult <- renderUI({
    other <- list(names(vals$result_objects))
    if(length(names(vals$result_objects)) > 1) {
      other <- names(vals$result_objects)
    }
    tagList(
      selectInput("PredictedResult", "Result to Predict",
                  choices = list("Cosmic" = list(
                    "Cosmic V3 SBS Signatures" = "cosmic_v3_sbs_sigs",
                    "Cosmic V3 DBS Signatures" = "cosmic_v3_dbs_sigs",
                    "Cosmic V3 INDEL Signatures" = "cosmic_v3_indel_sigs"),
                                 "Other" = other),
                  selected = "cosmic_v3_sbs_sigs"),
      bsTooltip("PredictedResult",
                "Result object containing the signatures to predict",
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  output$PrecitedSignatures <- renderUI({
    if(input$PredictedResult %in% names(cosmic_objects)) {
        vals$pSigs <- colnames(signatures(
          cosmic_objects[[input$PredictedResult]]))
        vals$pRes <- cosmic_objects[[input$PredictedResult]]
    }
    else {
        vals$pSigs <- colnames(signatures(
          vals$result_objects[[input$PredictedResult]]))
        vals$pRes <- vals$result_objects[[input$PredictedResult]]
    }
    tagList(
      checkboxGroupInput("PredSigs", "Sigantures to Predict",
                         choices = vals$pSigs, inline = T,
                         selected = vals$pSigs),
      bsTooltip("PredSigs",
                "Signatures to predict.",
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })
  
  observeEvent(input$PredictSigs, {
    if (input$PredictResultName == "" | length(input$PredSigs) < 2) {
      output$PredictWarning <- renderText({
        validate(
          need(input$PredictResultName != "",
               'You must provide a name for the new result object.'),
          need(length(c(input$PredSigs)) >= 2,
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
  
  observeEvent(input$PredictAlgorithm, {
    if(input$PredictAlgorithm == "deconstructSigs"){
      shinyjs::show(id = "PredictGenomeList")
    } else {
      shinyjs::hide(id = "PredictGenomeList")
    }
  })

  
  setResult <- function(x, y) {
    vals$result_objects[[x]] <- y
  }

  getPredict <- function(inputs, vals) {
    setResult(input$PredictResultName,
      predict_exposure(vals$musica, g = select_genome(input$PredictGenomeList),
                       table_name = input$SelectPredTable,
                       signature_res = vals$pRes,
                       algorithm = input$PredictAlgorithm,
                       signatures_to_use = input$PredSigs))
  }
  
  observeEvent(input$confirmPredictOverwrite, {
    removeModal()
    getPredict(inputs, vals)
    showNotification("Existing result overwritten.")
  })
  
  observeEvent(input$confirmOverwrite, {
    removeModal()
    add_tables(input, vals)
    #showNotification("Existing table overwritten.")
  })

  observeEvent(input$AnnotationDelimiter, {
    if(input$AnnotationDelimiter == "custom") {
      shinyjs::show(id = "CustomAnnotDelim")
    } else {
      shinyjs::hide(id = "CustomAnnotDelim")
    }
  })
  
  output$annotations <- renderDataTable({
    file <- input$AnnotationsFile
    ext <- tools::file_ext(file$datapath)
    req(file)
    delim = input$AnnotationDelimiter
    if (delim == "custom") {
      delim = input$CustomAnnotDelim
    }
    vals$annotations <- read.delim(file$datapath, 
                                   header = input$AnnotationHeader,
                                   sep = delim,
                                   as.is = TRUE)
    vals$annotations
    
  }, options = list(scrollX = T))
  
  output$AnnotationSamples <- renderUI({
    if(is.null(vals$annotations)) {
      return (NULL) 
    }
    tagList(
      selectInput("AnnotSampleColumn", "Sample Name Column", 
                  choices = colnames(vals$annotations))    
      )
  })
  
  #Add Annotations to Musica object
  observeEvent(input$AddAnnotation, {
    if (!is.null(getResult(input$AnnotationMusicaList))) {
      tryCatch( {
      new_annot <- merge(samp_annot(getResult(input$AnnotationMusicaList)),
                         vals$annotations, by.x = "Samples", 
                         by.y = input$AnnotSampleColumn,
                         all.x = T)
      sapply(names(new_annot), 
             FUN = function(a) {
          samp_annot(vals$result_objects[[input$AnnotationMusicaList]], a) <- 
            new_annot[,a]
        })
        showNotification("Annotations have been added")
        }, error = function(cond) {
          shinyalert::shinyalert(title = "Error", text = cond$message)
          return()
        })
    } else if (!is.null(vals$musica)) {
      tryCatch( {
      new_annot <- merge(samp_annot(vals$musica),
                         vals$annotations, by.x = "Samples",
                         by.y = input$AnnotSampleColumn,
                         all.x = T)
      sapply(names(new_annot),
             FUN = function(a) {
               samp_annot(vals$musica, a) <- new_annot[,a]
             })
        showNotification("Annotations have been added")
        }, error = function(cond) {
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
                  choices = list("cosmic signatures" =
                                   as.list(names(cosmic_objects)),
                                 "your signatures" = as.list(
                              names(vals$result_objects)))),
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
      other <- isolate(getResult(input$SelectResultB))
    }
    tryCatch({
      #shinybusy::show_spinner()
      isolate(vals$comparison <- compare_results(isolate(getResult(input$SelectResultA)),
                      other, threshold = as.numeric(input$Threshold),
                      metric = input$CompareMetric))
      #shinybusy::hide_spinner()
    }, error = function(cond) {
      #shinybusy::hide_spinner()
      shinyalert::shinyalert(title = "Error", text = cond$message)
    })
    colnames(vals$comparison) <- c(input$CompareMetric, 
                                   paste0(input$SelectResultA, "-Index"),
                                   paste0(input$SelectResultB, "-Index"),
                                   paste0(input$SelectResultA, "-Signature"),
                                   paste0(input$SelectResultB, "-Signature"))
                                   
    if(!is.null(isolate(vals$comparison))) {
      output$CompareTable <- renderDataTable({
        isolate(vals$comparison)
      }, options = list(scrollX = T))
      output$DownloadComparison <- renderUI({
        tagList(
          downloadButton("DownloadCompare", "Download"),
          bsTooltip("DownloadCompare",
                    "Download the comparison table", 
                    placement = "bottom", trigger = "hover", options = NULL)
        )
      })
    }
  })
  
  output$DownloadCompare <- downloadHandler(
    filename = function() { paste0("Sig-Compare-", Sys.Date(), ".csv")
      },
    content = function(file) {
      write.csv(vals$comparison, file)
    }
  )

  output$DiffAnalResult <- renderUI({
    tagList(
      selectInput("DiffAnalResult", "Result Object", 
                  choices = names(vals$result_objects)),
      bsTooltip("DiffResult", "Musica Result object")
    )
  })
  
  output$DiffAnalGroups <- renderUI({
    if(interactive()) {
      tagList(
        sortable::bucket_list(
          header = "Groups",
          #group_name = "diff_groups",
          orientation = "horizontal",
          add_rank_list(
            text = "Group 1",
            labels = unique(samp_annot(getResult(input$DiffAnalResult))[[input$DiffAnalAnnot]]),
            input_id = "DiffGroup1"
          ),
          add_rank_list(
            text = "Group2",
            input_id = "DiffGroup2"
          )
        )
      )
    }
  })
  
  output$DiffAnalAnnot <- renderUI({
    tagList(
      selectInput("DiffAnalAnnot", "Sample annotation", 
                choices = colnames(samp_annot(getResult(input$DiffAnalResult)))[-1],
                selected = 1),
      bsTooltip("DiffAnalAnnot", "Sample annotation used to run differential analysis.")
    )
  })
  
  
  observeEvent(input$DiffMethod, {
    method = input$DiffMethod
    if (method == "wilcox") {
      shinyjs::show(id = "DiffAnalGroups")
    } else {
      shinyjs::hide(id = "DiffAnalGroups")
    }
  })
  
  observeEvent(input$RunDiffAnal, {
    shinyjs::hide("DiffError")
    g1 <- input$DiffGroup1
    g2 <- input$DiffGroup2
    errors <- NULL
    if(!is.null(g1) & !is.null(g2)) {
      gMin <- min(length(g1), length(g2))
      g1 <- g1[1:gMin]
      g2 <- g2[1:gMin]
    }
    #shinybusy::show_spinner()
    tryCatch({
      vals$diff <- exposure_differential_analysis(getResult(input$DiffAnalResult),
                                   input$DiffAnalAnnot,
                                   method = input$DiffMethod,
                                   group1 = g1, 
                                   group2 = g2)
      output$DiffTable <- renderDataTable(
        vals$diff %>% tibble::rownames_to_column(var = "Signature"), 
        options = list(scrollX = T)
      )
      shinybusy::hide_spinner()
      # output$DownloadDiffAnal <- renderUI({
      #   tagList(
      #     downloadButton("DownloadDiff", "Download"),
      #     bsTooltip("DownloadDiff",
      #               "Download the differential exposure table",
      #               placement = "bottom", trigger = "hover", options = NULL)
      #   )
      # })
    }, error = function(cond) {
      shinybusy::hide_spinner()
      output$DiffTable <- renderDataTable({NULL})
      errors <- cond
    })
    output$DiffError <- renderText({
      errors
    })
  })
  
  output$DownloadDiff <- downloadHandler(
    filename = function() { paste0("Exp-Diff-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(vals$diff, file)
    }
  )

  getResult <- function(name) {
    return(vals$result_objects[[name]])
  }
  
###############################################################################
 
##################Visualization#################   
  output$select_res1 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res1",
        label = "Select Result",
        choices = c(names(vals$result_objects)),
        width = "50%"
      ),
      bsTooltip(id = "selected_res1", title = "Select one musica_result object to visualize signatures.",
                 placement = "right", options = list(container = "body"))
    )
  })
  
  output$select_res2 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res2",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      ),
      bsTooltip(id = "selected_res2", title = "Select one musica_result object to visualize exposures.",
                placement = "right", options = list(container = "body"))
    )
  })
  
  # observeEvent(input$get_res,{
  #   vals$data <- res_annot
  # })
  
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
    plotly <- input$plotly1
    options <- list(legend, text_size, facet_size, show_x_labels, same_scale, plotly)
    return(options)
  }
  
  observeEvent(input$get_plot1,{
    options <- get_sig_option(input)
    result <- vals$result_objects[[input$selected_res1]]
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
    height <- paste0(as.character(n * 75),"px")
    if(options[[6]]){
      removeUI(selector = "#sigplot_plot")
      removeUI(selector = "#sigplot_plotly")
      insertUI(
        selector = "#plotdiv1",
        ui = plotlyOutput(outputId = "sigplot_plotly", height = height)
      )
      output$sigplot_plotly <- renderPlotly(
        plot_signatures(
          result = result, 
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
      removeUI(selector = "#sigplot_plot")
      insertUI(
        selector = "#plotdiv1",
        ui = plotOutput(outputId = "sigplot_plot", height = height)
      )
      output$sigplot_plot <- renderPlot(
        plot_signatures(
          result = result, 
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
          numericInput(inputId = "pointsize", label = "Point Size", value = 2),
          bsTooltip(id = "addpoint", title = "If checked, then points for individual sample exposures will be 
                    plotted on top of the violin/box plots.",
                    placement = "right", options = list(container = "body")),
          bsTooltip(id = "pointsize", title = "Size of the points to be plotted on top of the violin/box plots.",
                    placement = "right", options = list(container = "body"))
        )
      )
      removeUI(selector = "#group1")
      insertUI(
        selector = "#insert_group",
        ui = tagList(
               radioButtons(
                 inputId = "group1", 
                 label = "Group By", 
                 choices = list("Signature" = "signature",
                                "Annotation" = "annotation"),
                 inline = TRUE,
                 selected = "signature"
               ),
               bsTooltip(id = "group1", title = "Determines how to group samples into the subplots. If set to \"annotation\", then a sample annotation must be supplied via the annotation parameter.",
                         placement = "right", options = list(container = "body"))
             )
      )
      vals$point_ind <- 1
    }
    else if(input$plottype == "bar"){
      removeUI(selector = "#point_opt")
      removeUI(selector = "#group1")
      insertUI(
        selector = "#insert_group",
        ui = tagList(
          radioButtons(
            inputId = "group1",
            label = "Group By",
            choices = list("None" = "none","Signature" = "signature",
                           "Annotation" = "annotation"),
            inline = TRUE,
            selected = "none"
          ),
          bsTooltip(id = "group1", title = "Determines how to group samples into the subplots. If set to \"annotation\", then a sample annotation must be supplied via the annotation parameter.",
                    placement = "right", options = list(container = "body"))
        )
      )
      vals$point_ind <- 0
    }
    else{
      return(NULL)
    }
  })
  
  addTooltip(session, id = "proportional", title = "If checked, the exposures will be normalized to 
             between 0 and 1 by dividing by the total number of counts for each sample.",
             placement = "right", options = list(container = "body"))
  addTooltip(session, id = "color", title = "Determines how to color the bars or box/violins. If set to \"annotation\", 
             then a sample annotation must be supplied via the annotation parameter", 
             placement = "right", options = list(container = "body"))
  
  observeEvent(input$group1,{
    if(input$group1 == "annotation" & input$color != "annotation"){
      if(ncol(samp_annot(vals$result_objects[[input$selected_res2]])) == 1){
        shinyalert::shinyalert(title = "Error", text = "Annotation not found. Please add annotation to the musica object.")
      }
      else{
        vals$annot <- as.list(colnames(samp_annot(vals$result_objects[[input$selected_res2]]))[-1])
        names(vals$annot) <- colnames(samp_annot(vals$result_objects[[input$selected_res2]]))[-1]
        insertUI(
          selector = "#insertannot",
          ui = tagList(
            selectInput(
              inputId = "annotation",
              label = "Annotation",
              choices = vals$annot
            ),
            bsTooltip(id = "annotation", title = "Sample annotation used to group the subplots or color the bars, boxes, or violins.",
                      placement = "right", options = list(container = "body"))
          )
        )
      }
    }
    else{
      if(input$color != "annotation"){
        removeUI(selector = "div:has(>> #annotation)")
      }
    }
  })
  
  observeEvent(input$color,{
    if(input$color == "annotation" & input$group1 != "annotation"){
      if(ncol(samp_annot(vals$result_objects[[input$selected_res2]])) == 1){
        shinyalert::shinyalert(title = "Error", text = "Annotation not found. Please add annotation to the musica object.")
      }
      else{
        insertUI(
          selector = "#insertannot",
          ui = tagList(
            selectInput(
              inputId = "annotation",
              label = "Annotation",
              choices = vals$annot
            ),
            bsTooltip(id = "annotation", title = "Sample annotation used to group the subplots or color the bars, boxes, or violins.",
                      placement = "right", options = list(container = "body"))
          )
        )
      }
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
          ),
          bsTooltip(id = "bucket", title = "Drag signatures from top bucket to the bottom bucket. 
                    Samples will be sorted in descending order by signatures. If multiple signatures are supplied, 
                    samples will be sorted by each signature sequentially",
                    placement = "right", options = list(container = "body"))
        )
      )
    }
    else{
      removeUI(selector = "#insertsig")
    }
  })
  
  observeEvent(input$selected_res2, {
    output$number <- renderUI(
      tagList(
        numericInput(inputId = "numsamp", label = "# of Top Samples", 
                     value = dim(vals$result_objects[[input$selected_res2]]@exposures)[2],
                     min = 1,
                     max = dim(vals$result_objects[[input$selected_res2]]@exposures)[2]),
        bsTooltip(id = "numsamp", title = "The top number of sorted samples to display.",
                  placement = "right", options = list(container = "body"))
      )
    )
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
    result <- vals$result_objects[[input$selected_res2]]
    if(options[[14]]){
      removeUI(selector = "#expplotly")
      removeUI(selector = "#expplot")
      insertUI(
        selector = "#plotdiv2",
        ui = plotlyOutput(outputId = "expplotly")
      )
      output$expplotly <- renderPlotly(
        plot_exposures(
          result = result, 
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
      removeUI(selector = "#expplot")
      removeUI(selector = "#expplotly")
      insertUI(
        selector = "#plotdiv2",
        ui = plotOutput(outputId = "expplot")
      )
      output$expplot <- renderPlot(
        plot_exposures(
          result = result, 
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
  # observeEvent(input$def_sigs,{
  #   insertUI(
  #     ui = tags$div("#sortbysigs"),
  #   radioButtons(
  #     inputId = "subset",
  #     label = "",
  #     choices = c("All Signatures" = "all_sigs","Selected Signatures" = "signature"),
  #     inline = TRUE,
  #     selected = ""
  #   )
  #   )
  #   })
  
  # observeEvent(input$def_sigs, {
  #   insertUI()
  # })
  
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
              input_id = "sort_sigs"
            )
          )
        )
      )
    #vals$sort_sigs <- "#sort_sigs"
    }
    else if(input$subset == "all_signatures"){
	      removeUI(selector = "#insertsig")
        #vals$sort_sigs <- NULL
      }
  })
  get_sigs <- function(input){
    req(input$subset)
    if(input$subset == "signature"){
      sig <- input$sort_sigs
    }
    else if(input$subset == "all_signatures"){
      sig <- NULL
    }
    return(sig)
  }
  
  observeEvent(input$subset_tum, {
    if(input$subset_tum == "tumors"){
      insertUI(
        selector = "#sortbytum",
        ui = tags$div(
          id = "inserttum",
          #selectInput("tum_val","",choices = as.list(unique(vals$result_objects[[input$select_res_heatmap]]@musica@sample_annotations$Tumor_Subtypes)))
          checkboxGroupInput("tum_val","Available Samples:",
                             as.list(unique(vals$result_objects[[input$select_res_heatmap]]@musica@sample_annotations$Tumor_Subtypes)))
          )
        )

    }
    else{
      removeUI(selector = "#inserttum")
    }
  })
  #subset_tumor = input$tum_val
  
  observeEvent(input$subset_annot, {
    if(input$subset_annot == "annotation"){
      insertUI(
        selector = "#sortbyannot",
        ui = tags$div(
          id = "#insertannot",
          checkboxGroupInput("annot_val","Available annotations:",
                             c(as.list(colnames(samp_annot(vals$result_objects[[input$select_res_heatmap]]))))
          #selectInput("annot_val","",choices = as.list(colnames(samp_annot(vals$result_objects[[input$select_res_heatmap]])))
        )
      ))
    }
    else{
      removeUI(selector = "#insertannot")
    }
  })
observeEvent(input$get_heatmap,{
  sigs <- get_sigs(input)
  output$heatmap <- renderPlot({
    req(input$select_res_heatmap)
    
    #paste0("Heatmap")
    # propor()
    # sel_col_names()
    # sel_row_names()
    # zscale()
    # input$sig_to
    # input$tum_val
    # input$annot_val
    input$get_heatmap
    isolate(plot_heatmap(res_annot = vals$result_objects[[input$select_res_heatmap]],proportional = propor(),show_row_names = sel_row_names(),show_column_names = sel_col_names(),scale = zscale(),subset_signatures = c(sigs),subset_tumor = input$tum_val,annotation = input$annot_val,column_title = paste0("Heatmap for ",input$select_res_heatmap)))
  })
 }) 
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("heatmap", ".png", sep = "")
    },
    content = function(file) {
      ggsave(file,plot = last_plot(),device = "png")
    }
  )
 
 ##############Clustering################
  output$select_res3 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res3",
        label = h3("Select Result"),
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
        selector = "#iter",
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
      ui = tags$div(
        DT::dataTableOutput("cluster_table"),
        downloadButton("download_cluster", "Download")
      )
    )
    annot <- samp_annot(vals$result_objects[[input$selected_res3]])
    row.names(annot) <- annot$Samples
    dat <- cbind(annot, vals$cluster)
    output$cluster_table <- DT::renderDataTable(
      DT::datatable(dat[-1])
    )
    output$download_cluster <- downloadHandler(
      filename = function() {
        paste0(input$selected_res3, "_cluster.txt")
      },
      content = function(file){
        write.table(dat[-1], file, sep = '\t', quote = F)
      }
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
    plotly <- input$plotly3
    if(plotly){
      removeUI(selector = "#cluster_plot")
      removeUI(selector = "#cluster_plotly")
      insertUI(
        selector = "#clusterplotdiv",
        ui = plotlyOutput(outputId = "cluster_plotly")
      )
      output$cluster_plotly <- renderPlotly(
        plot_cluster(result = result,
                     clusters = clusters,
                     group = group,
                     annotation = annotation,
                     plotly = plotly)
      )
    }
    else{
      removeUI(selector = "#cluster_plot")
      removeUI(selector = "#cluster_plotly")
      insertUI(
        selector = "#clusterplotdiv",
        ui = plotOutput(outputId = "cluster_plot")
      )
      output$cluster_plot <- renderPlot(
        plot_cluster(result = result,
                     clusters = clusters,
                     group = group,
                     annotation = annotation,
                     plotly = plotly)
      )
    }
  })
  ########################################
  ########################################Download########################################
  output$select_mus_obj_download <- renderUI({
    tagList(
      selectInput(
        inputId = "select_mus_obj_download",
        label = "Select Musica Object",
        choices = "musica"
      )
    )
  })
  
  output$select_res_download <- renderUI({
    tagList(
      selectInput(
        inputId = "select_res_download",
        label = "Select Musica Result Object",
        choices = c(names(vals$result_objects))
      )
    )
  })
  
  output$download_mus_obj <- downloadHandler(
    
    filename = function() {
      paste("musica_object", ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(vals$musica, file)
    }
  )
 output$download_res <- downloadHandler(
    
    filename = function() {
      paste("musica_results", ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(vals$result_objects[[input$select_res_download]], file)
    }
  )
}






