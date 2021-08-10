library(musicatk)
library(plotly)
library(sortable)
library(shinyBS)
library(shinyalert)
library(TCGAbiolinks)
library(shinyjqui)

options(shiny.maxRequestSize = 10000 * 1024 ^ 2)
source("server_tables.R", local = T)

server <- function(input, output, session) {
#################### GENERAL ##################################################
  
  #Deactivate all tabs except import
  shinyjs::addCssClass(selector = "a[data-value='musica']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='annotations']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='tables']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='discover']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='predict']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='visualization']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='compare']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='differentialanalysis']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='cluster']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='heatmap']",
                     class = "inactiveLink")
  shinyjs::addCssClass(selector = "a[data-value='download']",
                     class = "inactiveLink")

  #Initialize variables
  vals <- reactiveValues(
    genome = NULL,
    musica = NULL,
    files = NULL,
    result_objects = list(),
    p_sigs = NULL, #cosmic sig
    p_res = NULL, #cosmic result object
    annotations = NULL,
    diff = NULL,
    df = NULL,
    musica_upload = NULL,
    data = NULL,
    point_ind = 0, #indicator for point option in exposure
    annot = NULL,
    deleted_rows = NULL,
    deleted_row_indices = list(),
    cluster = NULL,
    var = NULL,
    musica_name_user = NULL,
    musica_message = NULL,
    sort_sigs = NULL
  )

  #Control flow to detect user's progress and active tabs
  observeEvent(input$menu, {
    if (input$menu == "musica") {
      if (is.null(vals$var)) {
        shinyalert::shinyalert("Error",
                               paste0("No data was uploaded. ",
                                      "Please go to \"Import\" ",
                                      "and upload your data.", "error"))
        updateTabItems(session, "menu", "import")
      }
      else{
        removeCssClass(selector = "a[data-value='musica']",
                       class = "inactiveLink")
      }
    }
    else if (input$menu == "annotations") {
      if (is.null(vals$musica) && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No musica object was found. ",
                                      "Please go to \"Import\" or ",
                                      "\"Create Musica Object\" to upload or ",
                                      "create an object.", "error"))
        updateTabItems(session, "menu", "import")
      }
    }
    else if (input$menu == "download") {
      if (is.null(vals$musica) && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No musica object was found. ",
                                      "Please go to \"Import\" or ",
                                      "\"Create Musica Object\" to upload or ",
                                      "create an object.", "error"))
        updateTabItems(session, "menu", "import")
      }
    }
    else if (input$menu == "tables") {
      if (is.null(vals$var) && is.null(vals$musica) &&
          is.null(vals$musica_upload)) {
        shinyalert::shinyalert("Error",
                               paste0("No data was uploaded. ",
                                      "Please go to \"Import\" and upload ",
                                      "your data.", "error"))
        updateTabItems(session, "menu", "import")
      }
      else if (!is.null(vals$var) && is.null(vals$musica) &&
               is.null(vals$musica_upload)) {
        shinyalert::shinyalert("Error",
                               paste0("No musica object was created. ",
                                      "Please go to \"Create Musica Object\", ",
                                      "to create a musica object or go to ",
                                      "\"Import\" -> \"Import Musica Result ",
                                      "Object\" to upload a muscia object."),
                               "error")
        updateTabItems(session, "menu", "musica")
      }
      else{
        removeCssClass(selector = "a[data-value='tables']",
                       class = "inactiveLink")
      }
    }
    else if (input$menu %in% c("discover", "predict")) {
      if (is.null(vals$var) && is.null(vals$musica) &&
          is.null(vals$musica_upload)) {
        shinyalert::shinyalert("Error",
                               paste0("No data was uploaded. ",
                                      "Please go to \"Import\" and upload ",
                                      "your data."), "error")
        updateTabItems(session, "menu", "import")
      }
      else if (!is.null(vals$var) && is.null(vals$musica) &&
               is.null(vals$musica_upload)) {
        shinyalert::shinyalert("Error",
                               paste0("No musica object was created.",
                               "Please go to \"Create Musica Object\",",
                               "to create a musica object or go to ",
                               "\"Import\" -> \"Import Musica Result Object\" ",
                               "to upload a muscia object."), "error")
        updateTabItems(session, "menu", "musica")
      }
      else if (!is.null(vals$var) && !is.null(vals$musica) &&
               length(tables(vals$musica)) == 0 &&
               is.null(vals$musica_upload)) {
        shinyalert::shinyalert("Error",
                               paste0("No mutation count table was created. ",
                                      "Please go to \"Build Tables\" to ",
                                      "create count table."),
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if (!is.null(vals$var) && is.null(vals$musica) &&
               !is.null(vals$musica_upload) &&
               length(tables(vals$musica_upload)) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No mutation count table was created. ",
                                      "Please go to \"Build Tables\" to ",
                                      "create count table."),
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else{
        removeCssClass(selector = "a[data-value='discover']",
                       class = "inactiveLink")
        removeCssClass(selector = "a[data-value='predict']",
                       class = "inactiveLink")
      }
    }
    else if (input$menu %in% c("visualization", "compare",
                               "differentialanalysis", "cluster", "heatmap")) {
      if (is.null(vals$var) && is.null(vals$musica) &&
          is.null(vals$musica_upload) && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No data was uploaded. ",
                                      "Please go to \"Import\" and upload ",
                                      "your data."), "error")
        updateTabItems(session, "menu", "import")
      }
      else if (!is.null(vals$var) && is.null(vals$musica) &&
               is.null(vals$musica_upload) &&
               length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No musica object was created. ",
                                      "Please go to \"Create Musica Object\" ",
                                      "to create a musica object or go to ",
                                      "\"Import\" -> \"Import Musica Result ",
                                      "Object\ to upload a muscia object."),
                               "error")
        updateTabItems(session, "menu", "musica")
      }
      else if (!is.null(vals$var) && !is.null(vals$musica) &&
               length(tables(vals$musica)) == 0 &&
               is.null(vals$musica_upload) &&
               length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No mutation count table was created. ",
                                      "Please go to \"Build Tables\" ",
                                      "to create count table."),
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if (!is.null(vals$var) && is.null(vals$musica)
               && !is.null(vals$musica_upload)
               && length(tables(vals$musica_upload)) == 0
               && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No mutation count table was created. ",
                                      "Please go to \"Build Tables\" to ",
                                      "create count table."),
                               "error")
        updateTabItems(session, "menu", "tables")
      }
      else if (!is.null(vals$var) && (!is.null(vals$musica) ||
                                      !is.null(vals$musica_upload))
               && (length(tables(vals$musica)) != 0 ||
                   length(tables(vals$musica_upload)) != 0)
               && length(vals$result_objects) == 0) {
        shinyalert::shinyalert("Error",
                               paste0("No musica_result object was created. ",
                                      "Please go to \"Signatures and ",
                                      "Exposures\" -> \"Discover Signatures ",
                                      "and Exposures\" to create ",
                                      "musica_result object."), "error")
        updateTabItems(session, "menu", "discover")
      }
      else{
        removeCssClass(selector = "a[data-value='visualization']",
                       class = "inactiveLink")
        removeCssClass(selector = "a[data-value='compare']",
                       class = "inactiveLink")
        removeCssClass(selector = "a[data-value='differentialanalysis']",
                       class = "inactiveLink")
        removeCssClass(selector = "a[data-value='cluster']",
                       class = "inactiveLink")
        removeCssClass(selector = "a[data-value='heatmap']",
                       class = "inactiveLink")
      }
    }
    else{
      return()
    }
  })

  #If variants uplodaed, enable musica tab
  observeEvent(vals$var, {
    if (!is.null(vals$var)) {
      removeCssClass(selector = "a[data-value='musica']",
                     class = "inactiveLink")
    }
  })

  #If musica object present, enable downstream tabs
  observeEvent(vals$musica, {
    if (!is.null(vals$musica)) {
      removeCssClass(selector = "a[data-value='tables']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='annotations']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='download']",
                     class = "inactiveLink")
    }
    if (length(tables(vals$musica)) != 0) {
      removeCssClass(selector = "a[data-value='discover']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='predict']",
                     class = "inactiveLink")
    }
  })

  #If there's uploaded musica object, enable downstream tabs
  observeEvent(vals$musica_upload, {
    if (!is.null(vals$musica_upload)) {
      removeCssClass(selector = "a[data-value='tables']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='annotations']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='download']",
                     class = "inactiveLink")
    }
    if (length(tables(vals$musica_upload)) != 0) {
      removeCssClass(selector = "a[data-value='discover']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='predict']",
                     class = "inactiveLink")
    }
  })

  #If there's result object, enable downstream tabs
  observeEvent(vals$result_objects, {
    if (length(vals$result_objects) > 0) {
      removeCssClass(selector = "a[data-value='annotations']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='visualization']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='compare']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='differentialanalysis']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='cluster']",
                     class = "inactiveLink")
      removeCssClass(selector = "a[data-value='heatmap']",
                     class = "inactiveLink")
    }
  })

###################### Zainab's Code ##########################################
  output$tcga_tumor <- renderUI({
    hr()
    #Extracting TCGA tumors and formatting it to display only the abbreviations
    p <- TCGAbiolinks:::getGDCprojects()$project_id
    t <- grep("TCGA", p)
    p <- p[t]
    p <- gsub(".*-", "", p)
    names(p) <- c("UCEC - Uterine Corpus Endometrial Carcinoma",
                  "LGG - Brain Lower Grade Glioma",
                  "SARC - Sarcoma", "PAAD - Pancreatic adenocarcinoma",
                  "ESCA - Esophageal carcinoma",
                  "PRAD - Prostate adenocarcinoma",
                  "LAML - Acute Myeloid Leukemia",
                  "KIRC - Kidney renal clear cell carcinoma",
                  "PCPG - Pheochromocytoma and Paraganglioma",
                  "HNSC - Head and Neck squamous cell carcinoma",
                  "OV - Ovarian serous cystadenocarcinoma",
                  "GBM - Glioblastoma multiforme",
                  "UCS - Uterine Carcinosarcoma", "MESO - Mesothelioma",
                  "TGCT - Testicular Germ Cell Tumors",
                  "KICH - Kidney Chromophobe", "READ - Rectum adenocarcinoma",
                  "UVM - Uveal Melanoma",
                  "THCA - Thyroid carcinoma",
                  "LIHC - Liver hepatocellular carcinoma", "THYM - Thymoma",
                  "CHOL - Cholangiocarcinoma",
                  "DLBC - Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
                  "KIRP - Kidney renal papillary cell carcinoma",
                  "BLCA - Bladder Urothelial Carcinoma",
                  "BRCA - Breast invasive carcinoma",
                  "COAD - Colon adenocarcinoma",
                  paste0("CESC - Cervical squamous cell carcinoma ",
                  "and endocervical adenocarcinoma"),
                  "LUSC - Lung squamous cell carcinoma",
                  "STAD - Stomach adenocarcinoma",
                  "SKCM - Skin Cutaneous Melanoma",
                  "LUAD - Lung adenocarcinoma",
                  "ACC - Adrenocortical carcinoma")
    textInput("tcga_tumor", "Enter TCGA tumor type")
    tags$style("#tcga_tumor {
                    font-size:8px;
                    height:10px;
           }")
    checkboxGroupInput("tcga_tumor", " ", choices = as.list(p))
  })
  observeEvent(input$import_tcga, {
    req(input$import_tcga)
    if (!is.null(input$tcga_tumor)) {
      shinybusy::show_spinner()
      maf <- TCGAbiolinks::GDCquery_Maf(input$tcga_tumor, pipelines = "mutect")
      vals$var <- extract_variants_from_maf_file(maf)
      showNotification("TCGA dataset successfully imported!")
      shinybusy::hide_spinner()
    }
    else if (is.null(input$tcga_tumor)) {
      shinyalert("Error: No tumor found. Please select a tumor from the list!")
      }
  })
#Displaying thr table for TCGA tumor variants
  output$tcga_contents <- renderDataTable({
    req(vals$var)
    req(input$tcga_tumor)
    return(head(vals$var))
    shinyjs::show(id = "tcga_contents")
    js$enableTabs()
  })

  observeEvent(input$import, {
    req(input$file)
    file_name <- vals$data$datapath
    if (all(tools::file_ext(file_name) != c("maf", "vcf"))) {
      shinyalert::shinyalert(paste0("Error: File format not supported! ",
                                    "Please upload .maf or .vcf files"))
    }
    else{
        shinybusy::show_spinner()
        vals$var <- extract_variants(c(file_name))
        shinybusy::hide_spinner()
        showNotification("Import successfully completed!")
        }
  })
#Displaying list of available genomes
  output$genome_list <- renderUI({
    g <- BSgenome::available.genomes()
    g <- strsplit(g, ",")
    gg <- gsub("^.*?\\.", "", g)
    selectInput("GenomeSelect", "Choose genome:",
                list("Common genomes" =
                        list("hg38", "hg19", "hg18", "mm9", "mm10"),
                      "Genomes" = gg),
                width = "100%")
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

#Adding the create musica funcionality
  tryCatch({
    observeEvent(input$get_musica_object, {
      shinybusy::show_spinner()
      vals$genome <- input$GenomeSelect
      if (!is.null(vals$var)) {
        vals$musica <- create_musica(x = vals$var,
                                   genome = select_genome(vals$genome),
                                   check_ref_chromosomes = check_chr(),
                                   check_ref_bases = check_bases(),
                            convert_dbs = convert_dbs(),
                            standardize_indels = stand_indels())

      if (req(input$get_musica_object)) {
        shinyjs::show(id = "download_musica")
      }
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
  error = function(cond) {
    shinyalert::shinyalert(title = "Error", text = cond$message)
  },
  warning = function(cond) {
    shinyalert::shinyalert(title = "Error", text = cond$message)
  })

 output$musica_console <- renderPrint({
   return(print(vals$musica_message))
 })
#Displaying musica object variants
  output$musica_contents <- renderDataTable({
    req(vals$var)
    return(head(vals$var))
    shinyjs::show(id = "musica_contents")
    js$enableTabs()
  })

  output$musica_contents_summary <- renderText({
    req(vals$musica)
    vt <- unique(vals$musica@variants$Variant_Type) #variant types
    ns <- length(vals$musica@variants$sample) #sample length
    mylist <- c("No. of Samples:\n", ns)
    return(mylist)
    shinyjs::show(id = "musica_contents_summary")
    js$enableTabs();
  })
  output$musica_contents_table <- renderDataTable({
    req(vals$musica)
    nvt <- as.data.frame(table(vals$musica@variants$Variant_Type))
    return(nvt)
    shinyjs::show(id = "musica_contents_table")
    js$enableTabs();
    })
#Adding resetting functionality
  observeEvent(input$reset, {
    removeUI("#musica_contents")
    removeUI("#musica_contents_summary")
    removeUI("#musica_contents_table")
    showNotification("Tables cleared!")
  })
#Adding upload musica tab functionality
  observe(
    if (!is.null(req(input$musica_file))) {
    vals$musica_name_user <- tools::file_path_sans_ext(input$musica_file$name)
  })

  output$musica_result_name <- renderUI(
    textInput("musica_result_name", value = paste0(vals$musica_name_user),
              h3("Name your musica result object:"))
  )
  observeEvent(input$musica_button, {
    if (input$musica_button == "result") {
      shinyjs::show(id = "musica_result_name")
    }
    else if (input$musica_button == "object") {
      shinyjs::hide(id = "musica_result_name")
    }
  })
  observeEvent(input$upload_musica, {
    req(input$musica_file)
    if (all(tools::file_ext(tolower(input$musica_file$name)) !=
            c("rda", "rds"))) {
      shinyalert::shinyalert(paste0("Error: File format not supported! ",
      "Please upload .rda or .rds files"))
    }
    else{
      if (input$musica_button == "result") {
        if (all(tools::file_ext(tolower(input$musica_file$name)) == "rda")) {
          vals$musica_upload <- load(input$musica_file$datapath)
          vals$musica_upload <- get(vals$musica_upload)
          vals$result_objects[[input$musica_result_name]] <- vals$musica_upload
        }
        else if (all(tools::file_ext(tolower(input$musica_file$name)) ==
                     "rds")) {
          vals$musica_upload <- readRDS(input$musica_file$datapath)
          vals$result_objects[[input$musica_result_name]] <- vals$musica_upload
        }
        showNotification("Musica Result Object successfully imported!")
      }
      else if (input$musica_button == "object") {
        if (all(tools::file_ext(tolower(input$musica_file$name)) ==
                "rda")) {
          vals$musica_upload <- load(input$musica_file$datapath)
          vals$musica_upload <- get(vals$musica_upload)
          vals$musica <- vals$musica_upload
        }
        else if (all(tools::file_ext(tolower(input$musica_file$name)) ==
                     "rds")) {
          vals$musica_upload <- readRDS(input$musica_file$datapath)
          vals$musica <- vals$musica_upload
        }
        showNotification("Musica Object successfully imported!")
      }}
 })
#Displaying musica result/object summary table
  output$musica_upload <- renderDataTable({
    req(vals$musica_upload)
    return(head(vals$musica_upload@musica@variants))
    shinyjs::show(id = "musica_upload")
    js$enableTabs();
  })
  output$musica_upload_summary <- renderText({
    req(vals$musica_upload)
    vt <- unique(vals$musica_upload@musica@variants$Variant_Type) #variant types
    nvt <- table(vals$musica_upload@musica@variants$Variant_Type)
    ns <- length(vals$musica_upload@musica@variants$sample) #sample length
    mylist <- c("No. of Samples:\n", ns, "\n", "Variant types", vt, "\n", nvt)
    return(mylist)
    shinyjs::show(id = "musica_upload_summary")
    js$enableTabs();
  })

  observeEvent(input$reset_musica, {
    removeUI("#musica_upload")
    removeUI("#musica_upload_summary")
    showNotification("Tables cleared!")
  })

  #Adding download feature
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
    # Clear the previous deletions in the import table in the Import files tab
    vals$files <- list(input$file[["name"]])
    vals$files <- unlist(vals$files)
    vals$files <- c(vals$files)
    dt <- list(input$file[["datapath"]])
    dt <- unlist(dt)
    dt <- c(dt)
    vals$files <- data.frame(files = vals$files, datapath = dt,
                             stringsAsFactors = FALSE)
    vals$data <- vals$files
    vals$deleted_rows <- NULL
    vals$deleted_row_indices <- list()
    })

    observeEvent(input$deletep_ressed, {
    row_num <- parse_delete_event(input$deletep_ressed)
    data_row <- vals$data[row_num, ]
    vals$deleted_rows <- rbind(data_row, vals$deleted_rows)
    vals$deleted_row_indices <- append(vals$deleted_row_indices,
                                       row_num, after = 0)
    # Delete the row from the data frame
    vals$data <- vals$data[-row_num, ]
  })

  observeEvent(input$undo, {
    if (nrow(vals$deleted_rows) > 0) {
      row <- vals$deleted_rows[1, ]
      vals$data <- add_row_at(vals$data, row, vals$deleted_row_indices[[1]])
      # Remove row
      vals$deleted_rows <- vals$deleted_rows[-1, ]
      # Remove index
      vals$deleted_row_indices <- vals$deleted_row_indices[-1]
    }
  })

  #Disable the undo button if we have not deleted anything
  output$undo_ui <- renderUI({
    if (!is.null(vals$deleted_rows) && nrow(vals$deleted_rows) > 0) {
      actionButton("undo", label = "Undo delete", icon("undo"))
    } else {
      actionButton("undo", label = "Undo delete", icon("undo"), disabled = TRUE)
    }
  })

  output$dtable <- DT::renderDataTable({
    # Add the delete button column
    req(vals$data)
    delete_button_column(vals$data, "delete_button")
  })
#Code taken from an online open source
#' Adds a row at a specified index
#'
#' @param df a data frame
#' @param row a row with the same columns as \code{df}
#' @param i the index we want to add row at.
#' @return the data frame with \code{row} added to \code{df} at index \code{i}
add_row_at <- function(df, row, i) {
  if (i > 1) {
    rbind(df[1:(i - 1), ], row, df[- (1:(i - 1)), ])
  } else {
    rbind(row, df)
  }
}
#' A column of delete buttons for each row in the data frame for the
#' first column
#'
#' @param df data frame
#' @param id id prefix to add to each actionButton. The buttons will be
#'  id'd as id_INDEX.
#' @return A DT::datatable with escaping turned off that has the delete
#'  buttons in the first column and \code{df} in the other
delete_button_column <- function(df, id, ...) {
  # function to create one action button as string
  f <- function(i) {
    # https://shiny.rstudio.com/articles/communicating-with-js.html
    as.character(actionButton(paste(id, i, sep = "_"), label = NULL,
                              icon = icon("trash"),
                              onclick =
                              paste('Shiny.setInputValue(\"deletep_ressed\",',
                              'this.id, {priority: "event"})')))
  }
  delete_col <- unlist(lapply(seq(nrow(df)), f))
  # Return a data table
  DT::datatable(cbind(delete = delete_col, df),
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
parse_delete_event <- function(idstr) {
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

  # Renaming COSMIC signatures.
  sbs_aet <- list("SBS1" =
                    "SBS1 - Spontaneous deamination of 5-methylcytosine",
                  "SBS2" = "SBS2 - APOBEC activity",
                  "SBS3" = "SBS3 - HR deficiency",
                  "SBS4" = "SBS4 - Tobacco smoking",
                  "SBS5" =
                    "SBS5 - Unknown (Aging / Tobacco smoking / NER deficiency)",
                  "SBS6" = "SBS6 -  MMR deficiency",
                  "SBS7a" = "SBS7a - UV light exposure",
                  "SBS7b" = "SBS7b - UV light exposure",
                  "SBS7c" = "SBS7c - UV light exposure",
                  "SBS7d" = "SBS7d - UV light exposure",
                  "SBS8" = "SBS8 - Unknown (HR deficiency / NER deficiency)",
                  "SBS9" =
                    "SBS9 - Unknown (Polymerase eta somatic hypermutation)",
                  "SBS10a" = "SBS10a - POLE exonuclease domain mutation",
                  "SBS10b" = "SBS10b - POLE exonuclease domain mutation",
                  "SBS10c" = "SBS10c - Unknown (Defective POLD1 proofreading)",
                  "SBS10d" = "SBS10d - Unknown (Defective POLD1 proofreading)",
                  "SBS11" =
                    paste0("SBS11 - Unknown (Temozolomide chemotherapy ",
                    "/ MMR deficiency + temozolomide)"),
                  "SBS12" = "SBS12 - Unknown",
                  "SBS13" = "SBS13 - APOBEC activity",
                  "SBS14" = "SBS14 - MMR deficiency + POLE mutation",
                  "SBS15" = "SBS15 - MMR deficiency",
                  "SBS16" = "SBS16 - Unknown",
                  "SBS17a" = "SBS17a - Unknown (Damage by ROS)",
                  "SBS17b" =
                    "SBS17b - Unknown (Damage by ROS / 5FU chemotherapy)",
                  "SBS18" = "SBS18 - Damage by ROS",
                  "SBS19" = "SBS19 - Unknown",
                  "SBS20" = "SBS20 - MMR deficiency + POLD1 mutation",
                  "SBS21" = "SBS21 - MMR deficiency",
                  "SBS22" = "SBS22 - Aristolochic acid exposure",
                  "SBS23" = "SBS23 - Unknown",
                  "SBS24" = "SBS24 - Aflatoxin exposure",
                  "SBS25" = "SBS25 - Unknown (Unknown chemotherapy)",
                  "SBS26" = "SBS26 - MMR deficiency",
                  "SBS27" = "SBS27 - Sequencing artifact",
                  "SBS28" =
                    "SBS28 - Unknown (POLE exonuclease domain mutation)",
                  "SBS29" = "SBS29 - Unknown (Tobacco chewing)",
                  "SBS30" = "SBS30 - BER deficiency",
                  "SBS31" = "SBS31 - Platinum chemotherapy",
                  "SBS32" = "SBS32 - Azathioprine exposure",
                  "SBS33" = "SBS33 - Unknown (Unknown)",
                  "SBS34" = "SBS34 - Unknown",
                  "SBS35" = "SBS35 - Platinum chemotherapy",
                  "SBS36" = "SBS36 - BER deficiency",
                  "SBS37" = "SBS37 - Unknown",
                  "SBS38" =
                    "SBS38 - Unknown (UV light exposure (indirect effect))",
                  "SBS39" = "SBS39 - Unknown",
                  "SBS40" = "SBS40 - Unknown",
                  "SBS41" = "SBS41 - Unknown",
                  "SBS42" = "SBS42 - Haloalkanes exposure",
                  "SBS43" = "SBS43 - Unknown (Possible sequencing artifact)",
                  "SBS44" = "SBS44 - MMR deficiency",
                  "SBS45" =
                    "SBS45 - 8-oxo-guanine introduced during sequencing",
                  "SBS46" =
                    "SBS46 - Sequencing artifact (early releases of TCGA)",
                  "SBS47" =
                    paste0("SBS47 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS48" =
                    paste0("SBS48 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS49" =
                    paste0("SBS49 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS50" =
                    paste0("SBS50 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS51" = "SBS51 - Unknown (Possible sequencing artifact)",
                  "SBS52" =
                    paste0("SBS52 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS53" =
                    paste0("SBS53 - Sequencing artifact ",
                    "(blacklisted cancer samples for poor quality)"),
                  "SBS54" = "SBS54 - Germline variants contamination",
                  "SBS55" = "SBS55 - Unknown (Possible sequencing artifact)",
                  "SBS56" = "SBS56 - Unknown (Possible sequencing artifact)",
                  "SBS57" = "SBS57 - Unknown (Possible sequencing artifact)",
                  "SBS58" = "SBS58 - Unknown (Possible sequencing artifact)",
                  "SBS59" = "SBS59 - Unknown (Possible sequencing artifact)",
                  "SBS60" = "SBS60 - Sequencing artifact",
                  "SBS84" = "SBS84 - AID activity",
                  "SBS85" = "SBS85 - AID activity",
                  "SBS86" = "SBS86 - Unkown (Unknown chemotherapy)",
                  "SBS87" = "SBS87 - Thiopurine chemotherapy",
                  "SBS88" = "SBS88 - Colibactin exposure",
                  "SBS89" = "SBS89 - Unknown",
                  "SBS90" = "SBS90 - Duocarmycin exposure",
                  "SBS91" = "SBS91 - Unknown",
                  "SBS92" = "SBS92 - Unknown (Tobacco smoking)",
                  "SBS93" = "SBS93 - Unknown",
                  "SBS94" = "SBS94 - Unknown"
  )
  dbs_aet <- list("DBS1" = "DBS1 - UV light exposure",
                  "DBS2" =
                    "DBS2 - Unknown (Tobacco smoking / Acetaldehyde exposure)",
                  "DBS3" = "DBS3 - POLE exonuclease domain mutation",
                  "DBS4" = "DBS4 - Unknown",
                  "DBS5" = "DBS5 - Platinum chemotherapy",
                  "DBS6" = "DBS6 - Unknown",
                  "DBS7" = "DBS7 - MMR deficiency",
                  "DBS8" = "DBS8 - Unknown",
                  "DBS9" = "DBS9 - Unknown",
                  "DBS10" = "DBS10 - MMR deficiency",
                  "DBS11" = "DBS11 - Unknown (APOBEC activity)"
  )
  indel_aet <- list("ID1" =
                      "ID1 - Slippage of nascent strand during DNA replication",
                    "ID2" =
                      paste0("ID2 - Slippage of template strand during ",
                      "DNA replication"),
                    "ID3" = "ID3 - Tobacco smoking",
                    "ID4" = "ID4 - Unknown",
                    "ID5" = "ID5 - Unknown",
                    "ID6" = "ID6 - HR deficiency",
                    "ID7" = "ID7 - MMR deficiency",
                    "ID8" =
                      "ID8 - Unknown (DSB repair by NHEJ / TOP2A mutation)",
                    "ID9" = "ID9 - Unknown",
                    "ID10" = "ID10 - Unknown",
                    "ID11" = "ID11 - Unknown",
                    "ID12" = "ID12 - Unknown",
                    "ID13" = "ID13 - UV light exposure",
                    "ID14" = "ID14 - Unknown",
                    "ID15" = "ID15 - Unknown",
                    "ID16" = "ID16 - Unknown",
                    "ID17" = "ID17 - TOP2A mutation",
                    "ID18" = "ID18 - Colibactin exposure"
  )
  colnames(signatures(cosmic_v3_sbs_sigs)) <-
    sbs_aet[colnames(signatures(cosmic_v3_sbs_sigs))]
  colnames(signatures(cosmic_v3_dbs_sigs)) <-
    dbs_aet[colnames(signatures(cosmic_v3_dbs_sigs))]
  colnames(signatures(cosmic_v3_indel_sigs)) <-
    indel_aet[colnames(signatures(cosmic_v3_indel_sigs))]

  cosmic_objects <- list("cosmic_v3_sbs_sigs" = cosmic_v3_sbs_sigs,
                      "cosmic_v3_dbs_sigs" = cosmic_v3_dbs_sigs,
                      "cosmic_v3_indel_sigs" = cosmic_v3_indel_sigs)

  # Update musica object list whenever a result or musica object is altered
  output$annotation_musica_list <- renderUI({
    result_names <- list(names(vals$result_objects))
    if (is.null(vals$musica)) {
      tagList(
        selectInput("annotation_musica_list", "Select object",
                    choices = list("result objects" =
                                     list(names(vals$result_objects))))
      )
    } else if (is.null(result_names[[1]])) {
      tagList(
        selectInput("annotation_musica_list", "Select object",
                    choices = list("musica object" = list("musica")))
      )
    }
    else {
      tagList(
        selectInput("annotation_musica_list", "Select object",
                    choices = list("musica object" = list("musica"),
                                   "result objects" =
                                     list(names(vals$result_objects))))
      )
    }
  })

  # Input to choose the annotation delimiter.
  observeEvent(input$annotation_delimiter, {
    if (input$annotation_delimiter == "custom") {
      shinyjs::show(id = "CustomAnnotDelim")
    } else {
      shinyjs::hide(id = "CustomAnnotDelim")
    }
  })

  # Read the annotation file, store it in vals$annotations, and display the
  # contents as a data table.
  output$annotations <- renderDataTable({
    file <- input$annotations_file
    ext <- tools::file_ext(file$datapath)
    req(file)
    delim <- input$annotation_delimiter
    if (delim == "custom") {
      delim <- input$CustomAnnotDelim
    }
    vals$annotations <- read.delim(file$datapath,
                                   header = input$annotation_header,
                                   sep = delim,
                                   as.is = TRUE)
    vals$annotations

  }, options = list(autoWidth = FALSE, scrollX = TRUE))

  # Input to choose the column containing the sample names.
  output$annotation_samples <- renderUI({
    if (is.null(vals$annotations)) {
      return(NULL)
    }
    tagList(
      selectInput("annot_sample_column", "Sample Name Column",
                  choices = colnames(vals$annotations))
    )
  })

  # Add annotations to the provided musica_object
  add_annotation <- function(musica_object) {
    new_annot <- merge(samp_annot(musica_object),
                       vals$annotations, by.x = "Samples",
                       by.y = input$annot_sample_column,
                       all.x = T)
    sapply(names(new_annot),
           FUN = function(a) {
             samp_annot(musica_object, a) <-
               new_annot[, a]
           })
    showNotification("Annotations have been added")
  }

  # Event that triggers add_annotation function.
  observeEvent(input$add_annotation, {
    # Add annotation to result object
    if (!is.null(get_result(input$annotation_musica_list))) {
      tryCatch({
        add_annotation(vals$result_objects[[input$annotation_musica_list]])
      }, error = function(cond) {
        shinyalert::shinyalert(title = "Error", text = cond$message)
        return()
      })
    # Add annotation to musica object
    } else if (!is.null(vals$musica)) {
      tryCatch({
        add_annotation(vals$musica)
      }, error = function(cond) {
        shinyalert::shinyalert(title = "Error", text = cond$message)
        return()
      })
    } else {
      print("Error: selected object does not exist")
    }

  })

  ####### Section for the Build Table Tab #######
  # Input to select count table
  output$discover_table <- renderUI({
    tagList(
      selectInput("select_discover_table", "Select Count Table",
                  choices = names(
                    extract_count_tables(vals$musica))),
      bsTooltip("select_discover_table",
                "Name of the table to use for signature discovery.",
                placement = "right", trigger = "hover", options = NULL)
    )
  })

  # UI for inputs to combine tables.
  output$combine_table <- renderUI({
    if (length(names(extract_count_tables(vals$musica))) > 1) {
      tagList(
        box(width = 6,
            helpText(paste0("Combine any 2 or more tables contained in your ",
            "musica object. This is optional.")),
        checkboxGroupInput("combine_tables", "Tables to Combine",
                    choices = names(extract_count_tables(vals$musica))),
        textInput("combined_table_name", "Name of combined table"),
        uiOutput("combine_warning"),
        actionButton("Combine", "Build Combined Table"),
        bsTooltip("combined_table_name",
                  "Combine tables into a single table that can be used for
                  discovery/prediction.",
                  placement = "bottom", trigger = "hover", options = NULL),
        bsTooltip("Combine",
                  paste0("Combines tables into a single table that can ",
                  "be used for discovery/prediction."),
                  placement = "bottom", trigger = "hover", options = NULL),
        bsTooltip("combine_tables",
                  "Tables to combine.",
                  placement = "left", trigger = "hover", options = NULL),
        bsTooltip("combined_table_name",
                  "Name for the combined table.",
                  placement = "bottom", trigger = "hover", options = NULL)
        )
      )
    }

  })

  #
  observeEvent(input$Combine, {
    if (input$combined_table_name == "" | length(input$combine_tables) < 2) {
      output$combine_warning <- renderText({
        validate(
          need(input$combined_table_name != "",
               "You must provide a name for the new result object."),
          need(length(input$combine_tables) < 2,
               "You must select two or more tables to combine.")
        )
      })
      return()
    }
    shinybusy::show_spinner()
    tryCatch({
      combine_count_tables(vals$musica, input$combine_tables,
                         input$combined_table_name)
    }, error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
      shinybuys::hide_spinner()
    })
    shinybusy::hide_spinner()
    showNotification("Table created.")

  })

  output$allow_table <- renderUI({
    if (!is.null(vals$musica)) {
      tagList(
      actionButton("add_table", "Create Table"),
      bsTooltip("add_table",
                paste0("Create a table containig the mutationl count ",
                "information of each sample."),
                placement = "bottom", trigger = "hover",
                options = NULL)
      )
    } else {
      tagList(
      helpText("You must first create or upload a musica object to generate
               count tables.")
      )
    }
  })

  # Event listener for add table button.
  observeEvent(input$add_table, {
    table_name <- input$select_table
    if (table_name == "SBS192 - Replication_Strand") {
      table_name <- "SBS192_Rep"
    }
    if (table_name == "SBS192 - Transcript_Strand") {
      table_name <- "SBS192_Trans"
    }
    if (table_name %in% names(extract_count_tables(vals$musica))) {
      # Modal to confirm overwrite of existing table
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

  # Confirm overwrite for existing table
  observeEvent(input$confirmOverwrite, {
    removeModal()
    add_tables(input, vals)
  })

  ####### Section for the Discover Signatures and Exposures Tab #######
  # Event listener for discover_signatures.
  observeEvent(input$discover_signatures, {
    if (input$discover_result_name == "" |
        input$number_of_signatures == "" |
        input$n_start == "" |
        dim(extract_count_tables(vals$musica)[[
          input$select_discover_table]]@count_table)[2] < 2 |
        input$number_of_signatures < 2) {
      output$discover_warning <- renderText({
        validate(
          need(input$discover_result_name != "",
               "You must provide a name for the new result object."),
          need(input$number_of_signatures != "",
               "You must specify the number of expected signatures."),
          need(input$n_start != "",
               "Please specify the number of random starts."),
          need(input$number_of_signatures >= 2,
               "Must specify 2 or more signatures."),
          need(dim(extract_count_tables(vals$musica)[[
            input$select_discover_table]]@count_table)[2] > 2,
            "You must provide 2 or more samples")
        )
      })
      return()
    }
    if (input$discover_result_name %in% names(vals$result_objects)) {
      # Confirm overwrite of result object
      showModal(modalDialog(
        title = "Existing Result Object.",
        "Do you want to overwrite the existing result object?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirm_result_overwrite", "OK"),
          modalButton("Cancel"))
      ))
    } else {
      disc_sigs(input, vals)
      showNotification(paste0("Musica Result object, ",
                              input$discover_result_name,
                            ", was created"))
    }
  })

  # Wrapper function for discover_signature
  disc_sigs <- function(input, vals) {
    set_result(input$discover_result_name, discover_signatures(
      vals$musica, table_name = input$select_discover_table,
      num_signatures = as.numeric(input$number_of_signatures),
      algorithm = input$Method,
      #seed = input$Seed,
      nstart = as.numeric(input$n_start)))
  }

  # Run discover_signatures if overwrite confirmed.
  observeEvent(input$confirm_result_overwrite, {
    removeModal()
    disc_sigs(input, vals)
    showNotification("Existing result object overwritten.")
  })

  # Discover Musica Result Object
  output$discover_result_name <- renderUI({
    name <- input$select_discover_table
    tagList(
      textInput("discover_result_name", "Name for musica result object",
                value = paste0(name, "-Result")),
      bsTooltip("discover_result_name",
                "Name for the newly created musica result object.",
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })

  # UI to select musica object.
  output$discover_musica_list <- renderUI({
    tagList(
      selectInput("discover_musica_list", h3("Select Musica Object"),
                  choices = names(vals$result_objects))
    )
  })

  # Select counts table for prediction.
  output$predict_table <- renderUI({
    tagList(
      selectInput("select_pred_table", "Select Counts Table",
                  choices = names(extract_count_tables(vals$musica))),
      bsTooltip("select_pred_table",
                "Name of the table used for posterior prediction",
                placement = "right", trigger = "hover", options = NULL)
    )
  })

  # Event listener alters UI based on selected COSMIC table.
  observeEvent(input$cosmic_count_table, {
    if (input$cosmic_count_table == "SBS") {
      shinyjs::show(id = "cosmic_SBS_sigs")
      shinyjs::hide(id = "cosmic_DBS_sigs")
      shinyjs::hide(id = "cosmic_INDEL_sigs")
    } else if (input$cosmic_count_table == "DBS") {
      shinyjs::hide(id = "cosmic_SBS_sigs")
      shinyjs::show(id = "cosmic_DBS_sigs")
      shinyjs::hide(id = "cosmic_INDEL_sigs")
    } else {
      shinyjs::hide(id = "cosmic_SBS_sigs")
      shinyjs::hide(id = "cosmic_DBS_sigs")
      shinyjs::show(id = "cosmic_INDEL_sigs")
    }
  })

  # UI to name Predict result object
  output$predict_result_name <- renderUI({
    name <- names(extract_count_tables(vals$musica))[1]
    tagList(
      textInput("predict_result_name", "Name for musica result object",
        value = paste0(name, "-Predict-Result")),
      bsTooltip("predict_result_name",
                "Name for the newly created musica result object.",
                placement = "bottom", trigger = "hover", options = NULL)
    )
  })

  # UI to select which result objects to predict.
  output$predicted_result <- renderUI({
    other <- list(names(vals$result_objects))
    if (is.null(other[[1]])) {
      tagList(
        selectInput("predicted_result", "Result to Predict",
                    choices = list("Cosmic" = list(
                      "Cosmic V3 SBS Signatures" = "cosmic_v3_sbs_sigs",
                      "Cosmic V3 DBS Signatures" = "cosmic_v3_dbs_sigs",
                      "Cosmic V3 INDEL Signatures" = "cosmic_v3_indel_sigs")),
                    selected = "cosmic_v3_sbs_sigs"),
        bsTooltip("predicted_result",
                  "Result object containing the signatures to predict",
                  placement = "right", trigger = "hover", options = NULL)
      )    }
    else {
    tagList(
      selectInput("predicted_result", "Result to Predict",
                  choices = list("Cosmic" = list(
                    "Cosmic V3 SBS Signatures" = "cosmic_v3_sbs_sigs",
                    "Cosmic V3 DBS Signatures" = "cosmic_v3_dbs_sigs",
                    "Cosmic V3 INDEL Signatures" = "cosmic_v3_indel_sigs"),
                                 "your signatures" = other),
                  selected = "cosmic_v3_sbs_sigs"),
      bsTooltip("predicted_result",
                "Result object containing the signatures to predict",
                placement = "right", trigger = "hover", options = NULL)
    )
    }
  })

  # UI to select signatures to predict
  output$precited_signatures <- renderUI({
    if (input$predicted_result %in% names(cosmic_objects)) {
        vals$p_sigs <- colnames(signatures(
          cosmic_objects[[input$predicted_result]]))
        vals$p_res <- cosmic_objects[[input$predicted_result]]
    }
    else {
        vals$p_sigs <- colnames(signatures(
          vals$result_objects[[input$predicted_result]]))
        vals$p_res <- vals$result_objects[[input$predicted_result]]
    }
    tagList(
      shinyWidgets::dropdownButton(circle = FALSE, label = "Signatures",
                                   div(style =
                                       "max-height:80vh; overflow-y: scroll",
                                       checkboxGroupInput("pred_sigs", "",
                         choices = vals$p_sigs, inline = FALSE,
                         selected = vals$p_sigs))),
      bsTooltip("pred_sigs",
                "Signatures to predict.",
                placement = "right", trigger = "hover", options = NULL)
    )
  })

  # Event listener for Predict signatures
  observeEvent(input$predict_sigs, {
    if (input$predict_result_name == "" | length(input$pred_sigs) < 2) {
      output$predict_warning <- renderText({
        validate(
          need(input$predict_result_name != "",
               "You must provide a name for the new result object."),
          need(length(c(input$pred_sigs)) >= 2,
               "You must select two or more signatures to predict.")
        )
      })
      return()
    }
    if (input$predict_result_name %in% names(vals$result_objects)) {
      showModal(modalDialog(
        title = "Existing Result Object.",
        "Do you want to overwrite the existing result object?",
        easyClose = TRUE,
        footer = list(
          actionButton("confirm_predict_overwrite", "OK"),
          modalButton("Cancel"))
      ))
    } else {
      get_predict(input, vals)
      showNotification(paste0("Musica result object, ",
                              input$predict_result_name,
                              ", was created"))
    }
  })

  # Event listener displays additional options for deconstructSigs algorithm.
  observeEvent(input$predict_algorithm, {
    if (input$predict_algorithm == "deconstructSigs") {
      shinyjs::show(id = "predict_genome_list")
    } else {
      shinyjs::hide(id = "predict_genome_list")
    }
  })

  # Wrapper function for predict_exposures
  get_predict <- function(inputs, vals) {
    set_result(input$predict_result_name,
      predict_exposure(vals$musica, g = select_genome(input$predict_genome_list),
                       table_name = input$select_pred_table,
                       signature_res = vals$p_res,
                       algorithm = input$predict_algorithm,
                       signatures_to_use = input$pred_sigs))
  }

  # Event triggers predict_exposures when user confirms overwrite.
  observeEvent(input$confirm_predict_overwrite, {
    removeModal()
    get_predict(inputs, vals)
    showNotification("Existing result overwritten.")
  })


  ####### Compare Tab ########
  output$compare_result_a <- renderUI({
    tagList(
      selectInput("select_result_a", "Select result object",
                  choices = c(names(vals$result_objects))),
      bsTooltip("select_result_a",
                "A musica result object",
                placement = "right", trigger = "hover", options = NULL)
    )
  })

  output$compare_result_b <- renderUI({
    tagList(
      selectInput("select_result_b", "Select comparison result object",
                  choices = list("cosmic signatures" =
                                   as.list(names(cosmic_objects)),
                                 "your signatures" = as.list(
                              names(vals$result_objects)))),
      bsTooltip("select_result_b",
                "A second musica result object",
                placement = "right", trigger = "hover", options = NULL)
    )
  })

  # Event listener triggers comparison.
  observeEvent(input$compare_results, {
    if (is.null(input$select_result_a) | input$select_result_a == "" |
        input$Threshold == "") {
      output$compare_validate <- renderText({
        validate(
          need(input$select_result_a != "",
               "Please select a result object to compare."),
          need(input$Threshold == "",
               "Please provide a similarity threshold from 0 to 1.")
        )
      })
      return()
    }
    # Retreive either cosmic or custom result objects.
    if (input$select_result_b %in% names(cosmic_objects)) {
      other <- cosmic_objects[[input$select_result_b]]
    } else {
      other <- isolate(get_result(input$select_result_b))
    }
    # Attempt to compare signatures
    tryCatch({
      isolate(vals$comparison <-
                compare_results(isolate(get_result(input$select_result_a)),
                      other, threshold = as.numeric(input$Threshold),
                      metric = input$compare_metric))
    }, error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
    })
    colnames(vals$comparison) <- c(input$compare_metric,
                                   paste0(input$select_result_a, "-Index"),
                                   paste0(input$select_result_b, "-Index"),
                                   paste0(input$select_result_a, "-Signature"),
                                   paste0(input$select_result_b, "-Signature"))
    # generate table containig comparison statistics.
    if (!is.null(isolate(vals$comparison))) {
      output$compare_table <- renderDataTable({
        isolate(vals$comparison)
      }, options = list(autoWidth = FALSE, scrollX = T))
      output$download_comparison <- renderUI({
        tagList(
          downloadButton("download_compare", "Download"),
          bsTooltip("download_compare",
                    "Download the comparison table",
                    placement = "bottom", trigger = "hover", options = NULL)
        )
      })
    }
  })

  output$download_compare <- downloadHandler(
    filename = function() {
      paste0("Sig-Compare-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(vals$comparison, file)
    }
  )

  ######## Differential Analysis Tab #########
  output$diff_anal_result <- renderUI({
    tagList(
      selectInput("diff_anal_result", "Result Object",
                  choices = names(vals$result_objects)),
      bsTooltip("diff_anal_result", "Musica Result object")
    )
  })

  # UI for Wilcoxon Rank Sum Test bucket list
  output$diff_anal_groups <- renderUI({
    if (interactive()) {
      tagList(
        sortable::bucket_list(
          header = "Groups",
          orientation = "horizontal",
          add_rank_list(
            text = "Group 1",
            labels = unique(samp_annot(
              get_result(input$diff_anal_result))[[input$diff_anal_annot]]),
            input_id = "diff_group1"
          ),
          add_rank_list(
            text = "Group2",
            input_id = "diff_group2"
          )
        )
      )
    }
  })

  # UI to select sample annotation.
  output$diff_anal_annot <- renderUI({
    tagList(
      selectInput("diff_anal_annot", "Sample annotation",
                choices = colnames(
                  samp_annot(get_result(input$diff_anal_result)))[-1],
                selected = 1),
      bsTooltip("diff_anal_annot",
                "Sample annotation used to run differential analysis.")
    )
  })

  # Event handler controls Wilcoxon bucket list.
  observeEvent(input$diff_method, {
    method <- input$diff_method
    if (method == "wilcox") {
      shinyjs::show(id = "diff_anal_groups")
    } else {
      shinyjs::hide(id = "diff_anal_groups")
    }
  })

  # Event handler for differential analysis.
  observeEvent(input$run_diff_anal, {
    shinyjs::hide("diff_error")
    g1 <- input$diff_group1
    g2 <- input$diff_group2
    errors <- NULL
    if (!is.null(g1) & !is.null(g2)) {
      g_min <- min(length(g1), length(g2))
      g1 <- g1[1:g_min]
      g2 <- g2[1:g_min]
    }
    # Run differential analysis
    tryCatch({
      vals$diff <-
        exposure_differential_analysis(get_result(input$diff_anal_result),
                                       input$diff_anal_annot,
                                       method = input$diff_method,
                                       group1 = g1,
                                       group2 = g2)
      output$diff_table <- renderDataTable(
        if (input$diff_method == "wilcox") {
          vals$diff
        } else {
          vals$diff %>% tibble::rownames_to_column(var = "Signature")
        },
        options = list(autoWidth = FALSE, scrollX = TRUE)
      )
    }, error = function(cond) {
      output$diff_table <- renderDataTable({
        NULL
        })
      errors <- cond
    })
    output$diff_error <- renderText({
      errors
    })
  })

  # Download differential analysis results.
  output$download_diff <- downloadHandler(
    filename = function() {
      paste0("Exp-Diff-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(vals$diff, file)
    }
  )

  ####### Helper Functions #######
  # Set a result object
  set_result <- function(x, y) {
    vals$result_objects[[x]] <- y
  }

  get_result <- function(name) {
    return(vals$result_objects[[name]])
  }

###############################################################################

##################Visualization#################
  
  #select box for signature
  output$select_res1 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res1",
        label = "Select Result",
        choices = c(names(vals$result_objects)),
        width = "50%"
      ),
      bsTooltip(id = "selected_res1",
                title =
                  "Select one musica_result object to visualize signatures.",
                placement = "right", options = list(container = "body"))
    )
  })

  #select box for exposure
  output$select_res2 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res2",
        label = "Select Result",
        choices = c(names(vals$result_objects))
      ),
      bsTooltip(id = "selected_res2",
                title =
                  "Select one musica_result object to visualize exposures.",
                placement = "right", options = list(container = "body"))
    )
  })

  #create text input for changing signature names
  observeEvent(input$rename, {
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
    for (i in 1:n) {
      id <- paste0("sig", i)
      if (input$rename) {
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

  #save plotting options for signatures in a list
  get_sig_option <- function(input) {
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
    if (input$rename) {
      ids <- vector()
      for (i in 1:n) {
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
    options <- list(legend, text_size, facet_size,
                    show_x_labels, same_scale, plotly)
    return(options)
  }

  #Make signature plot
  observeEvent(input$get_plot1, {
    options <- get_sig_option(input)
    result <- vals$result_objects[[input$selected_res1]]
    n <- ncol(vals$result_objects[[input$selected_res1]]@signatures)
    height <- paste0(as.character(n * 90), "px")
    if (options[[6]]) {
      #disable resizable
      jqui_resizable("#sigplot_plotly", operation = "destroy")
      jqui_resizable("#sigplot_plot", operation = "destroy")
      #remove previous plot
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
      #enable resizable
      jqui_resizable("#sigplot_plotly")
    }
    else{
      jqui_resizable("#sigplot_plotly", operation = "destroy")
      jqui_resizable("#sigplot_plot", operation = "destroy")
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
      jqui_resizable("#sigplot_plot")
    }
  })

  observeEvent(input$plottype, {
    #If box or violin selected, generate options for plotting points;
    #group by and color by have two choices
    if (input$plottype %in% c("box", "violin") & vals$point_ind == 0) {
      removeUI(selector = "#point_opt")
      insertUI(
        selector = "#points",
        ui = tags$div(
          id = "point_opt",
          checkboxInput(inputId = "addpoint", label = "Add Points",
                        value = TRUE),
          numericInput(inputId = "pointsize", label = "Point Size", value = 2),
          bsTooltip(id = "addpoint",
          title =
          "If checked, then points for individual sample exposures will be
                    plotted on top of the violin/box plots.",
                    placement = "right", options = list(container = "body")),
          bsTooltip(id = "pointsize",
                    title = paste0("Size of the points to be plotted on ",
                    "top of the violin/box plots."),
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
               bsTooltip(id = "group1",
                         title =
                           paste0("Determines how to group samples into ",
                                  "the subplots. If set to \"annotation\", ",
                                  "then a sample annotation must be supplied ",
                                  "via the annotation parameter."),
                         placement = "right",
                         options = list(container = "body"))
             )
      )
      removeUI(selector = "#color")
      insertUI(
        selector = "#insert_color",
        ui = tagList(
          radioButtons(
            inputId = "color",
            label = "Color By",
            choices = list("Signature" = "signature",
                           "Annotation" = "annotation"),
            inline = TRUE,
            selected = "signature"
          ),
          bsTooltip(id = "color", title = "Determines how to color samples.",
                    placement = "right", options = list(container = "body"))
        )
      )
      vals$point_ind <- 1
    }
    #if scatter is selected, generate point size option, remove group by option,
    #color by has 3 choices
    else if(input$plottype == "scatter"){
      removeUI(selector = "#point_opt")
      removeUI(selector = "#group1")
      insertUI(
        selector = "#points",
        ui = tags$div(
          id = "point_opt",
          numericInput(inputId = "pointsize", label = "Point Size", value = 0.7),
          bsTooltip(id = "pointsize", title = "Size of the points on scatter plots.",
                    placement = "right", options = list(container = "body"))
        )
      )
      removeUI(selector = "#color")
      insertUI(
        selector = "#insert_color",
        ui = tagList(
          radioButtons(
            inputId = "color",
            label = "Color By",
            choices = list("None" = "none","Signature" = "signatures",
                           "Annotation" = "annotation"),
            inline = TRUE,
            selected = "none"
          ),
          bsTooltip(id = "color", title = "Determines how to color samples.",
                    placement = "right", options = list(container = "body"))
        )
      )
      vals$point_ind <- 0
    }
    #if bar is selected, remove all points related options, group by has 3 choices,
    #color by has 2 choices
    else if(input$plottype == "bar"){
      removeUI(selector = "#point_opt")
      removeUI(selector = "#group1")
      insertUI(
        selector = "#insert_group",
        ui = tagList(
          radioButtons(
            inputId = "group1",
            label = "Group By",
            choices = list("None" = "none", "Signature" = "signature",
                           "Annotation" = "annotation"),
            inline = TRUE,
            selected = "none"
          ),
          bsTooltip(id = "group1",
                    title =
                      paste0("Determines how to group samples into the ",
                             "subplots. If set to \"annotation\", then a ",
                             "sample annotation must be supplied via ",
                             "the annotation parameter."),
                    placement = "right", options = list(container = "body"))
        )
      )
      removeUI(selector = "#color")
      insertUI(
        selector = "#insert_color",
        ui = tagList(
          radioButtons(
            inputId = "color",
            label = "Color By",
            choices = list("Signature" = "signature",
                           "Annotation" = "annotation"),
            inline = TRUE,
            selected = "signature"
          ),
          bsTooltip(id = "color", title = "Determines how to color samples.",
                    placement = "right", options = list(container = "body"))
        )
      )
      vals$point_ind <- 0
    }
    else{
      return(NULL)
    }
  })

  addTooltip(session, id = "proportional",
  title = "If checked, the exposures will be normalized to
  between 0 and 1 by dividing by the total
  number of counts for each sample.",
             placement = "right", options = list(container = "body"))
  addTooltip(session, id = "color",
  title = "Determines how to color the bars or box/violins.
  If set to \"annotation\", then a sample annotation must be
  supplied via the annotation parameter",
             placement = "right", options = list(container = "body"))

  #if annotation is selected for group by, enable annotation option
  observeEvent(input$group1, {
    if (input$group1 == "annotation" & input$color != "annotation") {
      if (ncol(samp_annot(vals$result_objects[[input$selected_res2]])) == 1) {
        shinyalert::shinyalert(title = "Error",
                               text = paste0("Annotation not found. ",
                               "Please add annotation to the musica object."))
      }
      else{
        vals$annot <-
          as.list(colnames(samp_annot(
            vals$result_objects[[input$selected_res2]]))[-1])
        names(vals$annot) <-
          colnames(samp_annot(vals$result_objects[[input$selected_res2]]))[-1]
        insertUI(
          selector = "#insertannot",
          ui = tagList(
            selectInput(
              inputId = "annotation",
              label = "Annotation",
              choices = vals$annot
            ),
            bsTooltip(id = "annotation",
                      title = paste0("Sample annotation used to group the ",
                      "subplots or color the bars, boxes, or violins."),
                      placement = "right", options = list(container = "body"))
          )
        )
      }
    }
    else{
      if (input$color != "annotation") {
        removeUI(selector = "div:has(>> #annotation)")
      }
    }
  })

  #if annotation is selected for color by, enable annotation option
  observeEvent(input$color, {
    if (input$color == "annotation" & input$group1 != "annotation") {
      if (ncol(samp_annot(vals$result_objects[[input$selected_res2]])) == 1) {
        shinyalert::shinyalert(title = "Error",
                               text = "Annotation not found.
                               Please add annotation to the musica object.")
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
            bsTooltip(id = "annotation",
                      title = "Sample annotation used to group the subplots
                      or color the bars, boxes, or violins.",
                      placement = "right", options = list(container = "body"))
          )
        )
      }
    }
    else{
      if (input$group1 != "annotation") {
        removeUI(selector = "div:has(>> #annotation)")
      }
    }
  })

  #If bar plot is sorted by signature exposure, generate a bucket_list
  observeEvent(input$sort, {
    if (input$sort == "signature") {
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
              labels = as.list(colnames(
                vals$result_objects[[input$selected_res2]]@signatures)),
              input_id = "sig_from"
            ),
            add_rank_list(
              text = "Selected Signatures:",
              labels = NULL,
              input_id = "sig_to"
            )
          ),
          bsTooltip(id = "bucket",
          title = "Drag signatures from top bucket to the bottom bucket.
          Samples will be sorted in descending order by signatures.
          If multiple signatures are supplied,
          samples will be sorted by each signature sequentially",
          placement = "right", options = list(container = "body"))
        )
      )
    }
    else{
      removeUI(selector = "#insertsig")
    }
  })

  #generate the option to determine number of top samples to include in bar plot
  observeEvent(input$selected_res2, {
    output$number <- renderUI(
      tagList(
        numericInput(inputId = "numsamp", label = "# of Top Samples",
                     value = dim(vals$result_objects[[
                       input$selected_res2]]@exposures)[2],
                     min = 1,
                     max = dim(vals$result_objects[[
                       input$selected_res2]]@exposures)[2]),
        bsTooltip(id = "numsamp",
                  title = "The top number of sorted samples to display.",
                  placement = "right", options = list(container = "body"))
      )
    )
  })

  #save plotting options for exposures in a list
  get_exp_option <- function(input) {
    plot_type <- input$plottype
    proportional <- input$proportional
    group_by <- input$group1
    color_by <- input$color
    if (input$group1 == "annotation" | input$color == "annotation") {
      annot <- input$annotation
    }
    else{
      annot <- NULL
    }
    if (!is.numeric(input$numsamp)) {
      num_samples <- NULL
    }
    else{
      num_samples <- input$numsamp
    }
    sort_by <- input$sort
    if (sort_by == "signature") {
      sort_samples <- input$sig_to
    }
    else{
      sort_samples <- sort_by
    }
    if (!is.numeric(input$theta)) {
      threshold <- NULL
    }
    else{
      threshold <- input$theta
    }
    same_scale <- input$scale2
    label_x_axis <- input$xlab2
    legend <- input$legend2
    if (length(input$addpoint) == 0) {
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

  #plot exposures
  observeEvent(input$get_plot2, {
    options <- get_exp_option(input)
    if(length(umap(vals$result_objects[[input$selected_res2]])) == 0){
      create_umap(vals$result_objects[[input$selected_res2]])
      result <- vals$result_objects[[input$selected_res2]]
    }
    else{
      result <- vals$result_objects[[input$selected_res2]]
    }
    if(options[[14]]){
      jqui_resizable("#expplotly", operation = "destroy")
      jqui_resizable("#expplot", operation = "destroy")
      removeUI(selector = "#expplotly")
      removeUI(selector = "#expplot")
      insertUI(
        selector = "#plotdiv2",
        ui = plotlyOutput(outputId = "expplotly")
      )
      if(options[[1]] == "scatter"){
        output$expplotly <- renderPlotly(
          plot_umap(
            result = result,
            color_by = options[[4]],
            proportional = options[[2]],
            same_scale = options[[9]],
            annotation = options[[5]],
            plotly = options[[14]],
            legend = options[[11]]
          )
        )
      }
      else{
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
      jqui_resizable("#expplotly")
    }
    else{
      jqui_resizable("#expplotly", operation = "destroy")
      jqui_resizable("#expplot", operation = "destroy")
      removeUI(selector = "#expplot")
      removeUI(selector = "#expplotly")
      insertUI(
        selector = "#plotdiv2",
        ui = plotOutput(outputId = "expplot")
      )
      if(options[[1]] == "scatter"){
        output$expplot <- renderPlot(
          plot_umap(
            result = result,
            color_by = options[[4]],
            proportional = options[[2]],
            same_scale = options[[9]],
            annotation = options[[5]],
            plotly = options[[14]],
            legend = options[[11]]
          )
        )
      }
      else{
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
      jqui_resizable("#expplot")
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
  #Subsetting my signatures
  observeEvent(input$subset, {
    if (input$subset == "signature") {
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
              labels = as.list(colnames(
                vals$result_objects[[input$select_res_heatmap]]@signatures)),
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
    }
    else if (input$subset == "all_signatures") {
      removeUI(selector = "#insertsig")
      }
  })
  get_sigs <- function(input) {
    req(input$subset)
    if (input$subset == "signature") {
      sig <- input$sort_sigs
    }
    else if (input$subset == "all_signatures") {
      sig <- NULL
    }
    return(sig)
  }

  #Subsetting by tumors
  observeEvent(input$subset_tum, {
    if (input$subset_tum == "tumors") {
      insertUI(
        selector = "#sortbytum",
        ui = tags$div(
          id = "inserttum",
          checkboxGroupInput(
            "tum_val", "Available Samples:",
            as.list(unique(vals$result_objects[[
              input$select_res_heatmap
              ]]@musica@sample_annotations$Tumor_Subtypes)))
          )
        )
    }
    else{
      removeUI(selector = "#inserttum")
    }
  })
  #Subsetting by annotation
  observeEvent(input$subset_annot, {
    if (input$subset_annot == "annotation") {
      insertUI(
        selector = "#sortbyannot",
        ui = tags$div(
          id = "#insertannot",
          checkboxGroupInput("annot_val", "Available annotations:",
                             c(as.list(colnames(samp_annot(
                               vals$result_objects[[input$select_res_heatmap
                                                    ]]))))
        )
      ))
    }
    else{
      removeUI(selector = "#insertannot")
    }
  })
   #Add heatmap functionality
  observeEvent(input$get_heatmap, {
    sigs <- get_sigs(input)
    output$heatmap <- renderPlot({
      req(input$select_res_heatmap)
      input$get_heatmap
      isolate(plot_heatmap(res_annot =
                           vals$result_objects[[input$select_res_heatmap]],
                         proportional = propor(),
                         show_row_names = sel_row_names(),
                         show_column_names = sel_col_names(),
                         scale = zscale(), subset_signatures = c(sigs),
                         subset_tumor = input$tum_val,
                         annotation = input$annot_val,
                         column_title = paste0("Heatmap for ",
                                               input$select_res_heatmap)))
  })
 })
 ##############Clustering################
  
  #select box for clustering
  output$select_res3 <- renderUI({
    tagList(
      selectInput(
        inputId = "selected_res3",
        label = h3("Select Result"),
        choices = c(names(vals$result_objects))
      ),
      bsTooltip(id = "selected_res3",
                title = "Select one musica_result
                object to visualize signatures.",
                placement = "right", options = list(container = "body"))
    )
  })

  #generate options for selecting number of clusters
  observeEvent(input$selected_res3, {
    output$no_cluster1 <- renderUI(
      tagList(
        numericInput(inputId = "numclust1", label = "Max Number of Clusters",
                     value = 10,
                     min = 2,
                     max = dim(vals$result_objects[[
                       input$selected_res3]]@exposures)[2])
      )
    )
    output$no_cluster2 <- renderUI(
      tagList(
        numericInput(inputId = "numclust2", label = "Max Number of Clusters",
                     value = 2,
                     min = 2,
                     max = dim(vals$result_objects[[
                       input$selected_res3]]@exposures)[2])
      )
    )
  })

  #make plot for exploratory analysis
  observeEvent(input$explore, {
    jqui_resizable("#explore_plot", operation = "destroy")
    removeUI(selector = "#explore_plot")
    insertUI(
      selector = "#insert_explore_plot",
      ui = plotOutput(outputId = "explore_plot")
    )
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
    jqui_resizable("#explore_plot")
  })

  #generate select box for dissimilarity matrix and options specific to
  #certain algorithms
  observeEvent(input$algorithm2, {
    choices <- list(hkmeans =
                      c("Euclidean" = "euclidean", "Manhattan" = "manhattan",
                        "Canberra" = "canberra"),
                    clara = c("Euclidean" = "euclidean",
                              "Manhattan" = "manhattan", "Jaccard" = "jaccard"),
                    kmeans = c("Euclidean" = "euclidean",
                               "Manhattan" = "manhattan", "Jaccard" = "jaccard",
                               "Cosine" = "cosine", "Canberra" = "canberra"),
                    hclust = c("Euclidean" = "euclidean",
                               "Manhattan" = "manhattan", "Jaccard" = "jaccard",
                               "Cosine" = "cosine", "Canberra" = "canberra"),
                    pam = c("Euclidean" = "euclidean",
                            "Manhattan" = "manhattan", "Jaccard" = "jaccard",
                            "Cosine" = "cosine", "Canberra" = "canberra"))
    output$diss <- renderUI(
      selectInput(
        inputId = "diss_method",
        label = "Method for Dissimilarity Matrix",
        choices = choices[[input$algorithm2]]
      )
    )
    if (input$algorithm2 == "hclust") {
      insertUI(
        selector = "#hclust",
        ui = selectInput(
          inputId = "hclust_method",
          label = "Hierarchical Clustering Method",
          choices = c("ward.D" = "ward.D", "ward.D2" = "ward.D2",
                      "single" = "single", "complete" = "complete",
                      "average" = "average", "mcquitty" = "mcquitty",
                      "median" = "median", "centroid" = "centroid")
        )
      )
    }
    else{
      removeUI(selector = "div:has(>> #hclust_method)")
    }
    if (input$algorithm2 == "clara") {
      insertUI(
        selector = "#clara",
        ui = numericInput(inputId = "clara_num",
                          label = "No. of Samples for CLARA",
                          value = 5,
                          min = 1,
                          max = dim(vals$result_objects[[
                            input$selected_res3]]@exposures)[2])
      )
    }
    else{
      removeUI(selector = "div:has(>> #clara_num)")
    }
    if (input$algorithm2 %in% c("kmeans", "hkmeans")) {
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

  #create select box to allow users to color plot by annotation
  observeEvent(input$group2, {
    if (input$group2 == "annotation") {
      vals$annot <- as.list(colnames(samp_annot(
        vals$result_objects[[input$selected_res3]]))[-1])
      names(vals$annot) <- colnames(samp_annot(
        vals$result_objects[[input$selected_res3]]))[-1]
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

  #perform clustering analysis
  observeEvent(input$cluster_calc, {
    result <- vals$result_objects[[input$selected_res3]]
    nclust <- input$numclust2
    proportional <- input$proportional3
    method <- input$algorithm2
    dis.method <- input$diss_method
    if (!is.null(input$hclust_method)) {
      hc.method <- input$hclust_method
    }
    else{
      hc.method <- "ward.D"
    }
    if (!is.numeric(input$clara_num)) {
      clara.samples <- 5
    }
    else{
      clara.samples <- input$clara_num
    }
    if (!is.numeric(input$max_iter)) {
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
      content = function(file) {
        write.table(dat[-1], file, sep = "\t", quote = F)
      }
    )
  })

  #visualize clustering result
  observeEvent(input$cluster_vis, {
    if (length(umap(vals$result_objects[[input$selected_res3]])) == 0) {
      create_umap(vals$result_objects[[input$selected_res3]])
      result <- vals$result_objects[[input$selected_res3]]
    }
    else{
      result <- vals$result_objects[[input$selected_res3]]
    }
    clusters <- vals$cluster
    group <- input$group2
    if (group == "annotation") {
      annotation <- input$annotation2
    }
    else{
      annotation <- NULL
    }
    plotly <- input$plotly3
    if (plotly) {
      jqui_resizable("#cluster_plot", operation = "destroy")
      jqui_resizable("#cluster_plotly", operation = "destroy")
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
      jqui_resizable("#cluster_plotly")
    }
    else{
      jqui_resizable("#cluster_plot", operation = "destroy")
      jqui_resizable("#cluster_plotly", operation = "destroy")
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
      jqui_resizable("#cluster_plot")
    }
  })
  ########################################
  ########################################Download#############################
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
  #Adding the download feature
  output$download_mus_obj <- downloadHandler(

    filename = function() {
      paste("musica_object", ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(vals$musica, file = file)
    }
  )
 output$download_res <- downloadHandler(

    filename = function() {
      paste("musica_results", ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(vals$result_objects[[input$select_res_download]], file = file)
    }
  )
}
