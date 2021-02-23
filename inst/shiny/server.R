library(musicatk)
options(shiny.maxRequestSize = 100*1024^2)
source("Tables.R", local = T)
vals <- reactiveValues(var = NULL,df = data.frame(),musica_contents = NULL,musica_uplaod = data.frame())

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
 
 

output$my_file_name <- renderDataTable({
 
   vals$df <- cbind(input$file[['name']],data.frame(delete = glue(
     "<button  rowid='{1:length(input$file[['name']])}' 
        onclick='Shiny.setInputValue(\"removeRow\",this.getAttribute(\"rowid\"))'>Delete</button>"),
     rowid = 1:length(input$file[['name']])
   ))
   
   datatable(
     vals$df,
     escape = FALSE,
     selection = "none",
     options = list(
     columnDefs = list(list(targets = ncol(vals$df),visible = FALSE))
     ))
   #return(vals$df)
  })

observeEvent(input$removeRow,{
  removeRow <- as.integer(input$removeRow)
  tblRowRemoved <- which(vals$df$rowid == removeRow)
  vals$df <- vals$df[-tblRowRemoved,]
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
  observeEvent(input$overwriteTable, { 
    overwrite <- T
  })
  observeEvent(input$keepTable, { 
    overwrite <- F 
  })
  observeEvent(input$AddTable, {
    add_tables(input)
  })

  # Test when musica code has been generated
  # observeEvent(input$MusicaResults, {
  #   musica_result <- discover_signatures(
  #     output$musica, table_name = input$SelectTable, 
  #     num_signatures = input$NumberOfSignatures,
  #     method = input$Methods,
  #     seed = input$Seed,
  #     nstart = input$nStart)
  # })
###############################################################################
  
}
