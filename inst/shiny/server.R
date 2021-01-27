library(musicatk)

source("Tables.R", local = T)

server <- function(input, output, session) {
  observeEvent(input$get_musica, {
    maf <-  GDCquery_Maf("BRCA", pipelines = "mutect")
  })

  data(res_annot)
    
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
 
##################Visualization#################   
  observeEvent(input$rename,{
    n <- ncol(res_annot@signatures)
    for(i in 1:n){
      id <- paste0("sig", i)
      if(input$rename == TRUE){
        insertUI(
          selector = '#signame',
          ui = textInput(inputId = id, paste0("Signature", i))
        )
      }
      else{
        removeUI(selector = paste0("div:has(> #", id, ")"))
      }
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$get_plot1,{
    if(input$rename){
      ids <- vector()
      for(i in 1:ncol(res_annot@signatures)){
        ids <- c(ids, input[[paste0("sig", i)]])
      }
      name_signatures(result = res_annot, ids)
    }
    output$sigplot <- renderPlot({
      plot_signatures(result = res_annot)
    })
  })
  observeEvent(input$get_plot2,{
    prop <- input$proportional
    output$expplot <- renderPlot(
      plot_exposures(result = res_annot, plot_type = "bar", proportional = prop)
    )
  })
################################################
}
