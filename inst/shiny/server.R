library(musicatk)
library(plotly)
library(sortable)

source("Tables.R", local = T)

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
