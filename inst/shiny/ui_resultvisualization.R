library(plotly)
library(shinyBS)
musicaresultvisualization <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),
  h2("Visualize Signatures & Exposures"),
  #actionButton(inputId = "get_res", label = "Get Results"),
  fluidRow(
    tabBox(
      title = "",
      id = "tabset1",
      width = NULL,
      tabPanel(style = "position: relative;",
        title = "Signatures",
        h3("Settings"),
        uiOutput(outputId = "select_res1"),
        uiOutput(outputId = "select_modality1"),
        uiOutput(outputId = "select_model1"),
        checkboxInput(inputId = "rename", label = "Rename Signatures", value = FALSE, width = "100%"),
        tags$div(id = "signame"),
        numericInput(inputId = "text_size1", label = "Text Size", value = 10, width = "50%"),
        #numericInput(inputId = "facet_size", label = "Facet Size", value = 10, width = "50%"),
        #checkboxInput(inputId = "legend1", label = "Legend", value = TRUE),
        checkboxInput(inputId = "xlab1", label = "X-axis Label", value = TRUE),
        checkboxInput(inputId = "ylab1", label = "Y-axis Label", value = TRUE),
        checkboxInput(inputId = "scale1", label = "Same Y-axis Scale", value = TRUE),
        checkboxInput(inputId = "percent1", label = "Percent", value = TRUE),
        checkboxInput(inputId = "plotly1", label = "Plotly", value = TRUE),
        #checkboxInput(inputId = "y_max1", label = "Manual Y-axis Maximums", value = FALSE, width = "100%"),
        #tags$div(id = "y_maxs"),
        #checkboxInput(inputId = "annotation1", label = "Add Signature Annotations", value = FALSE, width = "100%"),
        #tags$div(id = "annotations"),
        actionButton(inputId = "get_plot1", label = "Make Plot"),
        tags$a(href = "#", 
               tags$i(class = "fas fa-question-circle"),
               title = "Need help?", 
               `data-toggle` = "popover", 
               `data-trigger` = "focus", 
               `data-content` = "Bar plots can be used to display the probability of each type of mutation within each signature.
               By default, signatures are named by numbers, but an option of renaming signatures is provided if you want to rename 
               them according to their possible etiology (e.g. Smoking) or closest correlation to COSMIC (e.g. SBS4).",
               `data-html` = "true",
               `data-placement` = "left",
               style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
        tags$div(id = "plot_div1"),
        bsTooltip(id = "rename", title = "If checked, the names of signatures can be customized.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "text_size1", title = "Size of axis text.", placement = "right", options = list(container = "body")),
        #bsTooltip(id = "facet_size", title = "Size of facet text.", placement = "right", options = list(container = "body")),
        #bsTooltip(id = "legend1", title = "If checked, the legend for mutation types will be included in the plot.",
                  #placement = "right", options = list(container = "body")),
        bsTooltip(id = "xlab1", title = "If checked, the labels for the mutation types on the x-axis will be shown.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "ylab1", title = "If checked, the labels and tick marks on the y-axis will be shown.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "scale1", title = "If checked, the scale of the probability for each signature will be the same.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "plotly1", title = "If checked, the plot will be made interactive using plotly.",
                  placement = "right", options = list(container = "body"))
        #bsTooltip(id = "ymax1", title = "If supplied, the maximum value on the y-axis of each signature will be adjusted.", 
                  #placement = "right", options = list(container = "body")),
        #bsTooltip(id = "annotation1", title = "If supplied, annotation text will be displayed in the top right corner of each signature.", 
                  #placement = "right", options = list(container = "body"))
      ),
      tabPanel(title = "Exposures",
               fluidRow(
                 column(
                   width = 6,
                   box(
                     width = 12,
                     h3("General Settings"),
                     uiOutput(outputId = "select_res2"),
                     uiOutput(outputId = "select_modality2"),
                     uiOutput(outputId = "select_model2"),
                     radioButtons(
                       inputId = "plot_type",
                       label = "Plot Type",
                       choices = list("Bar Plot" = "bar", "Box Plot" = "box", 
                                      "Violin Plot" = "violin", "Scatter Plot" = "scatter"),
                       inline = TRUE,
                       selected = "bar"
                     ),
                     checkboxInput(inputId = "proportional", label = "Proportional", value = TRUE),
                     tags$div(id = "insert_group"),
                     radioButtons(
                       inputId = "group1",
                       label = "Group By",
                       choices = list("Signature" = "signature",
                                      "Annotation" = "annotation"),
                       inline = TRUE,
                       selected = "signature"
                     ),
                     tags$div(id = "insert_color"),
                     radioButtons(
                       inputId = "color",
                       label = "Color By",
                       choices = list("Signature" = "signature", "Annotation" = "annotation"),
                       inline = TRUE,
                       selected = "signature"
                     ),
                     tags$div(id = "insert_annot"),
                     tags$a(href = "#", 
                            tags$i(class = "fas fa-question-circle"),
                            title = "Need help?", 
                            `data-toggle` = "popover", 
                            `data-trigger` = "focus", 
                            `data-content` = "Bar plots, box plots, violin plots and scatter plots can be used to visualize the exposure levels of each signature in each sample.",
                            `data-html` = "true",
                            `data-placement` = "left",
                            style = "position: absolute; top: 5px; right: 5px; cursor: pointer;")
                   )
                 )
               ),
               fluidRow(
                 column(
                   width = 6,
                   box(
                     width = 12,
                     h3("Sorting"),
                     radioButtons(
                       inputId = "sort",
                       label = "Sort By",
                       choices = list("Total Counts" = "total", "Sample Name" = "name",
                                      "Signatures" = "signature"),
                       inline = TRUE,
                       selected = "total"
                     ),
                     tags$div(id = "sort_by_sig"),
                     uiOutput(outputId = "number"),
                     numericInput(inputId = "theta", label = "Threshold", value = NULL),
                     bsTooltip(id = "sort", title = "Used to sort bar plot from left to right.", 
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "group1", title = "Determines how to group samples into the subplots. 
                               If set to \"annotation\", then a sample annotation must be supplied via the annotation parameter.",
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "theta", title = "Exposures less than this threshold will be set to 0.",
                               placement = "right", options = list(container = "body"))
                   )
                 ),
                 column(
                   width = 6,
                   box(
                     width = 12,
                     h3("Aesthetic Settings"),
                     checkboxInput(inputId = "scale2", label = "Same Y-axis Scale", value = TRUE),
                     checkboxInput(inputId = "xlab2", label = "X-axis Label", value = FALSE),
                     checkboxInput(inputId = "legend2", label = "Legend", value = TRUE),
                     tags$div(id = "points"),
                     checkboxInput(inputId = "plotly2", label = "Plotly", value = TRUE),
                     bsTooltip(id = "scale2", title = "If checked, then all subplots will have the same scale.",
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "xlab2", title = "If checked, x-axis labels will be displayed at the bottom of the plot.",
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "legend2", title = "If checked, the legend will be displayed.",
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "plotly2", title = "If checked, the plot will be made interactive using plotly.",
                               placement = "right", options = list(container = "body"))
                   )
                 )
               ),
               actionButton(inputId = "get_plot2", label = "Make Plot"),
               tags$div(id = "plot_div2")
      )
    )
  )
)
