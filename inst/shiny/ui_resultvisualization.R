library(plotly)
library(shinyBS)
musicaresultvisualization <- fluidPage(
  h2("Visualize Signatures & Exposures"),
  #actionButton(inputId = "get_res", label = "Get Results"),
  fluidRow(
    tabBox(
      title = "",
      id = "tabset1",
      width = NULL,
      tabPanel(
        title = "Signatures",
        h3("Settings"),
        uiOutput(outputId = "select_res1"),
        checkboxInput(inputId = "rename", label = "Rename Signatures", value = FALSE, width = "100%"),
        tags$div(id = "signame"),
        numericInput(inputId = "textsize1", label = "Text Size", value = 10, width = "50%"),
        numericInput(inputId = "facetsize", label = "Facet Size", value = 10, width = "50%"),
        checkboxInput(inputId = "legend1", label = "Legend", value = TRUE),
        checkboxInput(inputId = "xlab1", label = "X-axis Label", value = TRUE),
        checkboxInput(inputId = "scale1", label = "Same Y-axis Scale", value = TRUE),
        checkboxInput(inputId = "plotly1", label = "Plotly", value = TRUE),
        actionButton(inputId = "get_plot1", label = "Make Plot"),
        tags$div(id = "plotdiv1"),
        bsTooltip(id = "rename", title = "If checked, the names of signatures can be customized.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "textsize1", title = "Size of axis text.", placement = "right", options = list(container = "body")),
        bsTooltip(id = "facetsize", title = "Size of facet text.", placement = "right", options = list(container = "body")),
        bsTooltip(id = "legend1", title = "If checked, the legend for mutation types will be included in the plot.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "xlab1", title = "If checked, the labels for the mutation types on the x-axis will be shown.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "scale1", title = "If checked, the scale of the probability for each signature will be the same.",
                  placement = "right", options = list(container = "body")),
        bsTooltip(id = "plotly1", title = "If checked, the plot will be made interactive using plotly.",
                  placement = "right", options = list(container = "body"))
      ),
      tabPanel(title = "Exposures",
               fluidRow(
                 column(
                   width = 6,
                   box(
                     width = 12,
                     h3("General Settings"),
                     uiOutput(outputId = "select_res2"),
                     radioButtons(
                       inputId = "plottype",
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
                     tags$div(id = "insertannot")
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
                     tags$div(id = "sortbysig"),
                     uiOutput(outputId = "number"),
                     numericInput(inputId = "theta", label = "Threshold", value = NULL),
                     bsTooltip(id = "sort", title = "Used to sort bar plot from left to right.", 
                               placement = "right", options = list(container = "body")),
                     bsTooltip(id = "group1", title = "Determines how to group samples into the subplots. If set to \"annotation\", then a sample annotation must be supplied via the annotation parameter.",
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
               tags$div(id = "plotdiv2")
      )
    )
  )
)
