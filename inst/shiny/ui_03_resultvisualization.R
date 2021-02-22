musicaresultvisualization <- fluidPage(
  h2("Visualize Signatures & Exposures"),
  actionButton(inputId = "get_res", label = "Get Results"),
  fluidRow(
    tabBox(
      title = "",
      id = "tabset1",
      tabPanel(
        title = "Signatures",
        h3("Settings"),
        checkboxInput(inputId = "rename", label = "Rename Signatures", value = FALSE),
        tags$div(id = "signame"),
        numericInput(inputId = "textsize1", label = "Text Size", value = 10),
        numericInput(inputId = "facetsize", label = "Facet Size", value = 10),
        checkboxInput(inputId = "legend1", label = "Legend", value = TRUE),
        checkboxInput(inputId = "xlab1", label = "X-axis Label", value = TRUE),
        checkboxInput(inputId = "scale1", label = "Same Y-axis Scale", value = TRUE),
        checkboxInput(inputId = "plotly1", label = "Plotly", value = TRUE),
        actionButton(inputId = "get_plot1", label = "Make Plot"),
        tags$div(id = "plotdiv1")
      ),
      tabPanel(title = "Exposures",
               fluidRow(
                 box(
                   h3("General Settings"),
                   radioButtons(
                     inputId = "plottype",
                     label = "Plot Type",
                     choices = list("Bar Plot" = "bar", "Box Plot" = "box", 
                                    "Violin Plot" = "violin"),
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
                   radioButtons(
                     inputId = "color",
                     label = "Color By",
                     choices = list("Signature" = "signature", "Annotation" = "annotation"),
                     inline = TRUE,
                     selected = "signature"
                   ),
                   tags$div(id = "insertannot")
                 )
               ),
               fluidRow(
                 box(
                   h3("Sorting & Downsampling"),
                   radioButtons(
                     inputId = "sort",
                     label = "Sort By",
                     choices = list("Total Counts" = "total", "Sample Name" = "name",
                                    "Signatures" = "signature"),
                     inline = TRUE,
                     selected = "total"
                   ),
                   tags$div(id = "sortbysig"),
                   numericInput(inputId = "numsamp", label = "# of Top Samples", value = NULL),
                   numericInput(inputId = "theta", label = "Threshold", value = NULL)
                 )
               ),
               fluidRow(
                 box(
                   h3("Aesthetic Settings"),
                   checkboxInput(inputId = "scale2", label = "Same Y-axis Scale", value = TRUE),
                   checkboxInput(inputId = "xlab2", label = "X-axis Label", value = TRUE),
                   checkboxInput(inputId = "legend2", label = "Legend", value = TRUE),
                   tags$div(id = "points"),
                   checkboxInput(inputId = "plotly2", label = "Plotly", value = TRUE)
                 )
               ),
               actionButton(inputId = "get_plot2", label = "Make Plot"),
               tags$div(id = "plotdiv2")
      )
    )
  )
)
