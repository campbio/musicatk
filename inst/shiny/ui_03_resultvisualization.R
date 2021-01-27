musicaresultvisualization <- fluidPage(
  h2("Visualize Signatures & Exposures"),
  fluidRow(
    tabBox(
      title = "",
      id = "tabset1",
      height = "400px",
      tabPanel(title = "Signatures",
               checkboxInput(inputId = "rename", label = "Rename Signatures", value = FALSE),
               tags$div(id = 'signame'),
               actionButton(inputId = "get_plot1", label = "Make Plot"),
               plotOutput(outputId = "sigplot")),
      tabPanel(title = "Exposures",
               checkboxInput(inputId = "proportional", label = "Proportional", value = TRUE),
               actionButton(inputId = "get_plot2", label = "Make Plot"),
               plotOutput(outputId = "expplot"))
    )
  )
)