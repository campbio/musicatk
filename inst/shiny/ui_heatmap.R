shinyPanelheatmap <- fluidPage(
  fluidRow(column(width = 12,
                  uiOutput(outputId = "select_res_heatmap"),
                  h3("Settings"),
                  checkboxInput("prop", "Proportional", FALSE),
                  checkboxInput("col_names", "Show column names", FALSE),
                  checkboxInput("row_names", "Show row names", TRUE),
                  checkboxInput("scale", "Z-score normalization", TRUE),
                  hr(),
                  radioButtons(
                    inputId = "subset",
                    label = "Subset By",
                    choices = list("Signature" = "signature"),
                    inline = TRUE,
                    selected = ""
                  ),
                  tags$div(id = "sortbysigs"),
                  radioButtons(
                    inputId = "subset_tum",
                    label = "",
                    choices = list("Tumors" = "tumors"),
                    inline = TRUE,
                    selected = ""
                  ),
                  tags$div(id = "sortbytum"),
                  radioButtons(
                    inputId = "subset_annot",
                    label = "Annotate by",
                    choices = list("Annotation" = "annotation"),
                    inline = TRUE,
                    selected = ""
                  ),
                  tags$div(id = "sortbyannot"),
                  hr(),
                  actionButton("get_heatmap", "Plot"),
                  downloadButton("download_heatmap", "Download"),
                  hr(),
                  plotOutput("heatmap")
                  
)))