shiny_panel_annotations <- fluidPage(
  box(width = 6,
    # Get list of musica and result objects
    uiOutput("annotation_musica_list"),
    wellPanel(id = "musica_data",
            h3("Upload Sample Annotations"),
            fileInput("annotations_file", "Annotations file:",
                      multiple = TRUE),
            checkboxInput("annotation_header", "Header", T),
            radioButtons("annotation_delimiter", "Delimiter",
                         choices = list("comma" = ",",
                                        "tab" = "\t",
                                        "space" = " ",
                                        "pipe" = "|",
                                        "semicolon" = ";",
                                        "custom" = "custom"),
                         selected = ",", inline = T),
            hidden(textInput("CustomAnnotDelim", "Delimiter")),
            uiOutput("annotation_samples")
  ),
  actionButton("add_annotation", "Add Annotation"),
  bsTooltip("add_annotation",
            "Add annotations to your existing data for downstream analysis.",
            placement = "bottom", trigger = "hover", options = NULL),
  bsTooltip("annotation_delimiter",
            "Choose the delimiter to parse your text file.",
            placement = "right", trigger = "hover", options = NULL)),
  box(width = 12,
    dataTableOutput("annotations")
  )
)
