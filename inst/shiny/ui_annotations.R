shinyPanelAnnotations <- fluidPage(
  box( width = 12,
    # column(width = 12,
    #        h3("Upload Sample Annotations"),
    #        actionButton("AddAnnotation", "Import an annotation text file.")
    # ),
    helpText("Use this tab to add sample annotations for use in downstream analysis. This tab is optional. "),
    uiOutput("AnnotationMusicaList"),
    wellPanel(id = "musica_data",
            h3("Upload Sample Annotations"),
            fileInput("AnnotationsFile", "Annotations file:",
                      multiple = TRUE),
            checkboxInput("AnnotationHeader", "Header", T),
            radioButtons("AnnotationDelimiter", "Delimiter",
                         choices = list("comma" = ",",
                                        "tab" = "\t",
                                        "space" = " ",
                                        "pipe" = "|",
                                        "semicolon" = ";"),
                         selected = ",", inline = T)
  ),
  actionButton("AddAnnotation", "Add Annotation"),
  bsTooltip("AddAnnotation",
            "Add annotations to your existing data for downstream analysis.", 
            placement = "bottom", trigger = "hover", options = NULL),
  bsTooltip("AnnotationDelimiter",
            "Choose the delimiter to parse your text file.", 
            placement = "right", trigger = "hover", options = NULL), 
  dataTableOutput("annotations")

  )
)