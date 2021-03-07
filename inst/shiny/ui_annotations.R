shinyPanelAnnotations <- fluidPage(
  box( width = 12,
    # column(width = 12,
    #        h3("Upload Sample Annotations"),
    #        actionButton("AddAnnotation", "Import an annotation text file.")
    # ),
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
  dataTableOutput("annotations")

  )
)