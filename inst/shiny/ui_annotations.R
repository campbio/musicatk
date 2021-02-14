shinyPanelAnnotations <- fluidPage(
  box(
    # column(width = 12,
    #        h3("Upload Sample Annotations"),
    #        actionButton("AddAnnotation", "Import an annotation text file.")
    # ),
  
    wellPanel(id = "musica_data",
            h3("Upload Sample Annotations"),
            fileInput("AnnotationsFile", "Annotations file:",
                      multiple = TRUE),
            checkboxInput("AnnotationHeader", "Header", T),
#           tableOutput("musica_contents"))
            uiOutput("annotations"),
            actionButton("AddAnnotation", "Add Annotation")
  )
  )
)