shiny_panel_annotations <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),   
  div(style = "position: relative;",
      box(width = 6,
    # Get list of musica and result objects
    uiOutput("annotation_musica_list"),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "Sample annotations can be used to store information about each sample such as tumor type or treatment status. 
           They are optional and can be used in downstream plotting functions such as plot_exposures or plot_umap to group or color samples by a particular annotation.
           Annotations are sample-level variables such as Tumor Type or Age. These can be provided in a character delimited text file which must include a column containing the sample IDs.
           The sample IDs must match the sample names in the musica object. NA values will be given to any samples not present in the annotation file.
           After selecting the correct delimiter, a data table will appear below to show the annotations that will be added to the musica object. If your annotation file does not contain a header, unselect the “Header” radio button. Choose the column that contains the sample names from the “Sample Name Columns” dropdown and then click “Add Annotation",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
    wellPanel(id = "musica_data",
            h3("Upload Sample Annotations"),
            fileInput("annotations_file", "Annotations file:",
                      multiple = TRUE),
            checkboxInput("annotation_header", "Header", TRUE),
            radioButtons("annotation_delimiter", "Delimiter",
                         choices = list("comma" = ",",
                                        "tab" = "\t",
                                        "space" = " ",
                                        "pipe" = "|",
                                        "semicolon" = ";",
                                        "custom" = "custom"),
                         selected = ",", inline = TRUE),
            hidden(textInput("CustomAnnotDelim", "Delimiter")),
            uiOutput("annotation_samples")
  ),
  actionButton("add_annotation", "Add Annotation"),
  bsTooltip("add_annotation",
            "Add annotations to your existing data for downstream analysis.",
            placement = "bottom", trigger = "hover", options = NULL),
  bsTooltip("annotation_delimiter",
            "Choose the delimiter to parse your text file.",
            placement = "right", trigger = "hover", options = NULL))),
  box(width = 12,
    DT::DTOutput("annotations")
  )
)
