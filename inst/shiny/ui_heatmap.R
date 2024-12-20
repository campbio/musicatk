shiny_panel_heatmap <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),  
  #Adding box formating
  fluidRow(div(style = "position: relative;",box(
    uiOutput(outputId = "select_res_heatmap"),
    uiOutput(outputId = "select_modality_heatmap"),
    uiOutput(outputId = "select_model_heatmap"),
                  h3("Settings"),
                  #Adding checkboxees for different heatmap fucntion parameters
                  checkboxInput("prop", "Proportional", FALSE),
                  checkboxInput("col_names", "Show column names", FALSE),
                  checkboxInput("row_names", "Show row names", TRUE),
                  checkboxInput("scale", "Z-score normalization", TRUE)))),
                  tags$a(href = "#", 
                         tags$i(class = "fas fa-question-circle"),
                         title = "Need help?", 
                         `data-toggle` = "popover", 
                         `data-trigger` = "focus", 
                         `data-content` = "Select the result object you want to use for heatmap visualization from the 
                         dropdown menu in the start labeled ‘Select Result’. Choose from different settings and press 
                         the ‘Plot’ button to visualize the heatmap. Check the ‘Proportional’ option to normalize the 
                         levels of signatures within a sample to sum to one. Z-score normalization will scale each 
                         signature to have a mean of zero and a standard deviation of one.. You can also show or hide 
                         the tumor and signature names by checking or unchecking the ‘Show column names’ and ‘Show row
                         names’, respectively. If you want to display a subset by signatures, click the ‘Selected 
                         Signatures’ option and select the signatures to show. You can also select the ‘Samples’ and 
                         ‘Annotation’ options to subset by the available samples and annotations.",
                         `data-html` = "true",
                         `data-placement` = "left",
                         style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),
                 fluidRow(box(
                  radioButtons(
                    inputId = "subset",
                    label = "Subset Signatures",
                    choices = list("All Signatures" = "all_signatures",
                                   "Selected Signatures" = "signature"),
                    inline = TRUE,
                    selected = "all_signatures"
                  ),
                  tags$div(id = "sortbysigs"),
                  radioButtons(
                    inputId = "subset_tum",
                    label = "Subset Samples",
                    choices = list("All Samples" = "all_samples",
                                   "Selected Tumor Types" = "tumors"),
                    inline = TRUE,
                    selected = "all_samples"
                  ),
                  tags$div(id = "sortbytum"),
                  h3("Annotation"),
                  checkboxInput("subset_annot", "Add annotation", FALSE),
                  #radioButtons(
                  #  inputId = "subset_annot",
                  #  label = "Annotate by",
                  #  choices = list("Annotation" = "annotation"),
                  #  inline = TRUE,
                  #  selected = ""
                  #),
                  tags$div(id = "sortbyannot"),
                  actionButton("get_heatmap", "Plot"),
                  )),
                  plotOutput("heatmap"),
                  #Adding help tooltips
                  bsTooltip("select_res_heatmap",
                            "Select musica result object for plotting a heatmap",
                  placement = "bottom", trigger = "hover", options = NULL),
                  bsTooltip("prop",
                            "Check box to normalize exposures",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("col_names",
                            "Check box to show column names",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("row_names"
                            , "Check box to show row names",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("scale",
                            "Check box to normalize by the z-score",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset",
                            "Choose for subsetting data by signatures present",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset_tum",
                            "Choose for subsetting by available samples ",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset_annot",
                            "Choose for subsetting by available annotations",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("get_heatmap", 
                            "Press button to plot heatmap",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("download_heatmap",
                            "Press button to download plot",
                  placement = "bottom", trigger = "hover",
                  options = NULL),
)
