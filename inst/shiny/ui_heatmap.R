shinyPanelheatmap <- fluidPage(
  fluidRow(box(uiOutput(outputId = "select_res_heatmap"),
               #uiOutput(outputId = "sig_ui"),
                  h3("Settings"),
                  checkboxInput("prop", "Proportional", FALSE),
                  checkboxInput("col_names", "Show column names", FALSE),
                  checkboxInput("row_names", "Show row names", TRUE),
                  checkboxInput("scale", "Z-score normalization", TRUE))),
                
           
           fluidRow(box(
                  radioButtons(
                    inputId = "def_sigs",
                    label = "Subset by",
                    choices = list("All Signatures" = "default_signatures"),
                    inline = TRUE,
                    selected = ""
                  ),
                  radioButtons(
                    inputId = "subset",
                    label = "",
                    choices = list("Selected Signatures" = "signature"),
                    inline = TRUE,
                    selected = ""
                  ),
                  tags$div(id = "sortbysigs"),
                  
                  
                  radioButtons(
                    inputId = "subset_tum",
                    label = "",
                    choices = list("Samples" = "tumors"),
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
                  tags$div(id = "sortbyannot")
                  )),
                  hr(),
                  actionButton("get_heatmap", "Plot"),
                  downloadButton("download_heatmap", "Download"),
                  hr(),
                  plotOutput("heatmap"),
                  bsTooltip("select_res_heatmap", "Select musica result object for plotting a heatmap", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("prop", "Check box to normalize exposures", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("col_names", "Check box to show column names", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("row_names", "Check box to show row names", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("scale", "Check box to normalize by the z-score", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset", "Choose for subsetting data by signatures present", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset_tum", "Choose for subsetting by available samples ", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("subset_annot", "Choose for subsetting by available annotations", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("get_heatmap", "Press button to plot heatmap", placement = "bottom", trigger = "hover",
                  options = NULL),
                  bsTooltip("download_heatmap", "Press button to download plot", placement = "bottom", trigger = "hover",
                  options = NULL),
)