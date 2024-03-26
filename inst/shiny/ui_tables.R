shiny_panel_tables <- fluidPage(
  fluidRow(            tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://use.fontawesome.com/releases/v5.15.4/css/all.css"),
    tags$script(HTML('
            $(document).ready(function(){
             $("[data-toggle=\'popover\']").popover();
            });
             '))
  ),
  div(style = "position: relative;",
      box(width = 6,
    selectInput("select_table", "Select Count Table",
                choices = list("SBS96", "SBS192 - Transcript_Strand",
                "SBS192 - Replication_Strand", "DBS", "Indel"),
                selected = 1),
    selectInput("table_genome_list", "Reference genome:",
                list("hg38", "hg19", "mm9", "mm10"),
                width = "100%"),
    textOutput("table_genome_warning"),
    uiOutput("allow_table"),
    shinybusy::use_busy_spinner(spin = "double-bounce"),
    tags$a(href = "#", 
           tags$i(class = "fas fa-question-circle"),
           title = "Need help?", 
           `data-toggle` = "popover", 
           `data-trigger` = "focus", 
           `data-content` = "The “Build Tables” tab generates count tables for different mutation schemas. 
           These schemas are the input to the mutational signature discovery and prediction functions. To build the count tables, 
           the user must select 1 of the 5 standard motifs in the “Select Count Table” dropdown.",
           `data-html` = "true",
           `data-placement` = "left",
           style = "position: absolute; top: 5px; right: 5px; cursor: pointer;"),    
    bsTooltip("select_table",
              "Name of the standard table to build.",
              placement = "right", trigger = "hover", options = NULL)
  ))),
  fluidRow(
  uiOutput("combine_table")),
  shinybusy::use_busy_spinner(spin = "double-bounce")
)
