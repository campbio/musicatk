shinyPanelImport <- fluidPage(
          fluidRow(
            useShinyjs(),
            column(width = 12,
                   # h3("Step 1: Select File")
            )
          ),
          hr(),
          fluidRow(box(width = 12,
            useShinyalert(),
            add_busy_spinner(spin = "fading-circle"),
            fileInput("file", "Select file:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          actionButton("upload", "Add samples"))),
          
          fluidRow(box(width = 12,h3("Added files"),
          #tags$div(id = "file_id",dataTableOutput("my_file_name"),style =  "font-size:40%"),
          tags$div(id = "file_id",DT::dataTableOutput("dtable"),style = "font-size:40%; overflow-y: scroll;overflow-x: scroll;"),
          uiOutput('undoUI'))),
          hr(),
          actionButton("import", "Import"),
          uiOutput("spinner"),
          fluidRow(box(width = 12,div(dataTableOutput("musica_contents")))),
          downloadButton("download_musica", "Download Variants"),
          bsTooltip("upload", "Press button to add your uploaded files to Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("file_id", "Table of files that have been added by you", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("import", "Press button to import the files in the Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("undoUI", "Press button to undo any deletes you made in the Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          
          
)
          