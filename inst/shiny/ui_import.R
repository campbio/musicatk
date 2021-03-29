shinyPanelImport <- fluidPage(
          fluidRow(
            useShinyjs(),
            column(width = 12,
                    h3("Step 1: Select File")
            )
          ),
          
          fileInput("file", "File:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          actionButton("upload", "Add samples"),
          h3("Added samples"),
          #tags$div(id = "file_id",dataTableOutput("my_file_name"),style =  "font-size:40%"),
          tags$div(id = "file_id",DT::dataTableOutput("dtable"),style = "font-size:40%"),
          uiOutput('undoUI'),
          hr(),
          actionButton("import", "Import"),
          uiOutput("spinner"),
          div(dataTableOutput("musica_contents")),
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
          
