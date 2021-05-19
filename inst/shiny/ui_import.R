shinyPanelImport <- fluidPage(
          fluidRow(
            
                useShinyjs(),
                   # h3("Step 1: Select File")
            )
          ,
          hr(),
          fluidRow(box(width = 6,
            useShinyalert(),
            add_busy_spinner(spin = "fading-circle"),
            tags$style("#file {
                    font-size:3px;
                    height:10px;
           }"),
            fileInput("file", "Select file:",
                  multiple = TRUE,
                  accept = c(".maf",".vcf")),
          actionButton("upload", "Add samples"))),
          
          fluidRow(box(width = 6,h3("Added files"),
          #tags$div(id = "file_id",dataTableOutput("my_file_name"),style =  "font-size:40%"),
          tags$div(id = "file_id",DT::dataTableOutput("dtable"),style = "font-size:40%; overflow-y: scroll;overflow-x: scroll;"),
          uiOutput('undoUI')),
          hr(),
          
          uiOutput("spinner")),
          
          fluidRow(box(width = 6,actionButton("import", "Import"),downloadButton("download_musica", "Download Variants"),div(dataTableOutput("musica_contents"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;") 
                       )),
          bsTooltip("upload", "Press button to add your uploaded files to Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("file_id", "Table of files that have been added by you", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("import", "Press button to import the files in the Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("undoUI", "Press button to undo any deletes you made in the Sample List", placement = "bottom", trigger = "hover",
                    options = NULL),
          
          
)
          
