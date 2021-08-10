shiny_panel_import <- fluidPage(
   fluidRow(box(width = 6, #setting box width
          fluidRow(box(width = 6,
              div(fileInput("file", "Select file:",
                  multiple = TRUE,
                  accept = c(".maf", ".vcf")),
          style = "font-size:50%;"),
          actionButton("upload", "Add samples"))),
          fluidRow(box(width = 6, h3("Added files"),
          tags$div(id = "file_id", DT::dataTableOutput("dtable"),
                   style = "font-size:40%;
                   overflow-y: scroll;overflow-x: scroll;"),
          uiOutput("undo_ui"),
          hr(),
          uiOutput("spinner"))),
          #Adding the download button
          fluidRow(box(width = 6, actionButton("import", "Import"),
                       downloadButton("download_musica", "Download Variants"),
                       hr(),
                       div(dataTableOutput("musica_contents"),
                           style = "font-size:40%; height:500px;
                           overflow-y: scroll;overflow-x: scroll;"),
                       )),
          #Adding help tootltips
          bsTooltip("upload",
                    "Press button to add your uploaded files to Sample List",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("file_id", "Table of files that have been added by you",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("import",
                    "Press button to import the files in the Sample List",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
          bsTooltip("undo_ui",
                    "Press button to undo any deletions
                    you made in the Sample List",
                    placement = "bottom", trigger = "hover",
                    options = NULL),
)))
