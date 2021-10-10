shiny_panel_musica <- fluidPage(
 fluidRow(box(uiOutput("genome_list"),
    textOutput("genome_select"),
    column(width = 12,
           h3("Settings"),
           #Adding checkbox input ui
           checkboxInput("ref_chr", "Check Reference Chromosomes", TRUE),
           checkboxInput("ref_bases", "Check Reference Bases", TRUE),
           checkboxInput("convert_dbs", "Convert DBS", TRUE),
           checkboxInput("stand_indels", "Standardize Indels", TRUE),
           actionButton("get_musica_object", "Create Musica Object"),
           #Adding download button
           downloadButton("download_musica_object", "Download Musica Object"),
           hr(),
           div(textOutput("musica_contents_summary")),
           hr(),
           div(dataTableOutput("musica_contents_table"),
           style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
           hr(),
           #Adding help tooltips
           bsTooltip("genome_list",
                     "Full genome build versions of different organisms",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("ref_chr",
                     "Perform a check to ensure
                     that the chromosomes in the variant
                     object match the reference
                     chromosomes in the genome object.",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("ref_bases",
                     "Check if the reference bases in the variant object
                     match the reference bases in the genome object.",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("convert_dbs",
                     "Convert adjacent SBS into DBS (original SBS are removed)",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("stand_indels",
                     "Convert indel style (e.g. 'C > CAT' becomes
                     '- > AT' and 'GCACA > G' becomes 'CACA > -')",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("get_musica_object",
                     "Press button to create your musica object.",
                     placement = "bottom", trigger = "hover",
                     options = NULL),
           bsTooltip("download_musica",
                     "Press button to dwownload the musica variants.",
                     placement = "bottom", trigger = "hover",
                     options = NULL)
           )))
    )
