shiny_panel_tables <- fluidPage(
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
    bsTooltip("select_table",
              "Name of the standard table to build.",
              placement = "right", trigger = "hover", options = NULL)
  ),
  uiOutput("combine_table"),
  shinybusy::use_busy_spinner(spin = "double-bounce")
)
