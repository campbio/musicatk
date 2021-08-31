# Function used in server.R to add counts tables.
add_tables <- function(input, vals) {
  table_name <- input$select_table
  strand_type <- NULL
  if (input$select_table != "Custom") {
    # Check inputs for SBS192
    if (input$select_table == "SBS192 - Transcript_Strand") {
        annotate_transcript_strand(vals$musica, input$table_genome_list,
                                   build_table = F)
        table_name <- "SBS192"
        strand_type <- "Transcript_Strand"
    }
    if (input$select_table == "SBS192 - Replication_Strand") {
        annotate_replication_strand(vals$musica, rep_range, build_table = F)
        table_name <- "SBS192"
        strand_type <- "Replication_Strand"
    }
    tryCatch({
      build_standard_table(vals$musica, select_genome(input$table_genome_list),
                         table_name = table_name,
                         strand_type = strand_type,
                         overwrite = T)
      shiny::showNotification("Table created.")
    }, error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
    }
    )
    return()
  }
  shinyalert::shinyalert(title = "Oops",
                         text = "Custom tables are not yet supported.")
}
