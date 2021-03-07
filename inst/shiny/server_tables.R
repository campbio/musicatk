add_tables <- function (input, vals) { 
  table_name <- input$SelectTable
  strand_type <- NULL
  if (input$SelectTable != "Custom") {
    # Check inputs for SBS192
    if (input$SelectTable == "SBS192 - Transcript_Strand") {
        annotate_transcript_strand(vals$musica, "19", build_table = F)
        table_name <- "SBS192"
        strand_type = "Transcript_Strand"
    }
    if (input$SelectTable == "SBS192 - Replication_Strand") {
        annotate_replication_strand(vals$musica, rep_range, build_table = F)
        table_name <- "SBS192"
        strand_type = "Replication_Strand"
    }
    tryCatch( {
      build_standard_table(vals$musica, vals$genome,
                         table_name = table_name,
                         strand_type = strand_type,
                         overwrite = T)
    },error = function(cond) {
      shinyalert::shinyalert(title = "Error", text = cond$message)
    }, warning = function(cond) {
      print(cond$message)
    }
    )
    return()
  }
  shinyalert::shinyalert(title = "Oops",
                         text = "Custom tables are not yet supported.")
}
