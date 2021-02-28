add_tables <- function (input, vals) { 
  if (input$SelectTable != "Custom") {
    # Check inputs for SBS192
    if (input$SelectTable == "SBS192 - Transcript_Strand") {
        annotate_transcript_strand(vals$musica, "19", build_table = T)
        return()
    }
    if (input$SelectTable == "SBS192 - Replication_Strand") {
        annotate_replication_strand(vals$musica, rep_range, build_table = T)
        return()
    }
    tryCatch( {
      build_standard_table(vals$musica, genome,
                         table_name = input$SelectTable, 
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
