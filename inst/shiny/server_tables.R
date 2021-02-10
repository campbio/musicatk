add_tables <- function (input, vals) { 
  strand_type <- input$StrandType
  # Check it table already exists
  print("hello")
  if (input$SelectTable != "Custom") {
    # Check inputs for SBS192
    if (input$SelectTable == "SBS192") {
      if (strand_type == "") {
        shinyalert::shinyalert(title = "Oops", 
                               text = "You must select strand type for table SBS192.")
      } else if (strand_type == "Transcript_Strand") {
        annotate_transcript_strand(vals$musica, "19")
        return ()
      } else {
        annotate_replication_strand(vals$musica, "19")
        return ()
      }
    }
    tryCatch( {
      build_standard_table(vals$musica, genome,
                         table_name = input$SelectTable, 
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
