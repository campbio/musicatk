add_tables <- function (input) { 
  strand_type <- input$StrandType
  # Check it table already exists
  print("hello")
  if(input$SelectTable %in% names(extract_count_tables(musica)) &&
     input$OverwriteTable == F) {
    shinyalert::shinyalert(title = "Warning", 
                           text = "Table exists but was not overwritten. No changed were made.")
    return (musica)
  }
  if (input$SelectTable != "Custom") {
    # Check inputs for SBS192
    if (input$SelectTable == "SBS192") {
      if (strand_type == "") {
        shinyalert::shinyalert(title = "Oops", 
                               text = "You must select strand type for table SBS192.")
        return (musica)
      } else if (strand_type == "Transcript_Strand") {
        annotate_transcript_strand(musica, "19")
        return (musica)
      } else {
        annotate_replication_strand(musica, "19")
        return(musica)
      }
    }
    build_standard_table(musica, select_genome("19"),
                         table_name = input$SelectTable, 
                         #strand_type = strand_type,
                         overwrite = input$OverwriteTable == TRUE)
    return(musica)
  }
  shinyalert::shinyalert(title = "Oops",
                         text = "Custom tables are not yet supported.")
}
