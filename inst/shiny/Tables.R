add_tables <- function (input) { 
  strand_type <- input$StrandType
  # Check it table already exists
  if(input$SelectTable %in% names(extract_count_tables(musica)) &&
     input$OverwriteTable == F) {
    shinyalert::shinyalert(title = "Warning", 
                           text = "Table exists but was not overwritten. No changed were made.")
    return ()
  }
  if (input$SelectTable != "Custom") {
    # Check inputs for SBS192
    if (input$SelectTable == "SBS192") {
      if (strand_type == "") {
        shinyalert::shinyalert(title = "Oops", 
                               text = "You must select strand type for table SBS192.")
        return ()
      } else if (strand_type == "Transcript_Strand") {
        annotate_transcript_strand(musica, input$Genome)
      }
    }
    print(input$OverwriteTable)
    build_standard_table(musica, g, table_name = input$SelectTable, 
                         overwrite = input$OverwriteTable == TRUE)
    return()
  }
  shinyalert::shinyalert(title = "Oops",
                         text = "Custom tables are not yet supported.")
}
