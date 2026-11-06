# Writes a minimal VCF to a temporary file and returns its path.
# 'rows' are the body lines, each a character vector of the fields from
# CHROM to the last sample column. 'format' gives the FORMAT header lines;
# use character(0) for a VCF without FORMAT fields.
write_test_vcf <- function(rows, samples = c("S1", "S2"),
                           format = paste0(
                             "##FORMAT=<ID=GT,Number=1,Type=String,",
                             "Description=\"Genotype\">"
                           ),
                           info = character(0), fileext = ".vcf") {
  header <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    info,
    format,
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      if (length(samples) > 0) c("FORMAT", samples)
    ), collapse = "\t")
  )
  body <- vapply(rows, paste, character(1), collapse = "\t")
  path <- tempfile(fileext = fileext)
  writeLines(c(header, body), path)
  path
}
