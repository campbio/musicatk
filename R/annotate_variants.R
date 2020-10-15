#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom S4Vectors queryHits decode subjectHits
#' @importFrom GenomicRanges strand
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
NULL

#' Uses a genome object to find context and add it to the variant table
#'
#' @param musica Input samples
#' @param g A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from.
#' @param flank_start Start of flank area to add, can be positive or negative
#' @param flank_end End of flank area to add, can be positive or negative
#' @param build_table Automatically build a table using the annotation and add
#' @param overwrite Overwrite existing count table
#' @return None
#' it to the musica
#' @examples
#' data(musica_sbs96_tiny)
#' g <- select_genome("19")
#' add_flank_to_variants(musica_sbs96_tiny, g, 1, 2)
#' add_flank_to_variants(musica_sbs96_tiny, g, -2, -1)
#' @export
add_flank_to_variants <- function(musica, g, flank_start, flank_end,
                                  build_table = TRUE, overwrite = FALSE) {
  stopifnot(sign(flank_start) == sign(flank_end), flank_start < flank_end)

  direction <- ifelse(sign(flank_start) == 1, "r", "l")

  #Determine output and calculations based on selected context area
  output_column <- paste(direction, "flank_", abs(flank_start), "_to_",
                         abs(flank_end), sep = "")

  dat <- musica@variants
  mut_type <- paste(dat$ref, ">", dat$alt, sep = "")
  chr <- dat$chr

  if (sign(flank_start) == 1) {
    center <- dat$end
  } else {
    center <- dat$start
  }
  ref <- dat$ref
  type <- mut_type

  #Mutation Context
  flank <- VariantAnnotation::getSeq(g, chr, center + flank_start,
                                     center + flank_end, as.character = TRUE)
  final_mut_context <- rep(NA, length(ref))

  # Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- type %in% forward_change

  final_mut_context[ind] <- flank[ind]

  # Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- type %in% rev_change

  # Reverse complement the context so only 6 mutation categories instead of 12
  rev_flank <- flank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement()

  final_mut_context[ind] <- as.character(rev_flank)
  dat[[output_column]] <- final_mut_context
  eval.parent(substitute(musica@variants <- dat))
  if (build_table) {
    dat_musica <- methods::new("musica", variants = dat, count_tables =
                               musica@count_tables,
                             sample_annotations = musica@sample_annotations)
    tab <- build_custom_table(dat_musica, variant_annotation = output_column,
                         name = output_column, return_instead = FALSE,
                         overwrite = overwrite)
    eval.parent(substitute(musica@count_tables[[output_column]] <- tab))
  }
}

#' Adds an annotation to the input musica's variant table with length of each
#' variant
#'
#' @param musica Input samples
#' @return None
#' @examples
#' data(musica)
#' annotate_variant_length(musica)
#' musica
#' @export
annotate_variant_length <- function(musica) {
  dat <- musica@variants
  var_length <- rep(NA, nrow(dat))
  var_length[which(dat$Variant_Type == "SBS")] <- 1
  var_length[which(dat$Variant_Type == "DBS")] <- 2
  indels <- which(dat$Variant_Type == "indel")
  var_length[indels] <- nchar(dat$alt[indels]) -
    nchar(dat$ref[indels])
  dat[["Variant_Length"]] <- var_length
  eval.parent(substitute(musica@variants <- dat))
}

#' Drops a column from the variant table that the user no longer needs
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param column_name Name of column to drop
#' @return None
#' @examples
#' data(musica)
#' drop_annotation(musica, "Variant_Type")
#' @export
drop_annotation <- function(musica, column_name) {
  dat <- musica@variants
  stopifnot(column_name %in% colnames(dat))
  data.table::set(dat, j = column_name, value = NULL)
  eval.parent(substitute(musica@variants <- dat))
}

#' Generates a variant type table
#'
#' @param tab Input variant table
#' @return Returns the inputted variant table with variant type ("SBS", "DBS",
#' "INS", "DEL") added as an appended "Variant_Type" column
#' @examples
#' data(musica)
#' variants <- get_variants(musica)
#' musicatk:::add_variant_type(variants)
#' @keywords internal
add_variant_type <- function(tab) {
  type <- rep(NA, nrow(tab))
  type[which(nchar(tab$ref) == 1 &
               nchar(tab$alt) == 1)] <- "SBS"
  type[which(nchar(tab$ref) == 2 &
               nchar(tab$alt) == 2)] <- "DBS"
  type[which(tab$ref == "-")] <- "INS"
  type[which(tab$alt == "-")] <- "DEL"
  type[which(is.na(type))] <- "unknown"
  tab[["Variant_Type"]] <- type
  return(tab)
}

#' Annotate variants with variant type ("SBS", "INS", "DEl", "DBS")
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @return None
#' @examples
#' data(musica)
#' annotate_variant_type(musica)
#' @export
annotate_variant_type <- function(musica) {
  type_added <- add_variant_type(musica@variants)
  eval.parent(substitute(musica@variants <- type_added))
}

#' Subsets a variant table based on Variant Type
#'
#' @param tab Input variant table
#' @param type Variant type to return e.g. "SBS", "INS", "DEL", "DBS"
#' @return Returns the input variant table subsetted to only contain variants
#' of the specified variant type
#' @examples
#' data(musica)
#' annotate_variant_type(musica)
#' subset_variant_by_type(get_variants(musica), "SBS")
#' @export
subset_variant_by_type <- function(tab, type) {
  if (!"Variant_Type" %in% colnames(tab)) {
    stop(paste("No Variant_Type annotation found, ",
               "please run annotate_variant_type first."))
  }
  if (!any(tab$Variant_Type %in% type)) {
    stop(paste("No variants of type: ", type))
  }
  return(tab[which(tab$Variant_Type %in% type), ])
}

#' Add transcript strand annotation to SBS variants (defined in genes only)
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param genome_build Which genome build to use: hg19, hg38, or a custom TxDb
#' object
#' @param build_table Automatically build a table from this annotation
#' @return None
#' @examples
#' data(musica)
#' annotate_transcript_strand(musica, 19)
#' @export
annotate_transcript_strand <- function(musica, genome_build, build_table = TRUE) {
  if (genome_build %in% c("19", "hg19")) {
    genes <- genes(
      TxDb.Hsapiens.UCSC.hg19.knownGene)
  } else if (genome_build %in% c("38", "hg38")) {
    genes <- genes(
      TxDb.Hsapiens.UCSC.hg38.knownGene)
  } else if (methods::isClass(genome_build, "TxDb")) {
    genes <- genome_build
  } else {
    stop("Please select hg19, hg38, or provde a TxDb object")
  }

  dat <- musica@variants
  sbs_index <- which(dat$Variant_Type == "SBS")
  sbs <- subset_variant_by_type(dat, "SBS")

  #Create VRanges object to determine strand of variants within genes
  vrange <- VariantAnnotation::VRanges(seqnames = sbs$chr, ranges =
                                         IRanges(sbs$start,
                                                 sbs$end), ref =
                                         sbs$ref, alt =
                                         sbs$alt)
  overlaps <- findOverlaps(vrange, genes)
  transcribed_variants <- rep("NA", nrow(dat))
  transcribed_variants[sbs_index[queryHits(overlaps)]] <- as.character(decode(
    strand(genes[subjectHits(overlaps)])))

  #Match transcription and +, -, to account for reverse complement
  mut_type <- paste(dat$ref, ">", dat$alt,
                            sep = "")
  final_strand <- rep(NA, nrow(dat))
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- mut_type %in% forward_change
  final_strand[ind] <- ifelse(transcribed_variants[ind] == "+", "U", "T")
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- mut_type %in% rev_change
  final_strand[ind] <- ifelse(transcribed_variants[ind] == "-", "U", "T")

  dat[["Transcript_Strand"]] <- final_strand
  eval.parent(substitute(musica@variants <- dat))
  if (build_table) {
    dat_musica <- methods::new("musica", variants = drop_na_variants(
      dat, "Transcript_Strand"), count_tables = musica@count_tables,
      sample_annotations = musica@sample_annotations)
    tab <- build_custom_table(musica = dat_musica, variant_annotation =
                                  "Transcript_Strand", name =
                                  "Transcript_Strand", return_instead = TRUE)
    eval.parent(substitute(musica@count_tables[["Transcript_Strand"]] <- tab))
  }
}

#' Add replication strand annotation to SBS variants based on bedgraph file
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param rep_range A GRanges object with replication timing as metadata
#' @param build_table Automatically build a table from this annotation
#' @return None
#' @examples
#' data(musica)
#' data(rep_range)
#' annotate_replication_strand(musica, rep_range)
#' @export
annotate_replication_strand <- function(musica, rep_range, build_table = TRUE) {
  dat <- musica@variants
  sbs_index <- which(dat$Variant_Type == "SBS")
  sbs <- subset_variant_by_type(dat, "SBS")

  #Create GRanges object to determine strand of variants within genes
  dat_range <- GenomicRanges::GRanges(seqnames = sbs$chr, ranges =
                                         IRanges(sbs$start,
                                                 sbs$end), ref =
                                         sbs$ref, alt =
                                         sbs$alt)
  overlaps <- GenomicRanges::findOverlaps(dat_range, rep_range)
  repl_variants <- rep("NA", length(dat_range))
  repl_variants[sbs_index[S4Vectors::queryHits(overlaps)]] <-
    GenomicRanges::elementMetadata(rep_range)@listData[[1]][
      S4Vectors::subjectHits(overlaps)]
  dat[["Replication_Strand"]] <- repl_variants
  eval.parent(substitute(musica@variants <- dat))
  if (build_table) {
    dat_musica <- methods::new("musica", variants = drop_na_variants(
      dat, "Replication_Strand"), count_tables = musica@count_tables,
      sample_annotations = musica@sample_annotations)
    tab <- build_custom_table(dat_musica, variant_annotation =
                                "Replication_Strand", name =
                                "Replication_Strand", data_factor =
                                factor(c("leading", "lagging")),
                              return_instead = FALSE)
    eval.parent(substitute(musica@count_tables[["Replication_Strand"]] <- tab))
  }
}
