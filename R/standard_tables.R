#' Uses a genome object to find context and generate standard SBS96 tables
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param g A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from.
#' @param overwrite Overwrite existing count table
#' @return Returns the created SBS96 count table object
#' @keywords internal
create_sbs96_table <- function(musica, g, overwrite = FALSE) {
  dat <- subset_variant_by_type(variants(musica), type = "SBS")
  ref <- as.character(dat$ref)
  alt <- as.character(dat$alt)
  mut_type <- paste(ref, ">", alt, sep = "")

  #Mutation Context
  context <- BSgenome::getSeq(x = g, names = as.character(dat$chr),
                              start = dat$start - 1,
                              end = dat$end + 1,
                              as.character = TRUE)

  final_mut_type <- rep(NA, length(ref))
  final_mut_context <- rep(NA, length(ref))

  # Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- mut_type %in% forward_change
  final_mut_type[ind] <- paste(as.character(ref[ind]), ">",
                               as.character(alt[ind]), sep = "")
  final_mut_context[ind] <- context[ind]

  # Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- mut_type %in% rev_change

  # Reverse complement the context so only 6 mutation categories instead of 12
  rev_context <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(context[ind])))
  rev_refbase <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(ref[ind])))
  rev_altbase <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(alt[ind])))

  final_mut_type[ind] <- paste0(rev_refbase, ">", rev_altbase)
  final_mut_context[ind] <- rev_context
  final_motif <- paste0(final_mut_type, "_", final_mut_context)

  ###### Now we separate into samples
  ## Define all mutation types for 96 substitution scheme
  b1 <- rep(rep(c("A", "C", "G", "T"), each = 4), 6)
  b2 <- rep(c("C", "T"), each = 48)
  b3 <- rep(c("A", "C", "G", "T"), 24)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut <- rep(forward_change, each = 16)
  annotation <- data.frame("motif" = paste0(mut, "_", mut_trinuc),
                           "mutation" = mut,
                           "context" = mut_trinuc)
  rownames(annotation) <- annotation$motif


  mut_id <- apply(cbind(mut, mut_trinuc), 1, paste,
                  collapse = "_")
  mutation <- factor(final_motif, levels = annotation$motif)

  mut_table <- as.matrix(as.data.frame.matrix(xtabs(~ mutation + dat$sample)))
  #xtabs adds dat$sample as dimname[2] so we remove it
  dimnames(mut_table) <- list(rownames(mut_table), colnames(mut_table))

  # Use COSMIC color scheme
  color_mapping <- c("C>A" = "#5ABCEBFF",
                     "C>G" = "#050708FF",
                     "C>T" = "#D33C32FF",
                     "T>A" = "#CBCACBFF",
                     "T>C" = "#ABCD72FF",
                     "T>G" = "#E7C9C6FF")

# Need to think about this more carefully - what happens if indel are zero but
#  not SNV and the user wants to combine them?
#  zero_samps <- which(colSums(mut_table) == 0)
#  if (length(zero_samps) > 0) {
#    warning(paste0("Dropping the following zero count samples: ",
#                   paste(names(zero_samps), collapse = ", ")))
#    mut_table <- mut_table[, -zero_samps, drop = FALSE]
#  }
  tab <- .create_count_table(musica = musica, name = "SBS96",
                            count_table = mut_table,
                            annotation = annotation,
                            features = data.frame(mutation = final_motif),
                            type = as.character(dat$Variant_Type),
                            color_variable = "mutation",
                            color_mapping = color_mapping,
                     description = paste0("Single Base Substitution table with",
                     " one base upstream and downstream"), return_table = TRUE,
                     overwrite = overwrite)
  return(tab)
}

#' Uses a genome object to find context and generate standard SBS192 table
#' using transcript strand
#'
#' @param musica Input samples
#' @param g A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from.
#' @param strand_type Transcript_Strand or Replication_Strand
#' @param overwrite Overwrite existing count table
#' @return Returns the created SBS192 count table object built using either
#' transcript strand or replication strand
#' @keywords internal
create_sbs192_table <- function(musica, g, strand_type, overwrite = FALSE) {
  if (!strand_type %in% c("Transcript_Strand", "Replication_Strand")) {
    stop("Please select either Transcript_Strand or Replication_Strand")
  }
  dat <- variants(musica)
  dat <- subset_variant_by_type(dat, "SBS")
  dat <- drop_na_variants(dat, strand_type)

  chr <- dat$chr
  range_start <- dat$start
  range_end <- dat$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 1, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 1,
                                      as.character = TRUE)
  ref_context <- paste(lflank, dat$ref, rflank, sep = "")

  final_mut_type <- rep(NA, nrow(dat))
  final_mut_context <- rep(NA, nrow(dat))

  ## Get mutation type
  initial_maf_type <- paste(dat$ref, ">", dat$alt,
                           sep = "")

  ## Get mutation context info for those on "+" strand
  forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind <- dat$Variant_Type == "SBS" & initial_maf_type %in% forward_change

  final_mut_type[ind] <- initial_maf_type[ind]
  final_mut_context[ind] <- ref_context[ind]

  ## Get mutation context info for those on "-" strand
  rev_change <- c("A>G", "A>T", "A>C", "G>T", "G>C", "G>A")
  ind <- dat$Variant_Type == "SBS" & initial_maf_type %in% rev_change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev_context <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    ref_context[ind]))
  rev_refbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$ref[ind]))
  rev_altbase <- Biostrings::reverseComplement(Biostrings::DNAStringSet(
    dat$alt[ind]))

  final_mut_type[ind] <- paste(as.character(rev_refbase), ">",
                              as.character(rev_altbase), sep = "")
  final_mut_context[ind] <- rev_context

  maf_mut_id <- paste(final_mut_type, final_mut_context, dat[[strand_type]],
                     sep = "_")
  tumor_id <- dat$sample

  ## Define all mutation types for 196 substitution scheme
  b1 <- rep(rep(c("A", "C", "G", "T"), each = 24), 2)
  b2 <-  rep(rep(c("C", "T"), each = 12), 8)
  b3 <- rep(c("A", "C", "G", "T"), 48)
  mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
  mut_type <- rep(rep(rep(forward_change, each = 4), 4), 2)
  if (strand_type == "Transcript_Strand") {
    mut_strand <- rep(c("T", "U"), each = 96)
  } else if (strand_type == "Replication_Strand") {
    mut_strand <- rep(c("leading", "lagging"), each = 96)
  }
  mut_id <- apply(cbind(mut_type, mut_trinuc, mut_strand), 1, paste,
                 collapse = "_")

  mutation <- factor(maf_mut_id, levels = mut_id)

  mut_table <- xtabs(~ mutation + tumor_id)

  #Convert to table by dropping xtabs class and call
  attr(mut_table, "call") <- NULL
  attr(mut_table, "class") <- NULL

  annotation <- data.frame(motif = mut_id, mutation =
                             unlist(lapply(strsplit(mut_id, "_"), "[[", 1)),
                           context = paste(unlist(lapply(strsplit(mut_id, "_"),
                                                   "[[", 2)),
                                           unlist(lapply(strsplit(mut_id, "_"),
                                                         "[[", 3)), sep = "_"),
                           row.names = mut_id)
  color_mapping <- .gg_color_hue(length(unique(annotation$mutation)))
  names(color_mapping) <- unique(annotation$mutation)

  tab <- .create_count_table(musica = musica, name = paste0("SBS192_", ifelse(
    strand_type == "Transcript_Strand", "Trans", "Rep")),
                             count_table = mut_table,
                             annotation = annotation,
                             features = data.frame(mutation = maf_mut_id),
                             type = paste(as.character(dat$Variant_Type),
                                          ifelse(strand_type ==
                                                   "Transcript_Strand",
                                                 "Trans", "Rep"), sep = "_"),
                             color_variable = "mutation",
                             color_mapping = color_mapping,
                             description = paste0("Single Base Substitution ",
                                                  "table with one base ",
                                                  "upstream and downstream and",
                                                  " transcript strand"),
                             return_table = TRUE, overwrite = overwrite)
  return(tab)
}

#' Creates and adds a table for standard doublet base subsitution (DBS)
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param overwrite Overwrite existing count table
#' @return Returns the created DBS table object
#' @keywords internal
create_dbs78_table <- function(musica, overwrite = overwrite) {
  dbs <- subset_variant_by_type(variants(musica), "DBS")

  ref <- dbs$ref
  alt <- dbs$alt

  #Reverse Complement broad categories to the other strand
  rc_ref <- which(ref %in% c("GT", "GG", "AG", "GA", "CA", "AA"))
  ref[rc_ref] <- rc(ref[rc_ref])
  alt[rc_ref] <- rc(alt[rc_ref])

  #Reverse complement more specific categories to the other strand
  rc_alt <- NULL
  rc_ref_AT <- which(ref == "AT")
  rc_alt <- c(rc_alt, rc_ref_AT[which(alt[rc_ref_AT] %in% c("TG", "GG", "TC"))])

  rc_ref_CG <- which(ref == "CG")
  rc_alt <- c(rc_alt, rc_ref_CG[which(alt[rc_ref_CG] %in% c("AC", "GA", "AA"))])

  rc_ref_GC <- which(ref == "GC")
  rc_alt <- c(rc_alt, rc_ref_GC[which(alt[rc_ref_GC] %in% c("TT", "CT", "TG"))])

  rc_ref_TA <- which(ref == "TA")
  rc_alt <- c(rc_alt, rc_ref_TA[which(alt[rc_ref_TA] %in% c("AG", "CC", "AC"))])

  if (length(rc_alt) > 0) {
    alt[rc_alt] <- rc(alt[rc_alt])
  }

  full <- paste(ref, ">NN_", alt, sep = "")
  full_motif <- c(paste0("AC>NN", "_", c("CA", "CG", "CT", "GA", "GG", "GT",
                                         "TA", "TG", "TT")),
            paste0("AT>NN", "_", c("CA", "CC", "CG", "GA", "GC", "TA")),
            paste0("CC>NN", "_", c("AA", "AG", "AT", "GA", "GG", "GT", "TA",
                                   "TG", "TT")),
            paste0("CG>NN", "_", c("AT", "GC", "GT", "TA", "TC", "TT")),
            paste0("CT>NN", "_", c("AA", "AC", "AG", "GA", "GC", "GG", "TA",
                                   "TC", "TG")),
            paste0("GC>NN", "_", c("AA", "AG", "AT", "CA", "CG", "TA")),
            paste0("TA>NN", "_", c("AT", "CG", "CT", "GC", "GG", "GT")),
            paste0("TC>NN", "_", c("AA", "AG", "AT", "CA", "CG", "CT", "GA",
                                   "GG", "GT")),
            paste0("TG>NN", "_", c("AA", "AC", "AT", "CA", "CC", "CT", "GA",
                                   "GC", "GT")),
            paste0("TT>NN", "_", c("AA", "AC", "AG", "CA", "CC", "CG", "GA",
                                   "GC", "GG")))

  sample_names <- sample_names(musica)
  num_samples <- length(sample_names)
  variant_tables <- vector("list", length = num_samples)
  for (i in seq_len(num_samples)) {
    sample_index <- which(dbs$sample == sample_names[i])
    if (length(sample_index) > 0) {
      variant_tables[[i]] <- table(factor(full[sample_index],
                                          levels = full_motif))
    } else {
      variant_tables[[i]] <- rep(0, length(full_motif))
    }
  }
  mut_table <- do.call(cbind, variant_tables)
  colnames(mut_table) <- sample_names
  
  full_tab <- matrix(data = 0, nrow = length(full_motif), 
                     ncol = length(levels(dbs$sample)), 
                     dimnames = list(full_motif, levels(dbs$sample)))
  full_tab[, colnames(mut_table)] <- mut_table

  annotation <- data.frame(motif = full_motif, mutation =
                             unlist(lapply(strsplit(full_motif, "_"), "[[", 1)),
                           context = unlist(lapply(strsplit(full_motif, "_"),
                                                   "[[", 2)),
                           row.names = full_motif)
  color_mapping <- .gg_color_hue(length(unique(annotation$mutation)))
  names(color_mapping) <- unique(annotation$mutation)
  tab <- .create_count_table(musica = musica, name = "DBS78",
                             count_table = full_tab,
                             annotation = annotation,
                             features = data.frame(mutation = full),
                             type = as.character(dbs$Variant_Type),
                             color_variable = "mutation",
                             color_mapping = color_mapping,
                             description = paste0("Standard count table for ",
                                                  "double base substitutions", 
                                                  "using COSMIC v3 schema"),
                             return_table = TRUE, overwrite = overwrite)
  return(tab)
}

.gg_color_hue <- function(n) {
  hues <- base::seq(15, 375, length = n + 1)
  return(grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)])
}

#' Reverse complement of a string using biostrings
#'
#' @param dna Input DNA string
#' @return Returns the reverse compliment of the input DNA string
#' @examples
#' rc("ATGC")
#' @export
rc <- function(dna) {
  if (is(dna, "character") && length(dna) == 1) {
    rev_com <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAString(dna)))
  } else if (is(dna, "character") && length(dna) > 1) {
    rev_com <- vapply(dna, rc, FUN.VALUE = character(1))
    names(rev_com) <- NULL
  } else {
    stop("Must be character or character vector")
  }
  return(rev_com)
}

#' Builds count tables using various mutation type schemas
#' 
#' Generates count tables for different mutation type schemas which can be 
#' used as input to the mutational signature discovery or prediction functions.
#' \code{"SBS96"} generates a table for single base substitutions following the
#' standard 96 mutation types derived from the trinucleotide context.
#' \code{"SBS192"} is the 96 mutation type schema with the addition of
#' transcriptional strand or replication strand information added to each base.
#' \code{"DBS"} generates a table for the double base substitution schema 
#' used in COSMIC V3. \code{"Indel"} generates a table for insertions and 
#' deletions following the schema used in COSMIC V3.
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param g A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from.
#' @param table_name Name of standard table to build. One of \code{"SBS96"},
#' \code{"SBS192"}, \code{"DBS"}, or \code{"Indel"}.
#' @param strand_type Strand type to use in SBS192 schema. One of
#' \code{"Transcript_Strand"} or \code{"Replication_Strand"}.
#' Only used if \code{table_name = SBS192}.  
#' @param overwrite If \code{TRUE}, any existing count table with the same
#' name will be overwritten. If \code{FALSE}, then an error will be thrown
#' if a table with the same name exists within the \code{musica} object.
#' @param verbose Show progress bar for processed samples
#' @return No object will be returned. The count tables will be automatically
#' added to the \code{musica} object. 
#' @examples
#' g <- select_genome("19")
#'
#' data(musica)
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#'
#' data(musica)
#' annotate_transcript_strand(musica, "19")
#' build_standard_table(musica, g, "SBS192", "Transcript_Strand")
#'
#' data(musica)
#' data(rep_range)
#' annotate_replication_strand(musica, rep_range)
#' build_standard_table(musica, g, "SBS192", "Replication_Strand")
#'
#' data(dbs_musica)
#' build_standard_table(dbs_musica, g, "DBS", overwrite = TRUE)
#'
#' data(indel_musica)
#' build_standard_table(indel_musica, g, table_name = "INDEL")
#' @export
build_standard_table <- function(musica, g, table_name, strand_type = NULL,
                                 overwrite = FALSE, verbose = FALSE) {
  if (table_name %in% c("SNV96", "SNV", "96", "SBS", "SBS96")) {
    message("Building count table from SBS with SBS96 schema")
    .table_exists_warning(musica, "SBS96", overwrite)
    tab <- create_sbs96_table(musica, g, overwrite)
  } else if (table_name %in% c("SBS192", "192")) {
    message("Building count table from SBS and ", strand_type, 
            " with SBS192 schema")
    .table_exists_warning(musica, "SBS192", overwrite)
    tab <- create_sbs192_table(musica, g, strand_type, overwrite)
  } else if (table_name %in% c("DBS78", "DBS", "doublet")) {
    message("Building count table from DBS with DBS78 schema")
    .table_exists_warning(musica, "DBS78", overwrite)
    tab <- create_dbs78_table(musica, overwrite)
  } else if (table_name %in% c("IND83", "INDEL83", "INDEL", "IND", "indel", 
                               "Indel")) {
    message("Building count table from INDELs with IND83 schema")
    .table_exists_warning(musica, "IND83", overwrite)
    tab <- create_ind83_table(musica, g, overwrite, verbose)
  } else {
    stop(paste0("There is no standard table named: ", table_name,
               " please select from SBS96, SBS192, DBS78, IND83."))
  }
  tab_list <- tables(musica)
  tab_list[[get_tab_name(tab)]] <- tab
  eval.parent(substitute(tables(musica) <- tab_list))
}

.table_exists_warning <- function(musica, table_name, overwrite = FALSE) {
  if (table_name %in% names(musica@count_tables)) {
    if (!overwrite) {
      stop(paste0("Table: ", table_name,
                  " already exists, use overwrite to continue."))
    } else {
      warning(paste0("Overwriting counts table: ", table_name))
    }
  }
}

create_ind83_table <- function(musica, g, overwrite = FALSE, 
                                 verbose = FALSE) {
  var <- variants(musica)
  all_ins <- as.data.frame(subset_variant_by_type(var, "INS"))
  all_del <- as.data.frame(subset_variant_by_type(var, "DEL"))
  #samples <- unique(c(as.character(all_ins$sample),
  #                    as.character(all_del$sample)))
  all_samples <- sample_names(musica)
  dimlist <- list(row_names = c(.get_indel_motifs("bp1", 0, 0),
                                .get_indel_motifs("bp1", 1, 0),
                                .get_indel_motifs("del", NA, NA),
                                .get_indel_motifs("ins", NA, NA),
                                .get_indel_motifs("micro", NA, NA)),
                  column_names = all_samples)
  mut_table <- matrix(NA, nrow = 83, ncol = length(all_samples), 
                      dimnames = dimlist)
  if (isTRUE(verbose)) {
    pb <- utils::txtProgressBar(min = 0, max = length(all_samples), initial = 0,
                                style = 3)
    i <- 0
    
  }
  for (sample in all_samples) {
    ins <- all_ins[which(all_ins$sample == sample), ]
    del <- all_del[which(all_del$sample == sample), ]

    ins_len <- nchar(ins$alt)
    del_len <- nchar(del$ref)
    ins1 <- ins[which(ins_len == 1), ]
    ins2 <- ins[which(ins_len > 1), ]

    del1 <- del[which(del_len == 1), ]
    del2 <- del[which(del_len > 1), ]

    if (nrow(del1) == 0) {
      del1_counts <- setNames(rep(0, 12), .get_indel_motifs("bp1", 0, 0))
    } else {
      del1_counts <- .count1(mut = del1, type = del1$ref, ins = FALSE, g = g)
    }

    if (nrow(ins1) == 0) {
      ins1_counts <- setNames(rep(0, 12), .get_indel_motifs("bp1", 1, 0))
    } else {
      ins1_counts <- .count1(mut = ins1, type = ins1$alt, ins = TRUE, g = g)
    }

    if (nrow(ins2) == 0) {
      ins2_counts <- setNames(rep(0, 24), .get_indel_motifs("ins", NA, NA))
    } else {
      ins2_counts <- .count2_ins(mut = ins2, type = ins2$alt, g = g)
    }

    if (nrow(del2) == 0) {
      del2_counts <- list(del = setNames(rep(0, 24),
                                         .get_indel_motifs("del", NA, NA)),
                          micro = setNames(rep(0, 11),
                                           .get_indel_motifs("micro", NA, NA)))
    } else {
      del2_counts <- .count2_del(mut = del2, type = del2$ref, g)
    }
    mut_table[, sample] <- c(del1_counts, ins1_counts, del2_counts$del,
                             ins2_counts, del2_counts$micro)
    if (isTRUE(verbose)) {
      i <- i + 1
      utils::setTxtProgressBar(pb, i)
    }
  }

  motif <- rownames(mut_table)
  mutation <- c(substr(motif[seq_len(24)], 1, 5),
                paste(unlist(lapply(strsplit(motif[25:83], "_"), "[[", 1)),
                      unlist(lapply(strsplit(motif[25:83], "_"), "[[", 2)),
                      unlist(lapply(strsplit(motif[25:83], "_"), "[[", 3)),
                      sep = "_"))
  context <- c(substr(motif[seq_len(24)], 7, 9),
               unlist(lapply(strsplit(motif[25:83], "_"), "[[", 3)))
  annotation <- data.frame(motif = motif, mutation = mutation,
                           context = context)

  color_mapping <- .gg_color_hue(length(unique(annotation$mutation)))
  names(color_mapping) <- unique(annotation$mutation)

  #TODO error in counting, we're missing some
  #incorrect_features <- length(rep(rownames(mut_table), rowSums(mut_table)))
  dummy_features <- rep(NA, length(which(var$Variant_Type %in% c("INS", 
                                                                 "DEL"))))

  tab <- .create_count_table(musica = musica, name = "IND83",
                             count_table = mut_table,
                             annotation = annotation,
                             features = data.frame(mutation = dummy_features),
                             type = Rle(rep("INDEL", length(which(
                               var$Variant_Type %in% c("INS", "DEL"))))),
                             color_variable = "mutation",
                             color_mapping = color_mapping,
                             description = paste0("Standard count table for ",
                                                  "small insertions and",
                                                  " deletions"),
                             return_table = TRUE, overwrite = overwrite)
  return(tab)
}

.get_indel_motifs <- function(indel, ins, plus) {
  if (indel == "bp1") {
    return(paste(ifelse(ins, "INS", "DEL"), c("C", "C", "C", "C", "C", "C", "T",
                                              "T", "T", "T", "T", "T"), 1,
                 c(c(0, 1, 2, 3, 4, paste0(5 + plus, "+"), 0, 1, 2, 3, 4,
                     paste0(5 + plus, "+"))), sep = "_"))
  } else if (indel == "ins") {
    return(paste(paste("INS_repeats", c("2", "2", "2", "2", "2", "2", "3", "3", "3", "3",
                                        "3", "3", "4", "4", "4", "4", "4", "4", "5+",
                                        "5+", "5+", "5+", "5+", "5+"),
                       c("0", "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4",
                         "5+", "0", "1", "2", "3", "4", "5+", "0", "1", "2", "3",
                         "4", "5+"), sep = "_")))
  } else if (indel == "micro") {
    return(paste("DEL_MH",  c("2", "3", "3", "4", "4", "4", "5+", "5+",
                              "5+", "5+", "5+"),
                 c("1", "1", "2", "1", "2", "3", "1", "2",
                   "3", "4", "5+"), sep = "_"))
  } else if (indel == "del") {
    return(paste("DEL_repeats", c("2", "2", "2", "2", "2", "2", "3", "3", "3", "3",
                                  "3", "3", "4", "4", "4", "4", "4", "4", "5+",
                                  "5+", "5+", "5+", "5+", "5+"),
                 c("0", "1", "2", "3", "4", "5+", "0", "1", "2", "3", "4",
                   "5+", "0", "1", "2", "3", "4", "5+", "0", "1", "2", "3",
                   "4", "5+"), sep = "_"))
  } else {
    stop("Unrecognized indel type")
  }
}

.count1 <- function(mut, type, ins, g) {
  #ifelse(ins, plus <- 0, plus <- 0)
  ifelse(ins, plus <- 0, plus <- 0)
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 7, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 7,
                                      as.character = TRUE)
  final_lflank <- lflank
  final_rflank <- rflank
  ind <- which(type %in% c("A", "G"))
  final_lflank[ind] <- rflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  final_rflank[ind] <- lflank[ind] %>%
    Biostrings::DNAStringSet() %>%
    Biostrings::reverseComplement() %>%
    as.character()
  final_type <- type
  final_type[ind] <- final_type[ind] %>% Biostrings::DNAStringSet() %>% 
    Biostrings::reverseComplement() %>% as.character()
  repeats <- rep(NA, length(final_type))
  for (i in seq_len(length(type))) {
    repeats[i] <- .count_repeat(letter = final_type[i], string = final_rflank[i]) +
      .count_repeat(final_type[i], stringi::stri_reverse(final_lflank[i])) + plus
  }
  repeats[repeats >= 5 + plus] <- paste0(5 + plus, "+")
  bp1_motif <- .get_indel_motifs("bp1", ins, plus)
  return(table(factor(paste(ifelse(ins, "INS", "DEL"), final_type, 1,
                            repeats, sep = "_"), levels = bp1_motif)))
}

.count2_ins <- function(mut, type, g) {
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - 10, range_start - 1,
                                      as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1, range_end + 10,
                                      as.character = TRUE)
  repeats <- rep(NA, length(type))
  for (i in seq_len(length(type))) {
    repeats[i] <- .count_repeat(type[i], rflank[i]) +
      .count_repeat(type[i], stringi::stri_reverse(lflank[i]))
  }
  repeats[repeats >= 5] <- paste0(5, "+")
  len <- nchar(type)
  len[which(len >= 5)] <- "5+"
  ins_motif <- .get_indel_motifs("ins", NA, NA)
  return(table(factor(paste("INS_repeats", len, repeats, sep = "_"),
                      levels = ins_motif)))
}

.count2_del <- function(mut, type, g) {
  chr <- mut$chr
  range_start <- mut$start
  range_end <- mut$end

  len <- nchar(type)
  lflank <- VariantAnnotation::getSeq(g, chr, range_start - len * 5,
                                      range_start - 1, as.character = TRUE)
  rflank <- VariantAnnotation::getSeq(g, chr, range_end + 1,
                                      range_end + len * 5, as.character = TRUE)
  #has_repeat <- type == lflank | type == rflank #could be used for a speedup
  #maybe_micro <- which(!has_repeat) #doesn't get used from some reason #TODO

  micro <- rep(NA, length(type))
  for (i in seq_len(length(type))) {
    micro[i] <- max(.micro_left(type[i], lflank[i]),
                    .micro_right(type[i], rflank[i]))
  }

  repeats <- rep(NA, length(type))
  for (i in seq_len(length(type))) {
    #TEST removing the +1 for dels to bring it in line with Alexandrov counting
    #repeats[i] <- .count_repeat(type[i], rflank[i]) +
    #  .count_repeat(type[i], rev(lflank[i]))
    repeats[i] <- .count_repeat(type[i], rflank[i]) +
      .count_repeat(type[i], stringi::stri_reverse(lflank[i]))
  }
  #micro_ind <- which(repeats == 0 & micro > 0)
  micro_ind <- which(repeats == 0 & micro > 0)
  repeat_ind <- which(repeats > 0 | (repeats == 0 & micro == 0))
  final_micro <- micro[micro_ind]
  final_repeats <- repeats[repeat_ind]
  #final_repeats[final_repeats >= 5] <- paste0(5, "+")
  final_repeats[final_repeats >= 5] <- paste0(5, "+")
  final_len <- len
  final_len[which(final_len >= 5)] <- "5+"
  final_micro[which(final_micro >= 5)] <- "5+"
  micro_motif <- .get_indel_motifs("micro", NA, NA)
  del_motif <- .get_indel_motifs("del", NA, NA)
  del_tab <- table(factor(paste("DEL_repeats", final_len[repeat_ind],
                                final_repeats, sep = "_"), levels = del_motif))
  micro_tab <- table(factor(paste("DEL_MH", final_len[micro_ind], final_micro,
                                  sep = "_"), levels = micro_motif))
  return(list(del = del_tab, micro = micro_tab))
}

.micro_right <-  function(letter, string) {
  len <- nchar(letter)
  elem <- len - 1
  while(elem > 0) {
    if (substr(letter, 1, elem) == substr(string, 1, elem)) {
      return(elem)
    } else {
      elem <- elem - 1
    }
  }
  return(0)
}

.micro_left <-  function(letter, string) {
  len <- nchar(letter)
  elem <- 1
  while(elem < len) {
    if (substr(letter, elem, len) == substr(string, elem, len)) {
      return(elem)
    } else {
      elem <- elem + 1
    }
  }
  return(0)
}

.count_repeat <- function(letter, string) {
  len <- nchar(letter)
  if (letter != substr(string, 1, len)) {
    return(0)
  } else {
    next_matches <- TRUE
    counts <- 1
    while(next_matches) {
      #message(paste0("Counts: ", counts))
      if (letter != substr(string, counts * len + 1, counts * len + len)) {
        return(counts)
      } else {
        counts <- counts + 1
      }
    }
  }
}
