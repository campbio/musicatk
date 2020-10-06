sig_compare <- function(sig1, sig2, threshold=0.9) {
  sig1_names <- colnames(sig1)
  sig2_names <- colnames(sig2)
  if (nrow(sig1) != nrow(sig2)) {
    stop("Signatures must have the same motifs")
  }
  matches <- matrix(nrow = ncol(sig1), ncol = ncol(sig2))
  for (i in seq_len(ncol(sig1))) {
    for (j in seq_len(ncol(sig2))) {
      matches[i, j] <- 1 - jsd(sig1[, i], sig2[, j])
    }
  }
  comparison <- NULL
  for (row in seq_len(nrow(matches))) {
    line <- which(matches[row, ] > threshold)
    if (length(line) > 0) {
      for (match in line) {
        comparison <- rbind(comparison, c(matches[row, match], row, match,
                                          sig1_names[row], sig2_names[match]))
      }
    }
  }
  if (is.null(comparison)) {
    stop("No matches found, try lowering threshold.")
  }
  comparison <- data.frame(comparison, stringsAsFactors = FALSE)
  colnames(comparison) <- c("cor", "xindex", "yindex", "xcol", "ycol")
  comparison$cor <- as.numeric(comparison$cor)
  comparison$xindex <- as.numeric(comparison$xindex)
  comparison$yindex <- as.numeric(comparison$yindex)
  return(comparison)
}


#' Compare two result files to find similar signatures
#'
#' @param result A \code{\linkS4class{musica_result}} object.
#' @param other_result A second \code{\linkS4class{musica_result}} object.
#' @param threshold threshold for similarity
#' @param result_name title for plot of first result signatures
#' @param other_result_name title for plot of second result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' compare_results(res, res, threshold = 0.8)
#' @export
compare_results <- function(result, other_result,
                            threshold = 0.9, result_name =
                              deparse(substitute(result)), other_result_name =
                              deparse(substitute(other_result))) {
  signatures <- result@signatures
  comparison <- sig_compare(signatures, other_result@signatures, threshold)
  result_subset <- methods::new("musica_result",
                                signatures = result@signatures[, comparison$xindex,
                                                               drop = FALSE], exposures =
                                  matrix(), type = "NMF", musica = result@musica,
                                tables = result@tables)
  other_subset <- methods::new("musica_result",
                               signatures = other_result@signatures[, comparison$yindex,
                                                                    drop = FALSE],
                               exposures = matrix(), type = "NMF",
                               musica = other_result@musica, tables = other_result@tables)
  result_plot <- plot_signatures(result_subset)
  result_plot <- result_plot + ggplot2::ggtitle(result_name)
  cosmic_plot <- plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle(other_result_name)
  gridExtra::grid.arrange(result_plot, cosmic_plot, ncol = 2)
  return(comparison)
}

#' Compare a result object to COSMIC V3 Signatures; Select exome or genome for
#' SBS and only genome for DBS or Indel classes
#'
#' @param result A \code{\linkS4class{musica_result}} object.
#' @param variant_class Compare to SBS, DBS, or Indel
#' @param sample_type exome (SBS only) or genome
#' @param threshold threshold for similarity
#' @param result_name title for plot user result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' compare_cosmic_v3(res, "SBS", "genome", threshold = 0.8)
#' @export
compare_cosmic_v3 <- function(result, variant_class, sample_type,
                              threshold = 0.9, result_name =
                                deparse(substitute(result))) {
  if (sample_type == "exome") {
    if (variant_class %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
      cosmic_res <- cosmic_v3_sbs_sigs_exome
    } else {
      stop(paste("Only SBS class is available for whole-exome, please choose",
                 " `genome` for DBS or Indel", sep = ""))
    }
  } else if (sample_type == "genome") {
    if (variant_class %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
      cosmic_res <- cosmic_v3_sbs_sigs
    } else if (variant_class %in% c("DBS", "dbs", "doublet")) {
      cosmic_res <- cosmic_v3_dbs_sigs
    } else if (variant_class %in% c("INDEL", "Indel", "indel", "ind", "IND",
                                    "ID")) {
      cosmic_res <- cosmic_v3_indel_sigs
    } else {
      stop("Only SBS, DBS, and Indel classes are supported")
    }
  } else {
    stop("Sample type must be exome or genome")
  }
  signatures <- result@signatures
  comparison <- sig_compare(signatures, cosmic_res@signatures, threshold)
  result_subset <- methods::new(
    "musica_result", signatures = result@signatures[, comparison$xindex, drop = FALSE],
    exposures = matrix(), type = "NMF", tables = result@tables,
    musica = result@musica)
  other_subset <- methods::new("musica_result", signatures =
                                 cosmic_res@signatures[, comparison$yindex,
                                                       drop = FALSE],
                               exposures = matrix(), type = "NMF", tables = cosmic_res@tables,
                               musica = cosmic_res@musica)
  result_plot <- plot_signatures(result_subset)
  result_plot <- result_plot + ggplot2::ggtitle(result_name)
  cosmic_plot <- plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle(paste("COSMIC Signatures v3",
                                                      variant_class, " ",
                                                      sample_type, sep = ""))
  gridExtra::grid.arrange(result_plot, cosmic_plot, ncol = 2)
  return(comparison)
}

#' Compare a result object to COSMIC V2 SBS Signatures (combination whole-exome
#' and whole-genome)
#'
#' @param result A \code{\linkS4class{musica_result}} object.
#' @param threshold threshold for similarity
#' @param result_name title for plot user result signatures
#' @return Returns the comparisons
#' @examples
#' res <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' compare_cosmic_v2(res, threshold = 0.7)
#' @export
compare_cosmic_v2 <- function(result, threshold = 0.9, result_name =
                                deparse(substitute(result))) {
  signatures <- result@signatures
  comparison <- sig_compare(signatures, cosmic_v2_sigs@signatures, threshold)
  result_subset <- methods::new("musica_result",
                                signatures =
                                  result@signatures[, comparison$xindex, drop =
                                                      FALSE], exposures =
                                  matrix(), type = result@type, musica = result@musica,
                                tables = result@tables)
  other_subset <- methods::new("musica_result",
                               signatures =
                                 cosmic_v2_sigs@signatures[, comparison$yindex,
                                                           drop = FALSE],
                               exposures = matrix(), type = "NMF",
                               musica = cosmic_v2_sigs@musica,
                               tables = cosmic_v2_sigs@tables)
  result_plot <- plot_signatures(result_subset)
  legend <- cowplot::get_legend(result_plot)
  result_plot <- result_plot + ggplot2::ggtitle(result_name) +
    theme(legend.position = "none")
  cosmic_plot <- plot_signatures(other_subset)
  cosmic_plot <- cosmic_plot + ggplot2::ggtitle("COSMIC Signatures v2") +
    theme(legend.position = "none")
  gridExtra::grid.arrange(result_plot, cosmic_plot, legend, ncol = 3,
                          widths = c(0.4, 0.4, 0.2))
  return(comparison)
}

#' Input a cancer subtype to return a list of related COSMIC signatures
#'
#' @param tumor_type Cancer subtype to view related signatures
#' @return Returns signatures related to a partial string match
#' @examples cosmic_v2_subtype_map ("lung")
#' @export
cosmic_v2_subtype_map <- function(tumor_type) {
  subtypes <- c("adrenocortical carcinoma", "all", "aml", "bladder", "breast",
                "cervix", "chondrosarcoma", "cll", "colorectum", "glioblastoma",
                "glioma low grade", "head and neck", "kidney chromophobe",
                "kidney clear cell", "kidney papillary", "liver", "lung adeno",
                "lung  small cell", "lung  squamous", "lymphoma b-cell",
                "lymphoma hodgkin", "medulloblastoma", "melanoma", "myeloma",
                "nasopharyngeal carcinoma", "neuroblastoma", "oesophagus",
                "oral gingivo-buccal squamous", "osteosarcoma", "ovary",
                "pancreas", "paraganglioma", "pilocytic astrocytoma", "prostate",
                "stomach", "thyroid", "urothelial carcinoma", "uterine carcinoma"
                , "uterine carcinosarcoma", "uveal melanoma")
  present_sig <- list(
    c(1, 2, 4, 5, 6, 13, 18), c(1, 2, 5, 13), c(1, 5), c(1, 2, 5, 10, 13),
    c(1, 2, 3, 5, 6, 8, 10, 13, 17, 18, 20, 26, 30), c(1, 2, 5, 6, 10, 13, 26),
    c(1, 5), c(1, 2, 5, 9, 13), c(1, 5, 6, 10), c(1, 5, 11), c(1, 5, 6, 14),
    c(1, 2, 4, 5, 7, 13), c(1, 5, 6), c(1, 5, 6, 27), c(1, 2, 5, 13),
    c(1, 4, 5, 6, 12, 16, 17, 22, 23, 24), c(1, 2, 4, 5, 6, 13, 17),
    c(1, 4, 5, 15), c(1, 2, 4, 5, 13), c(1, 2, 5, 9, 13, 17), c(1, 5, 25),
    c(1, 5, 8), c(1, 5, 7, 11, 17), c(1, 2, 5, 13), c(1, 2, 5, 6, 13),
    c(1, 5, 18), c(1, 2, 4, 5, 6, 13, 17), c(1, 2, 5, 7, 13, 29),
    c(1, 2, 5, 6, 13, 30), c(1, 5), c(1, 2, 3, 5, 6, 13), c(1, 5), c(1, 5, 19),
    c(1, 5, 6), c(1, 2, 5, 13, 15, 17, 18, 20, 21, 26, 28), c(1, 2, 5, 13),
    c(1, 2, 5, 13, 22), c(1, 2, 5, 6, 10, 13, 14, 26), c(1, 2, 5, 6, 10, 13),
    c(1, 5, 6)
  )
  partial <- grep(tumor_type, subtypes)
  for (i in seq_len(length(partial))) {
    print(subtypes[partial[i]])
    print(present_sig[[partial[i]]])
  }
}
