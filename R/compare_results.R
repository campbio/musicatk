NULL


sig_compare <- function(sig1, sig2, metric = c("cosine", "jsd"),
                        threshold = 0.9) {
  metric <- match.arg(metric)

  sig1_names <- colnames(sig1)
  sig2_names <- colnames(sig2)
  if (nrow(sig1) != nrow(sig2)) {
    stop("Signatures must have the same number of motifs.")
  }
  if (!is.null(rownames(sig1)) && !is.null(rownames(sig2)) &&
        !all(rownames(sig1) == rownames(sig2))) {
    warning(
      "The names of the motifs in signature matrix one do not equal the ",
      "names of the motifs in signature matrix two."
    )
  }
  if (metric == "jsd") {
    matches <- .jsd(sig1, sig2)
    metric_name <- "1_minus_jsd"
  } else {
    matches <- .cosine(sig1, sig2)
    metric_name <- "cosine"
  }

  comparison <- NULL
  for (row in seq_len(nrow(matches))) {
    line <- which(matches[row, ] > threshold)
    if (length(line) > 0) {
      for (match in line) {
        comparison <- rbind(comparison, c(
          matches[row, match], row, match,
          sig1_names[row], sig2_names[match]
        ))
      }
    }
  }
  if (is.null(comparison)) {
    stop("No matches found, try lowering threshold.")
  }
  comparison <- data.frame(comparison, stringsAsFactors = FALSE)
  colnames(comparison) <- c(
    metric_name, "x_sig_index", "y_sig_index",
    "x_sig_name", "y_sig_name"
  )
  comparison[[metric_name]] <- as.numeric(comparison[[metric_name]])
  comparison$x_sig_index <- as.numeric(comparison$x_sig_index)
  comparison$y_sig_index <- as.numeric(comparison$y_sig_index)
  comparison <- comparison[order(comparison[[metric_name]],
                                 decreasing = TRUE), ]
  return(comparison)
}

#' Compare two result files to find similar signatures
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param model_id The name of the first model to compare.
#' @param other_model_id The name of the second model to compare.
#' @param modality Modality of results being compared. Default \code{"SBS96"}.
#' @param result_name Name of the result list entry for the first result to
#' compare. Default \code{"result"}.
#' @param other_musica A second \code{\linkS4class{musica}} object. If null, the
#' provided musica object is used twice. Default \code{NULL}.
#' @param other_result_name Name of the result list entry for the second result
#' to compare. Default \code{"result"}.
#' @param threshold threshold for similarity
#' @param metric One of \code{"cosine"} for cosine similarity or \code{"jsd"}
#' for 1 minus the Jensen-Shannon Divergence. Default \code{"cosine"}.
#' @param result_rename title for plot of first result signatures
#' @param other_result_rename title for plot of second result signatures
#' @param decimals Specifies rounding for similarity metric displayed. Default
#' \code{2}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{FALSE}.
#' @return Returns the comparisons
#' @examples
#' data(res)
#' compare_results(res,
#'   model_id = "res", other_model_id = "res",
#'   modality = "SBS96", threshold = 0.8
#' )
#' @export
compare_results <- function(musica, model_id, other_model_id,
                            modality = "SBS96", result_name = "result",
                            other_musica = NULL, other_result_name = "result",
                            threshold = 0.9, metric = "cosine", result_rename =
                              deparse(substitute(model_id)),
                            other_result_rename =
                              deparse(substitute(other_model_id)),
                            decimals = 2, same_scale = FALSE) {
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))) {
    stop(
      result_name, " does not exist in the result_list. Current names are: ",
      paste(names(result_list(musica)), collapse = ", ")
    )
  }

  # check if valid modality
  if (!(modality %in%
          names(get_result_list_entry(musica, result_name)@modality))) {
    stop(
      modality, " is not a valid modality. Current modalities are: ",
      paste(names(get_result_list_entry(musica, result_name)@modality),
            collapse = ", ")
    )
  }

  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))) {
    stop(
      model_id, " is not a valid model_id. Current model names are: ",
      paste(names(get_modality(musica, result_name, modality)), collapse = ", ")
    )
  }

  # if no second musica object provided, use the given one
  if (is.null(other_musica)) {
    other_musica <- musica
  }

  # check if valid other_result_name
  if (!(other_result_name %in% names(result_list(other_musica)))) {
    stop(
      other_result_name,
      " does not exist in the result_list. Current names are: ",
      paste(names(result_list(other_musica)), collapse = ", ")
    )
  }

  # check if valid modality
  if (!(modality %in%
          names(get_result_list_entry(other_musica, other_result_name)
                @modality))) {
    stop(
      modality, " is not a valid modality. Current modalities are: ",
      paste(names(get_result_list_entry(other_musica, other_result_name)
                  @modality), collapse = ", ")
    )
  }

  # check if valid model_id
  if (!(other_model_id %in% names(get_modality(other_musica, other_result_name,
                                               modality)))) {
    stop(
      other_model_id, " is not a valid model_id. Current model names are: ",
      paste(names(get_modality(other_musica, other_result_name, modality)),
            collapse = ", ")
    )
  }

  result <- get_model(musica, result_name, modality, model_id)
  other_result <- get_model(other_musica, other_result_name, modality,
                            other_model_id)

  signatures <- signatures(result)
  comparison <- sig_compare(
    sig1 = signatures, sig2 = signatures(other_result),
    threshold = threshold, metric = metric
  )
  result_subset <- methods::new("result_model",
    signatures =
      signatures(result)[, comparison$x_sig_index, drop = FALSE],
    exposures = matrix(),
    num_signatures = dim(signatures)[2],
    modality = modality, model_id = "result_subset"
  )
  other_subset <- methods::new("result_model",
    signatures =
      signatures(other_result)[, comparison$y_sig_index, drop = FALSE],
    exposures = matrix(), modality = modality,
    model_id = "result_subset"
  )

  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in seq_len(dim(comparison)[1])) {
    result_subset_maxes <- c(result_subset_maxes,
                             max(signatures(result_subset)[, index]))
  }
  for (index in seq_len(dim(comparison)[1])) {
    other_subset_maxes <- c(other_subset_maxes,
                            max(signatures(other_subset)[, index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100

  if (same_scale == TRUE) {
    maxes <- rep(max(maxes), length(maxes))
  }

  .plot_compare_result_signatures(result_subset, other_subset, comparison,
                                  musica, res1_name = result_rename,
                                  res2_name = other_result_rename,
                                  decimals = decimals, same_scale = same_scale,
                                  maxes = maxes)
  return(comparison)
}

#' Compare a result object to COSMIC V3 Signatures; Select exome or genome for
#' SBS and only genome for DBS or Indel classes
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param model_id The name of the model containing the signatures to compare.
#' @param sample_type exome (SBS only) or genome
#' @param modality Compare to SBS, DBS, or Indel. Default \code{"SBS96"}
#' @param result_name Name of the result list entry. Default \code{"result"}.
#' @param threshold threshold for similarity
#' @param metric One of \code{"cosine"} for cosine similarity or \code{"jsd"}
#' for 1 minus the Jensen-Shannon Divergence. Default \code{"cosine"}.
#' @param result_rename title for plot user result signatures
#' @param decimals Specifies rounding for similarity metric displayed. Default
#' \code{2}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{TRUE}.
#' @param genome_build Reference build the COSMIC signatures are normalized to,
#' either \code{"hg19"}/\code{"GRCh37"} or \code{"hg38"}/\code{"GRCh38"}. This
#' should match the build the variants were called against. Case is ignored,
#' and \code{"19"}, \code{"38"}, \code{19} and \code{38} also work. COSMIC
#' publishes one Indel catalogue per version, used for either build. The v3.0
#' exome catalogue is GRCh37 only. Default \code{"hg19"}.
#' @param cosmic_version COSMIC version to compare against: \code{"latest"},
#' \code{"3.0"} or \code{"3.6"} (also given as the number \code{3.6} or as
#' \code{"v3.6"}). \code{"latest"} uses the newest version that
#' has a catalogue for the sample type, modality and build, which is v3.6
#' for every catalogue musicatk ships. Use \code{"3.0"} to reproduce results
#' from musicatk 2.7.15 and earlier, which compared against v3.0 (GRCh37)
#' only. A version without a
#' matching catalogue is an error listing the catalogues musicatk ships;
#' there is no fallback to another version or build. Default \code{"latest"}.
#' @return Returns the comparisons. Its \code{"cosmic_catalogue"} attribute is
#' the name of the COSMIC data object compared against, such as
#' \code{"cosmic_v3_6_sbs_sigs_grch37"}, and a message names its version,
#' modality, sample type and build.
#' @seealso \code{\link{cosmic_v3_6_sbs_sigs_grch37}} and the other COSMIC
#' signature objects.
#' @examples
#' data(res)
#' # The example data are hg38, so compare against the GRCh38 catalogue
#' compare_cosmic_v3(res,
#'   model_id = "res", modality = "SBS96",
#'   sample_type = "genome", genome_build = "hg38", threshold = 0.8
#' )
#' @export
compare_cosmic_v3 <- function(musica, model_id, sample_type, modality = "SBS96",
                              result_name = "result", metric = "cosine",
                              threshold = 0.9,
                              result_rename = deparse(substitute(model_id)),
                              decimals = 2, same_scale = FALSE,
                              genome_build = "hg19",
                              cosmic_version = "latest") {
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))) {
    stop(
      result_name, " does not exist in the result_list. Current names are: ",
      paste(names(result_list(musica)), collapse = ", ")
    )
  }

  # check if valid modality
  if (!(modality %in%
          names(get_result_list_entry(musica, result_name)@modality))) {
    stop(
      modality, " is not a valid modality. Current modalities are: ",
      paste(names(get_result_list_entry(musica, result_name)@modality),
            collapse = ", ")
    )
  }

  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))) {
    stop(
      model_id, " is not a valid model_id. Current model names are: ",
      paste(names(get_modality(musica, result_name, modality)), collapse = ", ")
    )
  }

  if (!(length(sample_type) == 1 && sample_type %in% c("exome", "genome"))) {
    stop("Sample type must be exome or genome")
  }
  cosmic_modality <- .cosmic_modality(modality)
  if (sample_type == "exome" && !identical(cosmic_modality, "SBS96")) {
    stop(paste("Only SBS class is available for whole-exome, please choose",
      " `genome` for DBS or Indel",
      sep = ""
    ))
  }
  if (is.na(cosmic_modality)) {
    stop("Only SBS, DBS, and Indel classes are supported")
  }
  catalogue <- .select_cosmic_catalogue(
    sample_type, cosmic_modality, .cosmic_build(genome_build),
    .cosmic_version(cosmic_version)
  )
  message(catalogue$message)
  cosmic_title <- paste0("COSMIC Signatures (v", catalogue$version, ")")
  cosmic_res <- getExportedValue("musicatk", catalogue$object)

  signatures <- signatures(get_model(musica, result_name, modality, model_id))
  comparison <- sig_compare(
    sig1 = signatures, sig2 = signatures(cosmic_res),
    threshold = threshold, metric = metric
  )
  result_subset <- methods::new(
    "result_model",
    signatures = signatures[, comparison$x_sig_index,
      drop = FALSE
    ],
    exposures = matrix(), num_signatures = dim(signatures)[2],
    modality = modality, model_id = "result_subset"
  )
  other_subset <- methods::new("result_model",
    signatures =
      signatures(cosmic_res)[,
        comparison$y_sig_index,
        drop = FALSE
      ],
    exposures = matrix(), modality = modality,
    model_id = "other_subset"
  )

  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in seq_len(dim(comparison)[1])) {
    result_subset_maxes <- c(result_subset_maxes,
                             max(signatures(result_subset)[, index]))
  }
  for (index in seq_len(dim(comparison)[1])) {
    other_subset_maxes <- c(other_subset_maxes,
                            max(signatures(other_subset)[, index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100

  if (same_scale == TRUE) {
    maxes <- rep(max(maxes), length(maxes))
  }

  .plot_compare_result_signatures(result_subset, other_subset, comparison,
                                  musica, res1_name = result_rename,
                                  res2_name = cosmic_title,
                                  decimals = decimals, maxes = maxes,
                                  same_scale = same_scale)
  attr(comparison, "cosmic_catalogue") <- catalogue$object
  return(comparison)
}

# The COSMIC catalogues musicatk ships (ADR 0003). An NA build means one
# catalogue for all builds, as COSMIC publishes for Indel. A new COSMIC
# release adds rows here.
.cosmic_catalogues <- function() {
  data.frame(
    version = c(rep("3.0", 4), rep("3.6", 7)),
    modality = c(
      "SBS96", "SBS96", "DBS78", "IND83",
      "SBS96", "SBS96", "SBS96", "SBS96", "DBS78", "DBS78", "IND83"
    ),
    sample_type = c(
      "exome", rep("genome", 3),
      "exome", "exome", rep("genome", 5)
    ),
    build = c(
      "GRCh37", "GRCh37", "GRCh37", NA,
      "GRCh37", "GRCh38", "GRCh37", "GRCh38", "GRCh37", "GRCh38", NA
    ),
    object = c(
      "cosmic_v3_sbs_sigs_exome", "cosmic_v3_sbs_sigs", "cosmic_v3_dbs_sigs",
      "cosmic_v3_indel_sigs", "cosmic_v3_6_sbs_sigs_grch37_exome",
      "cosmic_v3_6_sbs_sigs_grch38_exome", "cosmic_v3_6_sbs_sigs_grch37",
      "cosmic_v3_6_sbs_sigs_grch38", "cosmic_v3_6_dbs_sigs_grch37",
      "cosmic_v3_6_dbs_sigs_grch38", "cosmic_v3_6_indel_sigs"
    ),
    stringsAsFactors = FALSE
  )
}

# The catalogue modality for a modality name, or NA if COSMIC has none
.cosmic_modality <- function(modality) {
  if (modality %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
    return("SBS96")
  }
  if (modality %in% c("DBS", "DBS78", "dbs", "doublet")) {
    return("DBS78")
  }
  if (modality %in% c("INDEL", "Indel", "indel", "ind", "IND", "ID", "IND83")) {
    return("IND83")
  }
  NA_character_
}

# Normalize a build label the same way select_genome() does
.cosmic_build <- function(genome_build) {
  if (is.numeric(genome_build)) {
    genome_build <- as.character(genome_build)
  }
  build <- if (is.character(genome_build) && length(genome_build) == 1) {
    tolower(genome_build)
  } else {
    ""
  }
  if (build %in% c("hg19", "19", "grch37")) {
    return("GRCh37")
  }
  if (build %in% c("hg38", "38", "grch38")) {
    return("GRCh38")
  }
  stop("genome_build must be one of hg19 (GRCh37) or hg38 (GRCh38)")
}

# "latest" or a shipped version such as "3.6"; also takes 3.6 and "v3.6"
.cosmic_version <- function(cosmic_version) {
  versions <- unique(.cosmic_catalogues()$version)
  if (is.numeric(cosmic_version) && length(cosmic_version) == 1 &&
        !is.na(cosmic_version)) {
    cosmic_version <- sprintf("%.1f", cosmic_version)
  }
  if (is.character(cosmic_version) && length(cosmic_version) == 1 &&
        !is.na(cosmic_version)) {
    version <- sub("^v", "", tolower(cosmic_version))
    if (version == "latest" || version %in% versions) {
      return(version)
    }
  }
  stop(
    "cosmic_version must be \"latest\" or one of: ",
    paste(versions, collapse = ", ")
  )
}

# Pick the catalogue for a sample type, modality, build and version. "latest"
# is the newest version with a matching catalogue; there is no fallback to
# another build or version.
.select_cosmic_catalogue <- function(sample_type, modality, build, version) {
  catalogues <- .cosmic_catalogues()
  offered <- catalogues[catalogues$sample_type == sample_type &
                          catalogues$modality == modality, ]
  matches <- offered[is.na(offered$build) | offered$build == build, ]
  if (version != "latest") {
    matches <- matches[matches$version == version, ]
  }
  if (nrow(matches) == 0) {
    available <- paste0(
      "v", offered$version, " ",
      ifelse(is.na(offered$build), "any build", offered$build)
    )
    requested <- if (version == "latest") "" else paste0("v", version, " ")
    stop(
      "musicatk has no COSMIC ", requested, modality, " ", sample_type,
      " catalogue for ", build, ". Available for ", modality, " ",
      sample_type, ": ", paste(available, collapse = ", ")
    )
  }
  newest <- order(numeric_version(matches$version), decreasing = TRUE)[1]
  chosen <- matches[newest, ]
  detail <- if (is.na(chosen$build)) {
    "one catalogue for all builds"
  } else {
    chosen$build
  }
  text <- paste0(
    "Comparing to COSMIC v", chosen$version, " ", modality, " ",
    sample_type, " signatures (", detail, ")."
  )
  list(object = chosen$object, version = chosen$version, message = text)
}

#' Compare a result object to COSMIC V2 SBS Signatures (combination whole-exome
#' and whole-genome)
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param model_id The name of the model containing the signatures to compare.
#' @param modality Compare to SBS, DBS, or Indel. Default \code{"SBS96"}
#' @param result_name Name of the result list entry. Default \code{"result"}.
#' @param threshold threshold for similarity
#' @param metric One of \code{"cosine"} for cosine similarity or \code{"jsd"}
#' for 1 minus the Jensen-Shannon Divergence. Default \code{"cosine"}.
#' @param result_rename title for plot user result signatures
#' @param decimals Specifies rounding for similarity metric displayed. Default
#' \code{2}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{TRUE}.
#' @return Returns the comparisons
#' @examples
#' data(res)
#' compare_cosmic_v2(res, model_id = "res", threshold = 0.7)
#' @export
compare_cosmic_v2 <- function(musica, model_id, modality = "SBS96",
                              result_name = "result", metric = "cosine",
                              threshold = 0.9,
                              result_rename = deparse(substitute(result)),
                              decimals = 2, same_scale = FALSE) {
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))) {
    stop(
      result_name, " does not exist in the result_list. Current names are: ",
      paste(names(result_list(musica)), collapse = ", ")
    )
  }

  # check if valid modality
  if (!(modality %in%
          names(get_result_list_entry(musica, result_name)@modality))) {
    stop(
      modality, " is not a valid modality. Current modalities are: ",
      paste(names(get_result_list_entry(musica, result_name)@modality),
            collapse = ", ")
    )
  }

  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))) {
    stop(
      model_id, " is not a valid model_id. Current model names are: ",
      paste(names(get_modality(musica, result_name, modality)), collapse = ", ")
    )
  }

  signatures <- signatures(get_model(musica, result_name, modality, model_id))
  comparison <- sig_compare(
    sig1 = signatures,
    sig2 = signatures(musicatk::cosmic_v2_sigs),
    threshold = threshold, metric = metric
  )
  result_subset <- methods::new(
    "result_model",
    signatures = signatures[, comparison$x_sig_index,
      drop = FALSE
    ],
    exposures = matrix(), num_signatures = dim(signatures)[2],
    modality = modality, model_id = "result_subset"
  )
  other_subset <- methods::new("result_model",
    signatures =
      signatures(musicatk::cosmic_v2_sigs)[,
        comparison$y_sig_index,
        drop = FALSE
      ],
    exposures = matrix(), modality = modality,
    model_id = "other_subset"
  )


  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in seq_len(dim(comparison)[1])) {
    result_subset_maxes <- c(result_subset_maxes,
                             max(signatures(result_subset)[, index]))
  }
  for (index in seq_len(dim(comparison)[1])) {
    other_subset_maxes <- c(other_subset_maxes,
                            max(signatures(other_subset)[, index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100

  if (same_scale == TRUE) {
    maxes <- rep(max(maxes), length(maxes))
  }

  .plot_compare_result_signatures(result_subset, other_subset, comparison,
                                  musica, res1_name = result_rename,
                                  res2_name = "COSMIC Signatures (V2)",
                                  decimals = decimals, maxes = maxes,
                                  same_scale = same_scale)

  return(comparison)
}


#' Input a cancer subtype to return a list of related COSMIC signatures
#'
#' @param tumor_type Cancer subtype to view related signatures
#' @return Returns signatures related to a partial string match
#' @examples cosmic_v2_subtype_map("lung")
#' @export
cosmic_v2_subtype_map <- function(tumor_type) {
  subtypes <- c(
    "adrenocortical carcinoma", "all", "aml", "bladder", "breast",
    "cervix", "chondrosarcoma", "cll", "colorectum", "glioblastoma",
    "glioma low grade", "head and neck", "kidney chromophobe",
    "kidney clear cell", "kidney papillary", "liver", "lung adeno",
    "lung  small cell", "lung  squamous", "lymphoma b-cell",
    "lymphoma hodgkin", "medulloblastoma", "melanoma", "myeloma",
    "nasopharyngeal carcinoma", "neuroblastoma", "oesophagus",
    "oral gingivo-buccal squamous", "osteosarcoma", "ovary",
    "pancreas", "paraganglioma", "pilocytic astrocytoma",
    "prostate", "stomach", "thyroid", "urothelial carcinoma",
    "uterine carcinoma", "uterine carcinosarcoma", "uveal melanoma"
  )
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
  for (i in seq_along(partial)) {
    message(subtypes[partial[i]])
    message(present_sig[[partial[i]]])
  }
}


.plot_compare_result_signatures <- function(res1, res2, comparison, musica,
                                            res1_name = "", res2_name = "",
                                            decimals, maxes, same_scale) {
  if (colnames(comparison)[1] == "cosine") {
    annotations <- paste("CS = ", round(comparison$cosine, decimals), sep = "")
  }

  if (colnames(comparison)[1] == "1_minus_jsd") {
    annotations <- paste("1 - JSD = ",
                         round(comparison$"1_minus_jsd", decimals), sep = "")
  }

  res1_plot <- .plot_result_model_signatures(res1, musica,
    annotation = annotations,
    y_max = maxes, same_scale = same_scale
  )
  res1_plot <- ggpubr::annotate_figure(res1_plot,
    top = ggpubr::text_grob(res1_name, face = "bold")
  )

  res2_plot <- .plot_result_model_signatures(res2, musica,
    y_max = maxes,
    show_y_labels = FALSE, same_scale = same_scale
  )
  res2_plot <- ggpubr::annotate_figure(res2_plot,
    top = ggpubr::
      text_grob(res2_name, face = "bold")
  )

  layout <- matrix(seq(2), ncol = 2, nrow = 9, byrow = TRUE)
  layout <- rbind(layout, c(3, 3))
  g <- gridExtra::grid.arrange(res1_plot, res2_plot,
    layout_matrix = layout, ncol = 2,
    widths = c(0.45, 0.4)
  )

  return(g)
}
