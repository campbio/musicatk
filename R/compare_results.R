NULL


sig_compare <- function(sig1, sig2, metric = c("cosine", "jsd"),
                        threshold=0.9) {
  metric <- match.arg(metric)
  
  sig1_names <- colnames(sig1)
  sig2_names <- colnames(sig2)
  if (nrow(sig1) != nrow(sig2)) {
    stop("Signatures must have the same number of motifs.")
  }
  if (!is.null(rownames(sig1)) && !is.null(rownames(sig2)) &&
     !all(rownames(sig1) == rownames(sig2))) {
    warning("The names of the motifs in signature matrix one do not equal the ",
    "names of the motifs in signature matrix two.")
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
        comparison <- rbind(comparison, c(matches[row, match], row, match,
                                          sig1_names[row], sig2_names[match]))
      }
    }
  }
  if (is.null(comparison)) {
    stop("No matches found, try lowering threshold.")
  }
  comparison <- data.frame(comparison, stringsAsFactors = FALSE)
  colnames(comparison) <- c(metric_name, "x_sig_index", "y_sig_index", 
                            "x_sig_name", "y_sig_name")
  comparison[[metric_name]] <- as.numeric(comparison[[metric_name]])
  comparison$x_sig_index <- as.numeric(comparison$x_sig_index)
  comparison$y_sig_index <- as.numeric(comparison$y_sig_index)
  comparison <- comparison[order(comparison[[metric_name]], decreasing = TRUE), 
                           ]
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
#' compare_results(res, model_id = "res", other_model_id = "res",
#' modality = "SBS96", threshold = 0.8)
#' @export
compare_results <- function(musica, model_id, other_model_id, 
                            modality = "SBS96", result_name = "result", 
                            other_musica = NULL, other_result_name = "result",
                            threshold = 0.9, metric = "cosine", result_rename =
                            deparse(substitute(model_id)), other_result_rename =
                            deparse(substitute(other_model_id)),
                            decimals = 2, same_scale = FALSE) {
  
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))){
    stop(result_name, " does not exist in the result_list. Current names are: ",
         paste(names(result_list(musica)), collapse = ", "))
  }
  
  # check if valid modality
  if (!(modality %in% names(get_result_list_entry(musica, result_name)@modality))){
    stop(modality, " is not a valid modality. Current modalities are: ", 
         paste(names(get_result_list_entry(musica, result_name)@modality), collapse = ", "))
  }
  
  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))){
    stop(model_id, " is not a valid model_id. Current model names are: ",
         paste(names(get_modality(musica, result_name, modality)), collapse = ", "))
  }
  
  # if no second musica object provided, use the given one
  if (is.null(other_musica)){
    other_musica <- musica
  }
  
  # check if valid other_result_name
  if (!(other_result_name %in% names(result_list(other_musica)))){
    stop(other_result_name, " does not exist in the result_list. Current names are: ",
         paste(names(result_list(other_musica)), collapse = ", "))
  }
  
  # check if valid modality
  if (!(modality %in% names(get_result_list_entry(other_musica, other_result_name)@modality))){
    stop(modality, " is not a valid modality. Current modalities are: ", 
         paste(names(get_result_list_entry(other_musica, other_result_name)@modality), collapse = ", "))
  }
  
  # check if valid model_id
  if (!(other_model_id %in% names(get_modality(other_musica, other_result_name, modality)))){
    stop(other_model_id, " is not a valid model_id. Current model names are: ",
         paste(names(get_modality(other_musica, other_result_name, modality)), collapse = ", "))
  }
  
  result <- get_model(musica, result_name, modality, model_id)
  other_result <- get_model(other_musica, other_result_name, modality, other_model_id)
  
  signatures <- signatures(result)
  comparison <- sig_compare(sig1 = signatures, sig2 = signatures(other_result),
                            threshold = threshold, metric = metric)
  result_subset <- methods::new("result_model",
                                signatures = 
                                  signatures(result)[, comparison$x_sig_index, drop = FALSE], 
                                exposures = matrix(), 
                                num_signatures = dim(signatures)[2],
                                modality = modality, model_id = "result_subset")
  other_subset <- methods::new("result_model",
                               signatures = 
                                 signatures(other_result)[, comparison$y_sig_index, drop = FALSE], 
                               exposures = matrix(), modality = modality, 
                               model_id = "result_subset")
  
  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in 1:dim(comparison)[1]){
    result_subset_maxes <- c(result_subset_maxes, max(signatures(result_subset)[,index]))
  }
  for (index in 1:dim(comparison)[1]){
    other_subset_maxes <- c(other_subset_maxes, max(signatures(other_subset)[,index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100
  
  if (same_scale == TRUE){
    maxes <- rep(max(maxes), length(maxes))
  }

  .plot_compare_result_signatures(result_subset, other_subset, comparison, musica,
                                  res1_name = result_rename,
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
#' @return Returns the comparisons
#' @examples
#' data(res)
#' compare_cosmic_v3(res, model_id = "res", modality = "SBS96",
#' sample_type = "genome", threshold = 0.8)
#' @export
compare_cosmic_v3 <- function(musica, model_id, sample_type, modality = "SBS96",
                              result_name = "result", metric = "cosine", 
                              threshold = 0.9,
                              result_rename = deparse(substitute(model_id)),
                              decimals = 2, same_scale = FALSE) {
  
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))){
    stop(result_name, " does not exist in the result_list. Current names are: ",
         paste(names(result_list(musica)), collapse = ", "))
  }
  
  # check if valid modality
  if (!(modality %in% names(get_result_list_entry(musica, result_name)@modality))){
    stop(modality, " is not a valid modality. Current modalities are: ", 
         paste(names(get_result_list_entry(musica, result_name)@modality), collapse = ", "))
  }
  
  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))){
    stop(model_id, " is not a valid model_id. Current model names are: ",
         paste(names(get_modality(musica, result_name, modality)), collapse = ", "))
  }
  
  if (sample_type == "exome") {
    if (modality %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
      cosmic_res <- musicatk::cosmic_v3_sbs_sigs_exome
    } else {
      stop(paste("Only SBS class is available for whole-exome, please choose",
                 " `genome` for DBS or Indel", sep = ""))
    }
  } else if (sample_type == "genome") {
    if (modality %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
      cosmic_res <- musicatk::cosmic_v3_sbs_sigs
    } else if (modality %in% c("DBS", "dbs", "doublet")) {
      cosmic_res <- musicatk::cosmic_v3_dbs_sigs
    } else if (modality %in% c("INDEL", "Indel", "indel", "ind", "IND",
                                    "ID")) {
      cosmic_res <- musicatk::cosmic_v3_indel_sigs
    } else {
      stop("Only SBS, DBS, and Indel classes are supported")
    }
  } else {
    stop("Sample type must be exome or genome")
  }
  
  signatures <- signatures(get_model(musica, result_name, modality, model_id))
  comparison <- sig_compare(sig1 = signatures, sig2 = signatures(cosmic_res),
                            threshold = threshold, metric = metric)
  result_subset <- methods::new(
    "result_model", signatures = signatures[, comparison$x_sig_index,
                                                     drop = FALSE],
    exposures = matrix(), num_signatures = dim(signatures)[2],
    modality = modality, model_id = "result_subset")
  other_subset <- methods::new("result_model", signatures =
                                 signatures(cosmic_res)[, 
                                                        comparison$y_sig_index, 
                                                        drop = FALSE],
                               exposures = matrix(), modality = modality, 
                               model_id = "other_subset")
  
  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in 1:dim(comparison)[1]){
    result_subset_maxes <- c(result_subset_maxes, max(signatures(result_subset)[,index]))
  }
  for (index in 1:dim(comparison)[1]){
    other_subset_maxes <- c(other_subset_maxes, max(signatures(other_subset)[,index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100
  
  if (same_scale == TRUE){
    maxes <- rep(max(maxes), length(maxes))
  }
  
  .plot_compare_result_signatures(result_subset, other_subset, comparison, musica,
                                  res1_name = result_rename,
                                  res2_name = "COSMIC Signatures (V3)",
                                  decimals = decimals, maxes = maxes, 
                                  same_scale = same_scale)
  return(comparison)
  
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
  if (!(result_name %in% names(result_list(musica)))){
    stop(result_name, " does not exist in the result_list. Current names are: ",
         paste(names(result_list(musica)), collapse = ", "))
  }
  
  # check if valid modality
  if (!(modality %in% names(get_result_list_entry(musica, result_name)@modality))){
    stop(modality, " is not a valid modality. Current modalities are: ", 
         paste(names(get_result_list_entry(musica, result_name)@modality), collapse = ", "))
  }
  
  # check if valid model_id
  if (!(model_id %in% names(get_modality(musica, result_name, modality)))){
    stop(model_id, " is not a valid model_id. Current model names are: ",
         paste(names(get_modality(musica, result_name, modality)), collapse = ", "))
  }
  
  signatures <- signatures(get_model(musica, result_name, modality, model_id))
  comparison <- sig_compare(sig1 = signatures, 
                            sig2 = signatures(musicatk::cosmic_v2_sigs),
                            threshold = threshold, metric = metric)
  result_subset <- methods::new(
    "result_model", signatures = signatures[, comparison$x_sig_index,
                                            drop = FALSE],
    exposures = matrix(), num_signatures = dim(signatures)[2],
    modality = modality, model_id = "result_subset")
  other_subset <- methods::new("result_model", signatures =
                                 signatures(musicatk::cosmic_v2_sigs)[, 
                                                        comparison$y_sig_index, 
                                                        drop = FALSE],
                               exposures = matrix(), modality = modality, 
                               model_id = "other_subset")
  
  
  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in 1:dim(comparison)[1]){
    result_subset_maxes <- c(result_subset_maxes, max(signatures(result_subset)[,index]))
  }
  for (index in 1:dim(comparison)[1]){
    other_subset_maxes <- c(other_subset_maxes, max(signatures(other_subset)[,index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100
  
  if (same_scale == TRUE){
    maxes <- rep(max(maxes), length(maxes))
  }
  
  .plot_compare_result_signatures(result_subset, other_subset, comparison, musica,
                                  res1_name = result_rename,
                                  res2_name = "COSMIC Signatures (V2)",
                                  decimals = decimals, maxes = maxes, same_scale = same_scale)
  
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
                "pancreas", "paraganglioma", "pilocytic astrocytoma", 
                "prostate", "stomach", "thyroid", "urothelial carcinoma", 
                "uterine carcinoma", "uterine carcinosarcoma", "uveal melanoma")
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
    message(subtypes[partial[i]])
    message(present_sig[[partial[i]]])
  }
}


.plot_compare_result_signatures <- function(res1, res2, comparison, musica,
                                            res1_name = "", res2_name = "", 
                                            decimals, maxes, same_scale) {
  
  if (colnames(comparison)[1] == "cosine"){
    annotations <- paste("CS = ", round(comparison$cosine, decimals), sep = "")
  }
  
  if (colnames(comparison)[1] == "1_minus_jsd"){
    annotations <- paste("1 - JSD = ", round(comparison$"1_minus_jsd", decimals), sep = "")
  }
  
  res1_plot <- .plot_result_model_signatures(res1, musica, annotation = annotations, 
                               y_max = maxes, same_scale = same_scale) 
  res1_plot <- ggpubr::annotate_figure(res1_plot, 
                               top = ggpubr::text_grob(res1_name, face = "bold"))
  
  res2_plot <- .plot_result_model_signatures(res2, musica, y_max = maxes, 
                               show_y_labels = FALSE, same_scale = same_scale)
  res2_plot <- ggpubr::annotate_figure(res2_plot, 
                                       top = ggpubr::
                                         text_grob(res2_name, face = "bold"))
  
  layout <- matrix(seq(2), ncol = 2, nrow = 9, byrow = TRUE)
  layout <- rbind(layout, c(3, 3))
  g <- gridExtra::grid.arrange(res1_plot, res2_plot,
                               layout_matrix = layout, ncol = 2,
                               widths = c(0.45, 0.4)) 
  
  return(g)
  
}