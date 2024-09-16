#' @import ggplot2 
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' @import tidyverse 
NULL


#' Plot reconstruction error for a sample
#' 
#' Displays the observed distribution of counts for each mutation type, 
#' the distribution of reconstructed counts for each mutation type using 
#' the inferred mutational signatures, and the difference between the 
#' two distributions. 

#' @param musica A \code{\linkS4class{musica}} object.
#' @param result_name Name of the result list entry to use.
#' If null, the first entry is selected. Default \code{NULL}.
#' @param result_modality The modality to use. Must be "SBS96", "DBS78", or
#' "IND83".
#' @param model_name The name of the \code{\linkS4class{result_model}} object to
#' use.
#' @param sample Name of the sample within the
#' \code{\linkS4class{musica_result}} object.
#' @param plotly If \code{TRUE}, the the plot will be made interactive
#' using \code{\link[plotly]{plotly}}. Default \code{FALSE}.
#' @return Generates a ggplot or plotly object
#' @examples
#' data(res)
#' plot_sample_reconstruction_error(res, "TCGA-ER-A197-06A-32D-A197-08")
#' @export
plot_sample_reconstruction_error <- function(musica, result_name = NULL, 
                                             result_modality, model_name, sample,
                                             plotly = FALSE) {
  
  # if no result_name supplied, use first entry
  if (is.null(result_name)){
    result_name <- names(result_list(musica))[1]
    note("Using result_list entry:", result_name)
  }
  
  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))){
    stop(result_name, " does not exist in the result_list.")
  }
  
  # check if valid modality
  if (!(result_modality %in% c("SBS96", "DBS78", "IND83"))){
    stop(result_modality, " is not a valid modality.")
  }
  
  # check if valid model_name
  if (!(model_name %in% names(get_modality(musica, result_name, result_modality)))){
    stop(model_name, " is not a valid model_name.")
  }
  
  signatures <- .extract_count_table(musica, result_modality)[, sample, 
                                                             drop = FALSE]
  sample_name <- colnames(signatures)
  result <- get_model(musica, result_name, modality, model_name)
  reconstructed <- reconstruct_sample(result, sample)
  sigs <- cbind(signatures, reconstructed, signatures - reconstructed)
  colnames(sigs) <- c("Counts", "Reconstructed", "Difference")
  
  recontruct_result <- methods::new("result_model",
                                    signatures = sigs,
                                    exposures = matrix(),
                                    modality = result_modality)
  .plot_result_model_signatures(recontruct_result, musica, percent = FALSE, same_scale = FALSE) +
      ggplot2::ggtitle("Reconstruction error", subtitle = sample_name) + ylab("")
}
