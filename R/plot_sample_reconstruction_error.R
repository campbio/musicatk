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
#' @param sample Name of the sample within the \code{\linkS4class{musica}}
#' object.
#' @param model_id The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78",
#' or "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing desired model.
#' Default \code{"result"}.
#' @param plotly If \code{TRUE}, the the plot will be made interactive
#' using \code{\link[plotly]{plotly}}. Default \code{FALSE}.
#' @return Generates a ggplot or plotly object
#' @examples
#' data(res)
#' plot_sample_reconstruction_error(res, "TCGA-ER-A197-06A-32D-A197-08",
#' model_id = "res")
#' @export
plot_sample_reconstruction_error <- function(musica, sample, model_id,
                                             modality = "SBS96",
                                             result_name = "result",
                                             plotly = FALSE) {
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

  signatures <- .extract_count_table(musica, modality)[, sample,
    drop = FALSE
  ]
  sample_name <- colnames(signatures)
  result <- get_model(musica, result_name, modality, model_id)
  reconstructed <- reconstruct_sample(result, sample)
  sigs <- cbind(signatures, reconstructed, signatures - reconstructed)
  colnames(sigs) <- c("Counts", "Reconstructed", "Difference")

  recontruct_result <- methods::new("result_model",
    signatures = sigs,
    exposures = matrix(),
    modality = modality
  )
  .plot_result_model_signatures(recontruct_result, musica, percent = FALSE,
                                same_scale = FALSE) +
    ggplot2::ggtitle("Reconstruction error", subtitle = sample_name) + ylab("")
}
