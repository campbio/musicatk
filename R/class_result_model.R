#' Object that contains results for a single model
#' 
#' @slot signatures A matrix of signatures by mutational motifs
#' @slot exposures A matrix of samples by signature weights
#' @slot num_signatures Number of signatures in the model
#' @slot other_parameters Parameters relevant to the model
#' @slot credible_intervals Credible intervals for parameters
#' @slot metrics Performance metrics for the model
#' @slot umap List of umap data.frames for plotting and analysis
#' @slot model_id Model identifier
#' @export
#' @exportClass result_model

setClass(
  "result_model",
  slots = list(
    signatures = "matrix",
    exposures = "matrix",
    num_signatures = "numeric",
    other_parameters = "SimpleList",
    credible_intervals = "SimpleList",
    metrics = "SimpleList",
    umap = "matrix",
    model_id = "character"
  )
)
