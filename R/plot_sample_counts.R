#' @import ggplot2 
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' @import tidyverse 
NULL

#' Plot distribution of sample counts
#' 
#' Displays the proportion of counts for each mutation type across one
#' or more samples.
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param sample_names Names of the samples to plot.
#' @param modality Name of table used for plotting counts. Default \code{"SBS96"}.
#' @param text_size Size of axis text. Default \code{10}.
#' @param show_x_labels If \code{TRUE}, the labels for the mutation types
#' on the x-axis will be shown. Default \code{TRUE}.
#' @param show_y_labels If \code{TRUE}, the y-axis ticks and labels will be 
#' shown. Default \code{TRUE}.
#' @param same_scale If \code{TRUE}, the scale of the y-axis for each
#' sample will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each sample. Default \code{TRUE}.
#' @param annotation Vector of annotations to be displayed in the top right
#' corner of each sample. Vector length must be equivalent to the number of
#' samples. Default \code{NULL}.
#' @return Generates a ggplot object
#' @examples
#' data(musica_sbs96)
#' plot_sample_counts(musica_sbs96, sample_names = 
#' sample_names(musica_sbs96)[1])
#' @export
plot_sample_counts <- function(musica, sample_names, modality = "SBS96", 
                               text_size = 10,
                               show_x_labels = TRUE, show_y_labels = TRUE,
                               same_scale = TRUE, annotation = NULL) {
    
  # check if valid modality
  if (!(modality %in% names(extract_count_tables(musica)))){
    stop(modality, " is not a valid modality. Current modalities are: ", 
         paste(names(extract_count_tables(musica)), collapse = ", "))
  }
  
  # Extract counts for specific samples
  tab <- .extract_count_table(musica, modality)
  ix <- match(sample_names, colnames(tab))
  if (all(is.na(ix))) {
      stop("The values in 'sample_names' did not match any sample IDs in table '",
           modality, "'.")
  }
  else if (anyNA(ix)) {
      warning("The following samples in 'sample_names' were not found  in ",
              "table '", modality, 
              "' and will ", "be exlcuded from the plot: ",
              paste(sample_names[is.na(ix)], collapse = ", "))
      ix <- ix[!is.na(ix)]
  }
  sample_counts <- tab[, ix, drop = FALSE]
  
  result <- methods::new("result_model",
                         signatures = sample_counts, exposures = matrix(),
                         modality = modality)
  g <- .plot_result_model_signatures(result, musica, percent = FALSE, text_size = text_size,
                       show_x_labels = show_x_labels, show_y_labels = show_y_labels,
                       same_scale = same_scale, annotation = annotation)
  return(g)
}

