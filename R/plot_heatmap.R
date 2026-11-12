#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @title Plot heatmaps using the exposures matrix
#' @description The exposures for different signatures can be
#' visualized using a heatmap with this function.
#' Heatmaps make it easier to visualize the data by
#' representing the magnitude of exposure values
#' as color in 2-dimensions. The variation in color
#' intensity can help see if the exposures are clustered
#' or how they vary over space. Exposures can be
#' normalized by providing the \code{proportional} argument.
#' Column annotations can also be seen by passing the \code{col_annot} argument.
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' discovery or prediction.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing desired model.
#' Default \code{"result"}.
#' @param proportional If \code{TRUE}, then the exposures will be normalized
#' to between 0 and 1 by dividing by the total number of counts for each sample.
#' Default \code{FALSE}.
#' @param annotation Column annotations extracted from the musicatk object
#' @param show_column_names Boolean check. If \code{True}, column names are
#' shown. Otherwise, they aren't.
#' Default \code{FALSE}
#' @param show_row_names Boolean check. If \code{True}, row names are shown.
#' Otherwise, they aren't.
#' Default \code{FALSE}
#' @param scale Boolean check. If \code{True}, each signature (row) is
#' scaled to a z-score across samples; a signature with the same value in
#' every sample gets 0. Proportions and z-scores use all samples, before
#' \code{subset_tumor} is applied. Otherwise, values aren't scaled.
#' Default \code{TRUE}
#' @param annotation Users have the option of plotting the exposure matrix
#' based on their given annotation like
#' Tumor_Subtypes or age. Error given if the user given annotation
#' doesn't exist in the res_annot annotation object.
#' @param subset_tumor One or more annotation values, such as tumor types.
#' Only samples whose value in any sample annotation column exactly matches
#' one of them are plotted. Unless \code{annotation} is given, the columns
#' where the values were found are shown as the heatmap annotation. Error
#' given if a value is not in any sample annotation. Default \code{NULL}.
#' @param subset_signatures Users can specify certain signatures on
#'  which they want to subset the exposure matrix
#' plotting the heatmap.
#' @param ... Ellipsis used for passing any arguments directly
#' to the ComplexHeatmap's heatmap function.
#' @return Generates a heatmap for using the exposure matrix.
#' @examples
#' data(res_annot)
#' plot_heatmap(
#'   musica = res_annot, model_name = "res_annot",
#'   proportional = TRUE, scale = TRUE, annotation = "Tumor_Subtypes"
#' )
#' @export
plot_heatmap <- function(musica, model_name,
                         modality = "SBS96",
                         result_name = "result",
                         proportional = FALSE,
                         show_column_names = FALSE,
                         show_row_names = TRUE,
                         scale = TRUE,
                         subset_tumor = NULL,
                         subset_signatures = NULL,
                         annotation = NULL,
                         ...) {
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

  # check if valid model_name
  if (!(model_name %in% names(get_modality(musica, result_name, modality)))) {
    stop(
      model_name, " is not a valid model_name. Current model names are: ",
      paste(names(get_modality(musica, result_name, modality)), collapse = ", ")
    )
  }

  # Get result object from musica object
  result <- get_model(musica,
    result = result_name,
    modality = modality,
    model = model_name
  )

  exp <- exposures(result) # extracting the exposures matrix
  if (isTRUE(proportional)) {
    exp <- .proportional_exposures(exp) # normalizing
  }
  if (isTRUE(scale)) {
    # z-score each signature (row) across samples; a constant signature
    # has no spread, so it gets 0
    exp <- t(scale(t(exp)))
    exp[is.nan(exp)] <- 0
  }
  annot <- samp_annot(musica) # Extracting annotations from musicatk object
  # checking if annotation argument correct
  if (!is.null(annotation)) {
    if (any(!annotation %in% colnames(annot))) {
      stop("The given annotations are not present in the data")
    }
  }
  if (!is.null(subset_signatures)) {
    if (any(!subset_signatures %in% rownames(exp))) {
      stop("The given signatures are not present in the data")
    }
    exp <- exp[rownames(exp) %in% subset_signatures, , drop = FALSE]
  }
  if (!is.null(subset_tumor)) {
    matched <- .samples_with_annotation(annot, subset_tumor)
    exp <- exp[, colnames(exp) %in% matched$samples, drop = FALSE]
    if (ncol(exp) == 0) {
      stop(
        "No samples in model '", model_name, "' match 'subset_tumor': ",
        paste(subset_tumor, collapse = ", ")
      )
    }
    # Without 'annotation', show the columns the values were found in
    if (is.null(annotation)) {
      annotation <- matched$columns
    }
  }
  heatmap <- NULL # initializing an empty heatmap annotation variable
  if (!is.null(annotation)) {
    # Match each annotation to the exposure columns by sample name
    annot <- lapply(annotation, function(column) {
      .annot_for_samples(musica, column, colnames(exp))
    })
    names(annot) <- annotation
    annot <- data.frame(annot, check.names = FALSE)
    heatmap <- ComplexHeatmap::HeatmapAnnotation(df = annot)
  }
  # If/else conditions to check if annotation object available
  if (is.null(heatmap)) {
    if (scale == FALSE) {
      ComplexHeatmap::Heatmap(exp, name = "exposures",
                              show_column_names = show_column_names,
                              show_row_names = show_row_names,
                              col = c("blue", "green", "yellow", "orange",
                                      "red"), ...)
    } else {
      ComplexHeatmap::Heatmap(exp, name = "exposures",
                              show_column_names = show_column_names,
                              show_row_names = show_row_names, ...)
    }
  } else if (!is.null(heatmap)) {
    if (scale == FALSE) {
      ComplexHeatmap::Heatmap(exp, name = "exposures",
                              top_annotation = heatmap,
                              show_column_names = show_column_names,
                              show_row_names = show_row_names,
                              col = c("blue", "green", "yellow", "orange",
                                      "red"), ...)
    } else {
      ComplexHeatmap::Heatmap(exp, name = "exposures", top_annotation = heatmap,
                              show_column_names = show_column_names,
                              show_row_names = show_row_names, ...)
    }
  }
}

# Samples whose value in any sample annotation column (other than Samples)
# is exactly one of 'values', and the columns where the values were found.
# Stops, naming them, if some values are in no column.
.samples_with_annotation <- function(annot, values) {
  values <- as.character(values)
  columns <- setdiff(colnames(annot), "Samples")
  column_values <- lapply(columns, function(column) {
    as.character(annot[[column]])
  })
  matches <- lapply(column_values, function(x) x %in% values)
  missing <- setdiff(values, unlist(column_values))
  if (length(missing) > 0) {
    stop(
      "These 'subset_tumor' values are not in any sample annotation: ",
      paste(missing, collapse = ", ")
    )
  }
  keep <- Reduce(`|`, matches)
  list(
    samples = as.character(annot$Samples[keep]),
    columns = columns[vapply(matches, any, logical(1))]
  )
}
