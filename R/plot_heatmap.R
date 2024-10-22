#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap HeatmapAnnotation
#' @importFrom dplyr select
#' @importFrom dplyr filter_all
#' @importFrom dplyr any_vars
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
#' @param show_column_names Boolean check. If \code{True}, column names are shown.
#' Otherwise, they aren't.
#' Default \code{FALSE}
#' @param show_row_names Boolean check. If \code{True}, row names are shown.
#' Otherwise, they aren't.
#' Default \code{FALSE}
#' @param scale Boolean check. If \code{True}, values are scaled by z-score.
#' Otherwise, they aren't.
#' Default \code{TRUE}
#' @param annotation Users have the option of plotting the exposure matrix
#' based on their given annotation like
#' Tumor_Subtypes or age. Error given if the user given annotation
#' doesn't exist in the res_annot annotation object.
#' @param subset_tumor Users can specify certain tumor types on which
#' they want to subset the exposure matrix for
#' plotting the heatmap.
#' @param subset_signatures Users can specify certain signatures on
#'  which they want to subset the exposure matrix
#' plotting the heatmap.
#' @param ... Ellipsis used for passing any arguments directly
#' to the ComplexHeatmap's heatmap function.
#' @return Generates a heatmap for using the exposure matrix.
#' @examples
#' data(res_annot)
#' plot_heatmap(musica = res_annot, model_name = "res_annot", 
#' proportional = TRUE, scale = TRUE, annotation = "Tumor_Subtypes")
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
                         ...
                         ){
  
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
  
  # check if valid model_name
  if (!(model_name %in% names(get_modality(musica, result_name, modality)))){
    stop(model_name, " is not a valid model_name. Current model names are: ",
         paste(names(get_modality(musica, result_name, modality)), collapse = ", "))
  }
  
  # Get result object from musica object
  result <- get_model(musica, result = result_name,
                      modality = modality,
                      model = model_name)
  
  exp <- exposures(result) #extracting the exposures matrix
  if (isTRUE(proportional)) {
    exp <- sweep(exp, 2, colSums(exp), FUN = "/") #normalizing
  }
  if (isTRUE(scale)) {
    exp <- scale(exp)
  }
  annot <- samp_annot(musica) #Extracting annotations from musicatk object
  heatmap <- NULL #initializing an empty heatmap annotation variable
  #checking if annotation argument correct
  if (!is.null(annotation)) {
    if (any(!annotation %in% colnames(annot))) { 
      stop("The given annotations are not present in the data") 
    }
    annot <- as.data.frame(annot)
    annot <- annot[, annotation]
    heatmap <- ComplexHeatmap::HeatmapAnnotation(df = annot)
  }
  if (!is.null(subset_signatures)) {
    if (any(!subset_signatures %in% rownames(exp))) { 
      stop("The given signatures are not present in the data") 
    }
    exp <- as.data.frame(exp)
    exp <- subset(exp, rownames(exp) %in% subset_signatures)
    exp <- as.matrix(exp)
  }
  if (!is.null(subset_tumor)) {
    . <- NULL  # Gets rid of NOTE: no visible binding for global variable ‘.’
    annot <- samp_annot(musica)
    samps <- annot %>% dplyr::filter_all(dplyr::any_vars(grepl(subset_tumor,.)))
    samps <- as.character(samps$Samples)
    exp <- as.data.frame(exp)
    exp <- dplyr::select(exp, samps) #Selecting columns that match tumor subtype
    exp <- as.matrix(exp)
    
    names <- names(annot)[which(annot == subset_tumor, arr.ind = TRUE)[, "col"]]
    annot <- annot[[names[1]]]
    annot <- annot[annot == subset_tumor]
    heatmap <- ComplexHeatmap::HeatmapAnnotation(df = annot)
  }
  #If/else conditions to check if annotation object available
  if (is.null(heatmap)) {
    if(scale == FALSE){
      ComplexHeatmap::Heatmap(exp, name = "exposures",show_column_names = show_column_names, show_row_names = show_row_names,col = c("blue", "green", "yellow", "orange", "red"), ...)
    }
    else{
      ComplexHeatmap::Heatmap(exp, name = "exposures", show_column_names = show_column_names, show_row_names = show_row_names, ...)
    }
  }
    else if (!is.null(heatmap)) {
      if(scale == FALSE){
        ComplexHeatmap::Heatmap(exp, name = "exposures", top_annotation = heatmap, show_column_names = show_column_names, show_row_names = show_row_names,col = c("blue", "green", "yellow", "orange", "red"), ...)
      }
      else{
        ComplexHeatmap::Heatmap(exp, name = "exposures", top_annotation = heatmap, show_column_names = show_column_names, show_row_names = show_row_names, ...)
      }
    }
}
 
  

