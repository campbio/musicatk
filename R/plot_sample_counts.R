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
#' @param table_name Name of table used for plotting counts. If \code{NULL},
#' then the first table in the \code{\linkS4class{musica}} object will be used.
#' Default \code{NULL}.
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
plot_sample_counts <- function(musica, sample_names, table_name = NULL, 
                               text_size = 10,
                               show_x_labels = TRUE, show_y_labels = TRUE,
                               same_scale = TRUE, annotation = NULL) {
    
    if (is.null(table_name)) {
        table_name <- names(tables(musica))[1]
    }  
    
    # Extract counts for specific samples
    tab <- .extract_count_table(musica, table_name)
    ix <- match(sample_names, colnames(tab))
    if (all(is.na(ix))) {
        stop("The values in 'sample_names' did not match any sample IDs in table '",
             table_name, "'.")
    }
    else if (anyNA(ix)) {
        warning("The following samples in 'sample_names' were not found  in ",
                "table '", table_name, 
                "' and will ", "be exlcuded from the plot: ",
                paste(sample_names[is.na(ix)], collapse = ", "))
        ix <- ix[!is.na(ix)]
    }
    sample_counts <- tab[, ix, drop = FALSE]
    
    result <- methods::new("musica_result",
                           signatures = sample_counts, exposures = matrix(),
                           algorithm = "sample", musica = musica,
                           table_name = table_name)
    g <- plot_signatures(result, percent = FALSE, text_size = text_size,
                         show_x_labels = show_x_labels, show_y_labels = show_y_labels,
                         same_scale = same_scale, annotation = annotation)
    return(g)
}

