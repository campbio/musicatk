#' @import ggplot2 
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' @import tidyverse 
NULL

# Plotting utility functions -------------------------------
.pivot_signatures <- function(signatures, tab, sig_names = NULL,
                              color_variable = NULL) {
  if (is.null(sig_names)) {
    sig_names <- colnames(signatures)  
  }
  annot <- tab@annotation
  rownames(annot) <- annot$motif
  
  # Ensure signature colnames are unique
  # They can not be unique in the sig_compare function if one signature
  # is matched up against several others in the second result object
  colnames(signatures) <- paste0(colnames(signatures),
                                 "-", seq(ncol(signatures)))
  names(sig_names) <- colnames(signatures)
    
  # Reformat signature matrix into long data.frame
  signatures %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "motif") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(sig_names)),
                        names_to = "signature",
                        values_to = "exposure",
                        names_repair = "minimal") -> df
  
  # Check for mutation color variable in annot table
  final_color_variable <- NULL
  if (is.null(color_variable) && !is.null(tab@color_variable)) {
    color_variable <- tab@color_variable
  }
  
  # Set up color variable if supplied as vector or the name of a column in
  # the table annotation
  if (length(color_variable) == 1 && color_variable %in% colnames(annot)) {
    final_color_variable <- annot[df$motif, tab@color_variable]
  } else if (length(color_variable) == nrow(signatures)) {
    final_color_variable <- color_variable
  } else {
    warning("'color_variable' must be a column in the table annotation: ",
            paste(colnames(annot), collapse = ", "), ". Or it must be the ",
            "same length as the number of motifs in the signatures: ",
            nrow(signatures))
  }
  
  # Save color variable to df if it was specified
  if (!is.null(final_color_variable)) {
    df <- cbind(df, mutation_color = final_color_variable)
  }

  # Make sure signature order is preserved using factor
  df$signature <- factor(df$signature, levels = names(sig_names))
  return(list(df = df, names = sig_names))
}

.addSmallLegend <- function(myPlot, pointSize = 2,
                            textSize = 10, spaceLegend = 0.5) {
  myPlot +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes =
                                                    list(size = pointSize)),
                    color = ggplot2::guide_legend(override.aes =
                                                    list(size = pointSize))) +
    ggplot2::theme(legend.title = element_text(size = textSize),
            legend.text  = element_text(size = textSize),
            legend.key.size = ggplot2::unit(spaceLegend, "lines"),
            legend.box.background = ggplot2::element_rect(colour = "black"),
            legend.spacing.x = ggplot2::unit(0.25, "cm"))
}

.gg_default_theme <- function(p, text_size = 10) {
  p <- p + theme_bw() + theme(
    panel.grid = element_blank(),
    text = element_text(size = text_size))
  
  if("mono" %in% names(grDevices::pdfFonts())) {
    p <- p + theme(text = element_text(family = "mono",
                                       size = text_size))
  }
  return(p)
}

.discrete_colors <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  colors <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  return(colors)
}
