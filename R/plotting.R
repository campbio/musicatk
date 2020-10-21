#' @importFrom ggplot2 ggplot aes geom_bar theme theme_set geom_point aes_string
#' @importFrom ggplot2 facet_grid theme_bw xlab ylab element_blank element_text
#' @importFrom rlang .data
#' @importFrom withr with_seed
NULL

#' Return sample from musica object
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param sample_names Names or indices of samples to plot.
#' @param table_name Name of table used for plotting counts. If \code{NULL},
#' then the first table in the \code{\linkS4class{musica}} object will be used.
#' Default \code{NULL}.
#' @return Generates sample plot {no return}
#' @examples
#' data(musica_sbs96)
#' plot_sample_counts(musica_sbs96, sample_names = 
#' sample_names(musica_sbs96)[1])
#' @export
plot_sample_counts <- function(musica, sample_names, table_name = NULL) {

  if(is.null(table_name)) {
   table_name <- names(tables(musica))[1]
  }  
  
  # Extract counts for specific samples
  tab <- .extract_count_table(musica, table_name)
  ix <- match(sample_names, colnames(tab))
  if(all(is.na(ix))) {
    stop("The values in 'sample_names' did not match any sample IDs in table '",
         table_name, "'.")
  }
  else if(anyNA(ix)) {
    warning("The following samples in 'sample_names' were not found  in ",
            "table '", table_name, 
            "' and will ", "be exlcuded from the plot: ",
            paste(sample_names[is.na(ix)], collapse = ", "))
    ix <- ix[!is.na(ix)]
  }
  sample_counts <- tab[, ix, drop = FALSE]
  
  result <- methods::new("musica_result",
                          signatures = sample_counts,exposures = matrix(),
                          type = "sample", musica = musica,
                          tables = table_name)
  g <- plot_signatures(result) + ggplot2::ylab("Mutation Counts")
  return(g)
}


#' Plotting Signature Motif counts/spectra
#'
#' @param result S4 Result Object
#' @param legend Whether to include the legend for mutation types in the plot.
#' @param plotly add plotly layer for plot interaction
#' @param color_variable Annotation column to use for coloring plotted motifs,
#' provided by counts table from input result's musica object
#' @param color_mapping Mapping from color_variable to color names, provided by
#' counts table from input result's musica object
#' @param text_size Size of axis text
#' @param facet_size Size of facet text
#' @param show_x_labels Toggle plotting of x-axis labels
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{TRUE}.
#' @return Generates plot {no return}
#' @examples
#' data(res)
#' plot_signatures(res)
#' @export
plot_signatures <- function(result, legend = TRUE, plotly = FALSE,
                            color_variable = NULL, color_mapping = NULL,
                            text_size = 10, facet_size = 10,
                            show_x_labels = TRUE,
                            same_scale = TRUE) {
  signatures <- signatures(result)
  sig_names <- colnames(signatures)
  table_name <- table_name(result)
  tab <- tables(result)[[table_name]]
  annot <- get_annot_tab(tab)

  if(is.null(color_mapping)) {
    color_mapping <- get_color_mapping(tab)
  }
  plot_dat <- .pivot_signatures(signatures, tab,
                                color_variable = color_variable)
  
  # Wether to rescale y axis
  scales <- ifelse(isTRUE(same_scale), "fixed", "free_y")
  
  plot_dat$df %>%
    ggplot(aes_string(y = "exposure", x = "motif", fill = "mutation_color")) +
    geom_bar(stat = "identity") +
    facet_grid(factor(signature, ordered = TRUE) ~ .,
               scales = scales,
               labeller = ggplot2::as_labeller(structure(plot_dat$names,
                                  names = names(plot_dat$names)))) +
    xlab("Motifs") + ylab("Probability") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = color_mapping) +
    ggplot2::scale_x_discrete(labels = annot$context) -> p

  # Adjust theme
  p <- .gg_default_theme(p, text_size = text_size, facet_size = facet_size) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  if (!isTRUE(show_x_labels)) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   panel.grid.major.x = element_blank())
  }
  if (!isTRUE(legend)) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- .addSmallLegend(p) + theme(legend.position = "bottom",
                                    legend.title = element_blank())
  }
  if (isTRUE(plotly)) {
    p <- plotly::ggplotly(p)
  }
  return(p)
}

#' Plotting Signature Motif counts/spectra
#'
#' @param result S4 Result Object
#' @param sample Name or index of sample within result object
#' @param plotly add plotly layer for plot interaction
#' @return Generates plot {no return}
#' @examples
#' data(res)
#' plot_sample_reconstruction_error(res, "TCGA-ER-A197-06A-32D-A197-08")
#' @export
plot_sample_reconstruction_error <- function(result, sample,
                                             plotly = FALSE) {
  signatures <- .extract_count_table(musica(result), 
                                     table_name(result))[, sample, drop = FALSE]
  sample_name <- colnames(signatures)
  reconstructed <- reconstruct_sample(result, sample)
  sigs <- cbind(signatures, reconstructed, signatures - reconstructed)
  colnames(sigs) <- c("Counts", "Reconstructed", "Difference")
  
  recontruct_result <- methods::new("musica_result",
                      signatures = sigs,
                      exposures = matrix(), type = "NMF",
                      musica = musica(result),
                      tables = table_name(result))
  plot_signatures(recontruct_result, same_scale = FALSE) +
    ggplot2::ggtitle("Reconstruction error", subtitle = sample_name) + ylab("")
}


# Utility functions -------------------------------
.pivot_signatures <- function(signatures, tab, sig_names = NULL,
                              color_variable = NULL) {
  if(is.null(sig_names)) {
    sig_names <- colnames(signatures)  
  }
  annot <- tab@annotation
  
  # Ensure signature colnames are unique
  # They can not be unique in the sig_compare function if one signature
  # is matched up against several others in the second result object
  colnames(signatures) <- paste0(colnames(signatures),
                                 "-", seq(ncol(signatures)))
  names(sig_names) <- colnames(signatures)
    
  # Rormat signature matrix into long data.frame
  signatures %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "motif") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(names(sig_names)),
                        names_to = "signature",
                        values_to = "exposure",
                        names_repair = "minimal") -> df
  
  # Check for mutation color variable in annot table
  final_color_variable <- NULL
  if(is.null(color_variable) && !is.null(tab@color_variable)) {
    color_variable <- tab@color_variable
  }
  
  # Set up color variable if supplied as vector or the name of a column in
  # the table annotation
  if(length(color_variable) == 1 && color_variable %in% colnames(annot)) {
    final_color_variable <- annot[df$motif,tab@color_variable]
  } else if (length(color_variable) == nrow(signatures)) {
    final_color_variable <- color_variable
  } else {
    warning("'color_variable' must be a column in the table annotation: ",
            paste(colnames(annot), collapse = ", "), ". Or it must be the ",
            "same length as the number of motifs in the signatures: ",
            nrow(signatures))
  }
  
  # Save color variable to df if it was specified
  if(!is.null(final_color_variable)) {
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
            legend.spacing.x = ggplot2::unit(0.25, 'cm'))
}

.gg_default_theme <- function(p, text_size = 10, facet_size = 10) {
  p <- p + theme_bw() + theme(
    strip.text.y = element_text(size = facet_size),
    panel.grid = element_blank(),
    text = element_text(family = "Courier",
                        size = text_size))
  return(p)
}

.discrete_colors <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  colors <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  return(colors)
}

