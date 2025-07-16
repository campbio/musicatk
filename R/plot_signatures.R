#' @import ggplot2
#' @importFrom rlang .data
#' @import dplyr
#' @import tidyr
#' @import tidyverse
#' @importFrom scales pretty_breaks
NULL

#' Plots the mutational signatures
#'
#' After mutational signature discovery has been performed, this function
#' can be used to display the distribution of each mutational signature. The
#' \code{color_variable} and \code{color_mapping} parameters can be used
#' to change the default color scheme of the bars.
#'
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' discovery or prediction.
#' @param model_id The name of the model to plot.
#' @param modality The modality of the signatures to plot. Must be
#' "SBS96", "DBS78", or "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing the signatures
#' to plot. Default \code{"result"}.
#' @param color_variable Name of the column in the variant annotation data.frame
#' to use for coloring the mutation type bars. The variant annotation data.frame
#' can be found within the count table of the \code{\linkS4class{musica}}
#' object. If \code{NULL}, then the default column specified in the count
#' table will be used. Default \code{NULL}.
#' @param color_mapping A character vector used to map items in the
#' \code{color_variable} to a color. The items in \code{color_mapping}
#' correspond to the colors. The names of the items in \code{color_mapping}
#' should correspond to the unique items in \code{color_variable}. If
#' \code{NULL}, then the default \code{color_mapping} specified in the count
#' table will be used. Default \code{NULL}.
#' @param text_size Size of axis text. Default \code{10}.
#' @param show_x_labels If \code{TRUE}, the labels for the mutation types
#' on the x-axis will be shown. Default \code{TRUE}.
#' @param show_y_labels If \code{TRUE}, the y-axis ticks and labels will be
#' shown. Default \code{TRUE}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{FALE}.
#' @param y_max Vector of maximum y-axis limits for each signature. One value
#' may also be provided to specify a constant y-axis limit for all signatures.
#' Vector length must be 1 or equivalent to the number of signatures. Default
#' \code{NULL}.
#' @param annotation Vector of annotations to be displayed in the top right
#' corner of each signature. Vector length must be equivalent to the number of
#' signatures. Default \code{NULL}.
#' @param percent If \code{TRUE}, the y-axis will be represented in percent
#' format instead of mutation counts. Default \code{TRUE}.
#' @param plotly If \code{TRUE}, the the plot will be made interactive
#' using \code{\link[plotly]{plotly}}. Default \code{FALSE}.

#' @return Generates a ggplot or plotly object
#' @examples
#' data(res)
#' plot_signatures(res, model_id = "res")
#' @export
plot_signatures <- function(musica,
                            model_id,
                            modality = "SBS96",
                            result_name = "result",
                            color_variable = NULL,
                            color_mapping = NULL,
                            text_size = 10,
                            show_x_labels = TRUE,
                            show_y_labels = TRUE,
                            same_scale = FALSE,
                            y_max = NULL,
                            annotation = NULL,
                            percent = TRUE,
                            plotly = FALSE) {
  # dummy variables
  loc_num <- NULL
  mutation_color <- NULL
  label <- NULL
  x <- NULL
  xend <- NULL
  y <- NULL
  yend <- NULL
  ymax <- NULL

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

  # Get result object from musica object
  result <- get_model(musica,
    result = result_name,
    modality = modality,
    model = model_id
  )

  g <- .plot_result_model_signatures(
    result, musica, color_variable,
    color_mapping, text_size, show_x_labels,
    show_y_labels, same_scale, y_max, annotation,
    percent, plotly
  )

  return(g)
}

.plot_result_model_signatures <- function(result, musica,
                                          color_variable = NULL,
                                          color_mapping = NULL,
                                          text_size = 10,
                                          show_x_labels = TRUE,
                                          show_y_labels = TRUE,
                                          same_scale = FALSE,
                                          y_max = NULL,
                                          annotation = NULL,
                                          percent = TRUE,
                                          plotly = FALSE) {
  # dummy variables
  mutation_color <- NULL
  loc_num <- NULL
  label <- NULL
  ymax <- NULL
  x <- NULL
  xend <- NULL
  y <- NULL
  yend <- NULL
  exposure <- NULL
  motif <- NULL
  exposure_null <- NULL

  signatures <- signatures(result)
  sig_names <- colnames(signatures)
  table_name <- modality(result)
  tab <- tables(musica)[[table_name]]
  annot <- get_annot_tab(tab)
  num_sigs <- length(sig_names)

  if (is.null(color_mapping)) {
    color_mapping <- get_color_mapping(tab)
  }
  plot_dat <- .pivot_signatures(signatures, tab,
    color_variable = color_variable
  )

  width <- 0.45
  motif_label_locations <-
    plot_dat$df[plot_dat$df$signature == plot_dat$df[1, 2], ] %>%
    ungroup() %>%
    mutate(loc_num =
             c(seq_len(dim(plot_dat$df[plot_dat$df$signature ==
                                         plot_dat$df[1, 2], ])[1]))) %>%
    group_by(mutation_color) %>%
    summarise(
      x = min(loc_num) - width, xend = max(loc_num) + width,
      y = 0, yend = 0.01
    )

  # Whether to re-scale y axis
  scales <- ifelse(isTRUE(same_scale), "fixed", "free_y")

  # If annotation supplied
  if (!is.null(annotation)) {
    annotation_text <- data.frame(
      label = annotation,
      signature = names(plot_dat$names),
      mutation_color = names(color_mapping)[length(names(color_mapping))]
    )
  }

  # Rename signature labels
  sig_name_labels <- data.frame(
    label = sig_names,
    signature = names(plot_dat$names),
    mutation_color = names(color_mapping)[1]
  )

  # Add potential forced y-axis max
  plot_dat$df$ymax <- rep(y_max, length(unique(plot_dat$df$motif)))

  # Convert exposure probabilities to percentages
  if (percent == TRUE) {
    plot_dat$df$exposure <- plot_dat$df$exposure * 100
    max_num_digits <- floor(log10(max(plot_dat$df$exposure) * 1.2)) + 1
    y_axis_label <- "Percent of Mutations"
    y_axis_spacing <- rep(strrep(" ", max_num_digits), 2)
  } else {
    y_axis_label <- "Mutation Counts"
    max_num_digits <- floor(log10(max(plot_dat$df$exposure) * 1.2)) + 1
    if (max(plot_dat$df$exposure) == 1) {
      max_num_digits <- 3
    }
    if (max(plot_dat$df$exposure) < 1) {
      max_num_digits <- 2
      y_axis_label <- "Mutation Probability"
    }

    y_axis_spacing <- rep(strrep(" ", max_num_digits), 2)
  }

  if (is.null(plot_dat$df$context)) {
    plot_dat$df$context <- annot$context
  }

  # Plot signatures
  plot_dat$df %>%
    ggplot(aes(y = exposure, x = motif, fill = mutation_color)) +
    geom_bar(stat = "identity") +
    facet_grid(factor(signature) ~ ., scales = scales) +
    ggplot2::xlab("Motifs") +
    ggplot2::ylab(y_axis_label) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::scale_fill_manual(values = color_mapping) +
    #ggplot2::scale_x_discrete(labels = annot$context) +
    ggplot2::scale_x_discrete(limits = plot_dat$df$motif, labels = plot_dat$df$context) +
    ggplot2::scale_y_continuous(
      expand = expansion(mult = c(0, 0.2)),
      limits = c(0, NA), n.breaks = 5,
      breaks = scales::pretty_breaks()
    ) +
    ggplot2::geom_text(
      data = sig_name_labels,
      mapping = aes(x = -Inf, y = Inf, label = label),
      hjust = -0.075, vjust = 1.5,
      fontface = "bold"
    ) -> p

  # Add annotations, if necessary
  if (exists("annotation_text") == TRUE) {
    p <- p + ggplot2::geom_text(
      data = annotation_text,
      mapping = aes(x = Inf, y = Inf, label = label),
      hjust = 1.075, vjust = 1.5,
      fontface = "bold"
    )
  }

  # Change y-axis scale, if necessary
  if (!is.null(y_max)) {
    p <- p + geom_blank(aes(y = ymax))
  }
  
  # Adjust theme
  p <- .gg_default_theme(p, text_size = text_size) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 1, 0, 1)) +
    theme(strip.background = element_blank(), strip.text.y = element_blank())
  

  # If SBS, need to change color of labels so one is white
  if (table_name == "SBS96") {
    label_colors <- c("black", "white", rep("black", 4))
  } else {
    label_colors <- c(rep("black", length(color_mapping)))
  }

  # Plot motif labels
  plot_dat$df$exposure_null <- rep(0, dim(plot_dat$df)[1])
  plot_dat$df %>%
    ggplot(aes(y = exposure_null, x = motif)) +
    geom_bar(stat = "identity") +
    ggplot2::scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
      limits = c(0, NA), breaks = c(0, 0.01),
      labels = y_axis_spacing, n.breaks = 4
    ) +
    ggplot2::ylab("") +
    ggplot2::geom_rect(
      #data = motif_label_locations,
      data = motif_label_locations %>% arrange(x),
      aes(xmin = x, xmax = xend, ymin = max(y), ymax = max(yend)),
      #fill = color_mapping, color = "black",
      #linewidth = 0.25, inherit.aes = FALSE
    #) +
    #ggplot2::geom_text(
      #data = motif_label_locations,
      #aes(
        #x = x + (xend - x) / 2, y = y + (yend - y) / 2,
        #label = stringr::str_to_title(mutation_color)
      #),
      #fontface = "bold", size = 4,
      #color = label_colors
    #) -> p2
    
    fill = factor(color_mapping, levels = color_mapping), color = "black", 
    linewidth = 0.25, inherit.aes = FALSE) -> p2 

  # adjust motif label direction if indel signature
  if (table_name %in% c("IND83", "INDEL83", "INDEL", "IND", "indel", 
                        "Indel")){
    p2 <- p2 + ggplot2::geom_text(data=motif_label_locations, 
                                  aes(x=x+(xend-x)/2, y=y+(yend-y)/2, 
                                      label = stringr::str_to_title(mutation_color)), 
                                  fontface = "bold", size = 4, 
                                  color = label_colors, angle = 90)
  }else{
    p2 <- p2 + ggplot2::geom_text(data=motif_label_locations, 
                                  aes(x=x+(xend-x)/2, y=y+(yend-y)/2, 
                                      label = stringr::str_to_title(mutation_color)), 
                                  fontface = "bold", size = 4, 
                                  color = label_colors)
  }
  
  # Adjust theme
  p2 <- .gg_default_theme(p2, text_size = text_size) +
    theme(plot.margin = margin(0, 1, 0, 1)) + # see function below
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )

  if (!isTRUE(show_x_labels)) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    )
  }

  if (!isTRUE(show_y_labels)) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
    p2 <- p2 + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
  }
  
  # adjust height of motif lables if indel signature
  if (table_name %in% c("IND83", "INDEL83", "INDEL", "IND", "indel", 
                        "Indel")){
    height <- 5
  }else{
    height <- 1
  }

  # combine signatures and motif labels
  figure <- ggpubr::ggarrange(p2, p, ncol = 1, nrow = 2, heights = c(height,15))

  if (isTRUE(plotly)) {
    figure <- plotly::ggplotly(p)
  }

  return(figure)
}
