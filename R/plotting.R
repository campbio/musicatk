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
  colnames(signatures) <- paste0(
    colnames(signatures),
    "-", seq_len(ncol(signatures))
  )
  names(sig_names) <- colnames(signatures)

  # Reformat signature matrix into long data.frame
  df <- signatures |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "motif") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(sig_names)),
      names_to = "signature",
      values_to = "exposure",
      names_repair = "minimal"
    )

  # Check for mutation color variable in annot table
  final_color_variable <- NULL
  if (is.null(color_variable) && !is.null(tab@color_variable)) {
    color_variable <- tab@color_variable
  }

  # Set up color variable if supplied as vector or the name of a column in
  # the table annotation
  if (length(color_variable) == 1 && color_variable %in% colnames(annot)) {
    final_color_variable <- annot[df$motif, color_variable]
  } else if (length(color_variable) == nrow(signatures)) {
    # One value per motif (row of signatures): match each row of df to it
    final_color_variable <-
      color_variable[match(df$motif, rownames(signatures))]
  } else {
    warning(
      "'color_variable' must be a column in the table annotation: ",
      paste(colnames(annot), collapse = ", "), ". Or it must be the ",
      "same length as the number of motifs in the signatures: ",
      nrow(signatures)
    )
  }

  # Save color variable to df if it was specified
  if (!is.null(final_color_variable)) {
    df <- cbind(df, mutation_color = final_color_variable)
  }

  # Make sure signature order is preserved using factor
  df$signature <- factor(df$signature, levels = names(sig_names))
  return(list(df = df, names = sig_names))
}

.add_small_legend <- function(my_plot, point_size = 2,
                              text_size = 10, space_legend = 0.5) {
  my_plot +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        override.aes =
          list(size = point_size)
      ),
      color = ggplot2::guide_legend(
        override.aes =
          list(size = point_size)
      )
    ) +
    ggplot2::theme(
      legend.title = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      legend.key.size = ggplot2::unit(space_legend, "lines"),
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.spacing.x = ggplot2::unit(0.25, "cm")
    )
}

.gg_default_theme <- function(p, text_size = 10) {
  p <- p + theme_bw() + theme(
    panel.grid = element_blank(),
    text = element_text(size = text_size)
  )

  if ("mono" %in% names(grDevices::pdfFonts())) {
    p <- p + theme(text = element_text(
      family = "mono",
      size = text_size
    ))
  }
  return(p)
}

# Colours for the values of the signature colour variable, named by value.
# 'values' are the values per motif; NULL color_mapping uses the table's
# default, or generated colours if that doesn't cover the values. An
# unnamed color_mapping is named in the order the values first appear.
.signature_color_mapping <- function(color_mapping, tab, values) {
  values <- unique(as.character(values))
  if (is.null(color_mapping)) {
    color_mapping <- get_color_mapping(tab)
    if (!all(values %in% names(color_mapping))) {
      color_mapping <- .discrete_colors(length(values))
      names(color_mapping) <- values
    }
    return(color_mapping)
  }
  if (is.null(names(color_mapping))) {
    if (length(color_mapping) != length(values)) {
      stop(
        "An unnamed 'color_mapping' needs one colour for each of the ",
        length(values), " values of 'color_variable'; it has ",
        length(color_mapping), "."
      )
    }
    names(color_mapping) <- values
  }
  missing <- setdiff(values, names(color_mapping))
  if (length(missing) > 0) {
    stop(
      "'color_mapping' has no colour for these values of ",
      "'color_variable': ", paste(missing, collapse = ", ")
    )
  }
  return(color_mapping)
}

# Black text on light fills, white text on dark fills
.label_text_colors <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  brightness <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  unname(ifelse(brightness < 0.25, "white", "black"))
}

# Divide each sample's (column's) exposures by its total. A sample with a
# total of 0 has no proportions, so stop and name it.
.proportional_exposures <- function(exposures) {
  total <- colSums(exposures)
  if (any(total == 0)) {
    stop(
      "Some samples have a total exposure of 0, so their proportions are ",
      "undefined. Remove these samples or use proportional = FALSE: ",
      paste(colnames(exposures)[total == 0], collapse = ", ")
    )
  }
  sweep(exposures, 2, total, FUN = "/")
}

.discrete_colors <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  colors <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  return(colors)
}
