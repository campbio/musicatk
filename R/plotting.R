#' @importFrom ggplot2 ggplot aes geom_bar theme theme_set geom_point aes_string
#' @importFrom ggplot2 facet_grid theme_bw xlab ylab element_blank element_text
#' @importFrom rlang .data
NULL

#' Return sample from musica object
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of table used for plotting counts
#' @param sample_name Sample name to plot counts
#' @return Generates sample plot {no return}
#' @examples
#' musica <- readRDS(system.file("testdata", "musica_sbs96.rds", package = "musicatk"))
#' plot_sample_counts(musica, "SBS96", get_sample_names(musica)[1])
#' @export
plot_sample_counts <- function(musica, table_name, sample_name) {
  if (length(sample_name) != 1) {
    stop("`please specify exactly one sample")
  }
  sample_number <- which(get_sample_names(musica = musica) == sample_name)
  sample <- .extract_count_table(musica, table_name)[, sample_number]
  plot_full(sample)
}

#' Complete Plotting of input data
#'
#' @param sample Single sample DataFrame
#' @return Generates sample plot {no return}
plot_full <- function(sample) {
  mut_summary <- table_96(sample)
  major <- table(mut_summary[, "Type"])
  df <- data.frame(major)
  mut <- data.frame(table(mut_summary$mutation))

  plot_major <- ggplot(df, aes_string(x = "Var1", y = "Freq",
                                               fill = "Var1")) + geom_bar(stat =
                                                               "identity") +
    theme(axis.text.x = element_text(angle = 70, vjust = 0.5),
          text = element_text(family = "Courier"), legend.position = "none") +
    xlab("Major Motif") + ylab("Counts")  + ggplot2::scale_fill_manual(
      values = c("#5ABCEBFF", "#050708FF", "#D33C32FF", "#CBCACBFF",
                 "#ABCD72FF", "#E7C9C6FF"))

  df2 <- data.frame(major = substr(mut$Var1, 1, 3), minor =
                      substr(mut$Var1, 5, 7), num = mut$Freq)
  df2$major <- as.character(df2$major)
  df2$minor <- as.character(df2$minor)
  plot_iris2 <- df2 %>% ggplot(aes_string(y = "num", x = "minor",
                                          fill = "major")) +
    ggplot2::facet_grid(~ major, scales = "free_x") +
    geom_bar(stat = "identity") + cowplot::background_grid(major = "y",
                                                           minor = "none") +
    cowplot::panel_border() + theme(axis.text.x = element_text(angle = 90,
                                                               vjust = 0.5,
                                                               hjust = 1),
                                    text = element_text(family = "Courier")) +
    xlab("Minor Motif") + ylab("Counts") + ggplot2::labs(fill = "SNP") +
    ggplot2::scale_fill_manual(values = c(
      "#5ABCEBFF", "#050708FF", "#D33C32FF", "#CBCACBFF", "#ABCD72FF",
      "#E7C9C6FF"))
  legend <- cowplot::get_legend(plot_iris2)
  plot_iris2 <- plot_iris2 + theme(legend.position = "none")

  theme_set(cowplot::theme_cowplot(font_size = 8)) # reduce default font size

  p <- cowplot::plot_grid(plot_iris2, legend, plot_major, labels = "AUTO",
                          ncol = 2, rel_widths = c(10, 1))
  return(p)
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
#' result <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' plot_signatures(result)
#' @export
plot_signatures <- function(result, legend = TRUE, plotly = FALSE,
                            color_variable = NULL, color_mapping = NULL,
                            text_size = 10, facet_size = 10,
                            show_x_labels = TRUE,
                            same_scale = TRUE) {
  signatures <- result@signatures
  sig_names <- colnames(signatures)
  table_name <- result@tables
  tab <- result@musica@count_tables[[table_name]]
  annot <- tab@annotation

  if(is.null(color_mapping)) {
    color_mapping <- tab@color_mapping
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
    ggplot2::scale_x_discrete(labels=annot$context) -> p

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
    p <- .addSmallLegend(p) + theme(legend.position="bottom",
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
#' result <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' plot_sample_reconstruction_error(result, "TCGA-ER-A197-06A-32D-A197-08")
#' @export
plot_sample_reconstruction_error <- function(result, sample,
                                             plotly = FALSE) {
  signatures <- .extract_count_table(result@musica, result@tables)[, sample,
                                                              drop = FALSE]
  sample_name <- colnames(signatures)
  reconstructed <- reconstruct_sample(result, sample)
  sigs <- cbind(signatures, reconstructed, signatures - reconstructed)
  colnames(sigs) <- c("Counts", "Reconstructed", "Difference")
  
  recontruct_result <- methods::new("musica_result",
                      signatures = sigs,
                      exposures = matrix(), type = "NMF",
                      musica = result@musica,
                      tables = result@tables)
  plot_signatures(recontruct_result, same_scale = FALSE) +
    ggplot2::ggtitle("Reconstruction error", subtitle = sample_name) + ylab("")
}

#' Plots signature weights for each sample
#'
#' @param result S4 Result Object
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @param label_samples Whether to display sample names or not
#' @param samples_plotted A subset of samples to plot
#' @param sort_samples Defaults to numerical, may be total 'counts', or a
#' specific signatures name (e.g. 'Signature1)
#' @param num_samples Number of sorted samples to plot
#' @param thresh_zero Max level to zero out for better plotting when sorting
#' by multiple signatures
#' @param legend Don't plot legend
#' @param plotly add plotly layer for plot interaction
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res.rds", package = "musicatk"))
#' plot_exposures(result)
#' @export
plot_exposures <- function(result, proportional = TRUE, label_samples = FALSE,
                           samples_plotted = colnames(result@exposures),
                           sort_samples = "numerical",
                           num_samples = length(colnames(result@exposures)),
                           thresh_zero = FALSE, legend = TRUE,
                           plotly = FALSE) {
  samples <- result@exposures
  if (thresh_zero) {
    samples[samples < thresh_zero] <- 0
  }
  if (!all(samples_plotted %in% colnames(samples))) {
    stop(strwrap(prefix = " ", initial = "", "Some samples specified are
    not in this dataset, available samples are as follows: "), "\n",
         paste0(colnames(samples), collapse = "\n"))
  }
  samples <- samples[, samples_plotted, drop = FALSE]
  y_label <- "counts"
  counts <- colSums(samples)
  if (proportional) {
    y_label <- "%"
    samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  }
  samples %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "make") %>%
    tidyr::gather("var", "val", -"make") -> plot_dat
  if (length(sort_samples == 1) && sort_samples == "numerical") {
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(gtools::mixedsort(plot_dat$var)))
  }else if (length(sort_samples == 1) && sort_samples == "counts") {
    #counts <- colSums(samples)
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(plot_dat$var)
                           [order(counts, decreasing = TRUE)])
  }else {
    if (!all(sort_samples %in% rownames(samples))) {
      stop("Signature is not present in this result, please choose from: \n",
           paste0(rownames(samples), collapse = "\n"))
    }
    if (length(sort_samples) == 1) {
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[order(samples[
                                 sort_samples, ], decreasing = TRUE)])
    }
    else{
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[do.call(order,
                                                            c(decreasing = TRUE,
                                               as.data.frame(
                                                 t(samples[sort_samples, ]))))])
    }
  }
  samples_to_use <- levels(plot_dat$var)[seq_len(num_samples)]
  plot_dat <- plot_dat[which(plot_dat$var %in% samples_to_use), ]


  #Testing
  plot_dat_table <- plot_dat %>% dplyr::arrange(.data$make)
  if (all(!sort_samples %in% c("numerical", "counts"))) {
    plot_levels <- c(setdiff(unique(plot_dat_table$make), sort_samples),
                     rev(sort_samples))
    plot_dat_table$make <- factor(plot_dat_table$make, levels = plot_levels)
  }

  plot_dat_table %>%
    ggplot() + geom_bar(aes_string(y = "val", x = "var", fill = "make"),
                        stat = "identity") + theme_bw() + theme(legend.title =
                                                  element_blank(), axis.text.x =
                         element_text(angle = 90, hjust = 1, vjust = 0.5),
                       panel.grid.major.x = element_blank(),
                       text = element_text(family = "Courier")) +
    xlab("Samples") + ylab(paste("Signatures (", y_label, ")", sep = "")) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) -> p
  if (!label_samples) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x =
                     element_blank())
  }
  if (isTRUE(legend)) {
    p <- p + theme(legend.position = "none")
  }
  if (plotly) {
    p <- plotly::ggplotly(p, tooltip = c("x", "y"))
  }
  return(p)
}

#' Plots signature weights for each sample gridded by single annotation
#'
#' @param result S4 Result Object
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @param label_samples Whether to display sample names or not
#' @param sort_samples Defaults to numerical, may be total 'counts', or a
#' specific signatures name (e.g. 'Signature1)
#' @param annotation Annotation to use for facet_wrap
#' @param by_group Plot subplots by annotation or by signature
#' @param num_samples Number of sorted samples to plot
#' @param no_legend Remove legend from the plot
#' @param thresh_zero Max level to zero out for better plotting when sorting
#' by multiple signatures
#' @param plotly add plotly layer for plot interaction
#' @return Generates plot {no return}
#' @examples
#' result <- readRDS(system.file("testdata", "res_annot.rds",
#' package = "musicatk"))
#' plot_exposures_by_annotation(result, "Tumor_Subtypes")
#' @export
plot_exposures_by_annotation <- function(result, annotation,
                                         proportional = TRUE,
                                         label_samples = FALSE,
                                         sort_samples = "numerical",
                                         by_group = TRUE, num_samples = FALSE,
                                         no_legend = FALSE, thresh_zero = FALSE,
                                         plotly = FALSE) {
  samples <- result@exposures
  if (thresh_zero) {
    samples[samples < thresh_zero] <- 0
  }
  y_label <- "counts"
  if (proportional) {
    y_label <- "%"
    samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  }
  samples %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "make") %>%
    tidyr::gather("var", "val", -"make") -> plot_dat
  annot_dt <- result@musica@sample_annotations
  plot_dat$annotation <- annot_dt[[annotation]][match(plot_dat$var,
                                                      annot_dt[["Samples"]])]
  #Normal alphabetical sorting of samples
  if (length(sort_samples == 1) && sort_samples == "numerical") {
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(gtools::mixedsort(plot_dat$var)))
    #Sorting of samples by counts
  }else if (length(sort_samples == 1) && sort_samples == "counts") {
    counts <- colSums(samples)
    plot_dat$var <- factor(plot_dat$var, levels =
                             unique(plot_dat$var)
                           [order(counts, decreasing = TRUE)])
    #Sorting of samples by counts of specific signature(s)
  }else {
    if (!all(sort_samples %in% rownames(samples))) {
      stop("Signature is not present in this result, please choose from: \n",
           paste0(rownames(samples), collapse = "\n"))
    }
    if (length(sort_samples) == 1) {
      plot_dat$var <- factor(plot_dat$var, levels =
                               unique(plot_dat$var)[order(samples[
                                 sort_samples, ], decreasing = TRUE)])
    }
    else {
      plot_dat$var <- factor(plot_dat$var, levels = unique(plot_dat$var)[
        do.call(order, c(decreasing = TRUE, as.data.frame(t(samples[
          sort_samples, ]))))])
    }
  }
  if (num_samples) {
    if (all(!sort_samples %in% c("numerical", "counts"))) {
      sub_dat <- plot_dat[plot_dat$make == sort_samples, ]
    } else {
      sub_dat <- plot_dat
    }
    data.table::setorderv(sub_dat, cols = "val", order = -1)
    annots <- unique(plot_dat$annotation)
    samples_to_use <- NULL
    for (annot in annots) {
      samples_to_use <- c(samples_to_use, as.character(utils::head(sub_dat[
        sub_dat$annotation == annot, "var"], num_samples)))
    }
    plot_dat <- plot_dat[which(plot_dat$var %in% samples_to_use), ]
  }

  plot_dat_table <- plot_dat %>% dplyr::arrange(.data$make)

  #If used, sorts indicated signature to the bottom so it's easier to see
  if (all(!sort_samples %in% c("numerical", "counts"))) {
    plot_levels <- c(setdiff(unique(plot_dat_table$make), sort_samples),
                     rev(sort_samples))
    plot_dat_table$make <- factor(plot_dat_table$make, levels = plot_levels)
  }
  if (by_group) {
    plot_dat_table %>%
      ggplot(aes_string(y = "val", x = "var", fill = "make",
                        text = "annotation")) + geom_bar(
                          stat = "identity") + theme_bw() + theme(
                            legend.title = element_blank(), axis.text.x =
              element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major.x = element_blank()) + xlab("Samples") +
      ylab(paste("Signatures (", y_label, ")", sep = "")) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) -> p
    p <- p + ggplot2::facet_wrap(~ annotation, drop = TRUE, scales = "free")
  } else {
    plot_dat_table %>%
      ggplot(aes_string(x = "annotation", y = "val", fill = "annotation",
                        text = "annotation")) +
      ggplot2::geom_boxplot() + theme_bw() + theme(legend.title =
                                                     element_blank(),
                                                   axis.text.x =
              element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid.major.x = element_blank()) + xlab("Annotation") +
      ylab(paste("Signatures (", y_label, ")", sep = "")) -> p
    #See if this makes plotly work
    p <- p + geom_point(alpha = 0)

    p <- p + ggplot2::facet_wrap(~ make, drop = FALSE, scales = "free")
  }
  p <- p + theme(text = element_text(family = "Courier", size = 10),
                 strip.text.y = element_text(size = 5))
  if (!label_samples) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x =
                     element_blank())
  }
  if (no_legend) {
    p <- p + theme(legend.position = "none")
  }
  if (plotly) {
    p <- plotly::ggplotly(p)#, hoverinfo = "legendgroup+x", hovermode = "x")
  }
  return(p)
}

#' Create a UMAP data.frame from a result object
#'
#' @param result S4 Result Object
#' @param annotation Annotation to use for umap coloring
#' @param n_neighbors The size of local neighborhood used for views of
#' manifold approximation. Larger values result in more global the manifold,
#' while smaller values result in more local data being preserved. Default 30.
#' See `?uwot::umap` for more information.
#' @param min_dist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where
#' nearby points on the manifold are drawn closer together, while larger
#' values will result on a more even dispersal of points. Default 0.2.
#' See `?uwot::umap` for more information.
#' @param spread The effective scale of embedded points. In combination with
#' ‘min_dist’, this determines how clustered/clumped the embedded points are.
#' Default 1.
#' See `?uwot::umap` for more information.
#' @param proportional Whether weights are normalized to sum to 1 or not
#' @param seed Use a seed for reproducible results
#' @return UMAP data.frame
#' @examples
#' result <- readRDS(system.file("testdata", "res_annot.rds",
#' package = "musicatk"))
#' create_umap(result = result, annotation = "Tumor_Subtypes",
#' n_neighbors = 5)
#' @export
create_umap <- function(result, annotation, n_neighbors = 30, min_dist = 0.75,
                        spread = 1, proportional = TRUE, seed = FALSE) {
  samples <- t(result@exposures)
  range_limit <- function(x) {
    (x / max(x))
    }
  if (proportional) {
    samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  }
  if (seed) {
    umap_out <- withr::with_seed(seed, uwot::umap(samples, n_neighbors =
                                                     n_neighbors, min_dist =
                                                     min_dist, spread = spread,
                                                   n_threads = 1, pca = NULL,
                                                   metric = "cosine"))
  } else {
    umap_out <- uwot::umap(samples, n_neighbors = n_neighbors,
                            min_dist = min_dist, spread = spread,
                            n_threads = 1, pca = NULL, metric = "cosine")
  }
  x <- umap_out[, 1]
  y <- umap_out[, 2]
  annot <- result@musica@sample_annotations
  samp_ind <- match(rownames(samples), annot$Samples)
  df <- data.frame(x = x, y = y, type = annot[[annotation]][samp_ind],
                   samp = rownames(samples))

  sig_df <- NULL
  n_sigs <- ncol(samples)
  n_samples <- nrow(samples)
  sig_names <- colnames(samples)
  for (i in 1:n_sigs) {
    sig_name <- sig_names[i]
    sig_df <- rbind(sig_df, data.frame(x = x, y = y, level =
                                                     range_limit(
                                                       samples[, sig_name]),
                                                   Signatures = rep(sig_name,
                                                                    n_samples)))
  }
  umaps <- list(umap_df = df, umap_df_sigs = sig_df, umap_type =
                     ifelse(proportional, "Proportional", "Counts"))
  eval.parent(substitute(result@umap <- umaps))
}

#' Plot a UMAP data.frame
#'
#' @param result Result object containing UMAP data.frame
#' @param point_size Scatter plot point size
#' @param legend Remove legend
#' @param label_clusters Add annotation labels to clusters (may not work well
#' for split or small clusters)
#' @param label_size Size of cluster labels
#' @param legend_size Set legend size
#' @param text_box Place a box around cluster labels for improved readability
#' @param plotly Create plotly version of plot
#' @return Returns a ggplot2 plot of the created umap, if plotly = TRUE the
#' ggplotly object is returned
#' @examples
#' result <- readRDS(system.file("testdata", "res_annot.rds",
#' package = "musicatk"))
#' create_umap(result, "Tumor_Subtypes", n_neighbors = 5)
#' plot_umap(result)
#' @export
plot_umap <- function(result, point_size = 0.7, legend = TRUE,
                      label_clusters = TRUE, label_size = 3, legend_size = 3,
                      text_box = TRUE, plotly = FALSE) {
  umap_df <- result@umap$umap_df
  cluster <- as.character(umap_df$type)
  if (plotly) {
    p <- ggplot(umap_df, aes_string(x = "x", y = "y", col = "type",
                                    text = "samp")) +
      geom_point(size = point_size) + ggplot2::ggtitle("UMAP")
    p <- plotly::ggplotly(p, tooltip = c("text", "x", "y", "type"))
  } else if (label_clusters) {
    p <- ggplot(umap_df, aes_string(x = "x", y = "y", col = "type")) +
      geom_point(size = point_size) + ggplot2::ggtitle("UMAP")
    if (isTRUE(legend)) {
      p <- p + theme(legend.position = "none")
    }
    centroid_list <- lapply(unique(cluster), function(x) {
      df_sub <- umap_df[umap_df$type == x, ]
      median_1 <- stats::median(df_sub[, "x"])
      median_2 <- stats::median(df_sub[, "y"])
      cbind(median_1, median_2, x)
    })
    centroid <- do.call(rbind, centroid_list)
    centroid <- data.frame(
      x = as.numeric(centroid[, 1]),
      y = as.numeric(centroid[, 2]),
      type = centroid[, 3]
    )
    p <- p + ggplot2::geom_point(data = centroid, mapping =
                                   ggplot2::aes_string(x = "x", y = "y"),
                                 size = 0, alpha = 0) +
      theme(legend.position = "none")
    if (text_box) {
      p <- p + ggrepel::geom_label_repel(data = centroid, mapping =
                                           ggplot2::aes_string(label = "type"),
                                         size = label_size)

    } else {
      p <- p + ggrepel::geom_text_repel(data = centroid, mapping =
                                          ggplot2::aes_string(label = "type"),
                                        size = label_size)
    }
  }
  return(p)
}

#' Plot a UMAP data.frame
#'
#' @param result Result object containing UMAP data.frame
#' @return Returns ggplot2 plot of the created umap
#' @examples
#' result <- readRDS(system.file("testdata", "res_annot.rds",
#' package = "musicatk"))
#' create_umap(result, "Tumor_Subtypes", n_neighbors = 5)
#' plot_umap_sigs(result)
#' @export
plot_umap_sigs <- function(result) {
  umap_df_sigs <- result@umap$umap_df_sigs
  ggplot(umap_df_sigs, aes_string(x = "x", y = "y", colour = "level")) +
    ggplot2::facet_wrap(~ Signatures, drop = TRUE, scales = "free") +
    geom_point(aes_string(alpha = "level", size = "level")) +
    ggplot2::scale_colour_gradientn(colours = c("grey", "red", "blue"),
                                    breaks = c(0, 0.0001, 0.1)) +
    ggplot2::scale_size_continuous(range = c(0.001, 1))

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

