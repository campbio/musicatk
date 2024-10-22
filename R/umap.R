
#' @title Create a UMAP from a model result
#' @description Proportional sample exposures will be used as input into the 
#' \code{\link[uwot]{umap}} function to generate a two dimensional UMAP.
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' signature discovery or prediction.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing the model.
#' Default \code{"result"}.
#' @param n_neighbors The size of local neighborhood used for views of
#' manifold approximation. Larger values result in more global the manifold,
#' while smaller values result in more local data being preserved.
#' If \code{n_neighbors} is larger than the number of samples,
#' then \code{n_neighbors} will automatically be set to the number of samples
#' in the \code{\linkS4class{musica}}. Default \code{30}.
#' @param min_dist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where
#' nearby points on the manifold are drawn closer together, while larger
#' values will result on a more even dispersal of points. Default \code{0.2}.
#' @param spread The effective scale of embedded points. In combination with
#' ‘min_dist’, this determines how clustered/clumped the embedded points are.
#' Default \code{1}.
#' @return A \code{\linkS4class{musica}} object with a new UMAP
#' stored in the \code{UMAP} slot of the \code{\linkS4class{result_model}}
#' object for the model.
#' @seealso See \link{plot_umap} to display the UMAP and
#' \code{\link[uwot]{umap}} for more information on the individual parameters
#' for generating UMAPs. 
#' @examples
#' data(res_annot)
#' create_umap(res_annot, model_name = "res_annot")
#' @export
create_umap <- function(musica, model_name, modality = "SBS96", result_name = "result",
                        n_neighbors = 30, min_dist = 0.75, spread = 1) {
  
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
  
  # dummy call to Matrix
  Matrix::Matrix()
  
  samples <- exposures(result)
  samples <- sweep(samples, 2, colSums(samples), FUN = "/")
  
  # n_neighbors cannot be bigger than the total number of samples
  if (n_neighbors > ncol(samples)) {
    n_neighbors <- ncol(samples)
    message("The parameter 'n_neighbors' cannot be bigger than the total ",
            "number of samples. Setting 'n_neighbors' to ", n_neighbors, ".")
  }
  
  umap_out <- uwot::umap(t(samples), n_neighbors = n_neighbors, min_dist = 
                           min_dist, spread = spread, n_threads = 1,
                         n_sgd_threads = 1, pca = NULL,
                         metric = "cosine")
  rownames(umap_out) <- colnames(samples)
  colnames(umap_out) <- c("UMAP_1", "UMAP_2")
  eval.parent(substitute(umap(musica, result_name, modality, model_name) <- umap_out))
}

#' @title Plot a UMAP from a musica result
#' @description Plots samples on a UMAP scatterplot. Samples can be colored by
#' the levels of mutational signatures or by a annotation variable. 
#'
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' signature discovery or prediction.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing the model.Default
#' \code{"result"}.
#' @param color_by One of \code{"signatures"}, \code{"annotation"}, or
#' \code{"none"}. If \code{"signatures"}, then one UMAP scatterplot will be
#' generated for each signature and points will be colored by the level of 
#' that signature in each sample. If \code{annotation}, a single UMAP will 
#' be generated colored by the annotation selected using the parameter 
#' \code{annotation}. If \code{"none"}, a single UMAP scatterplot will be
#' generated with no coloring. Default \code{"signature"}.
#' @param proportional If \code{TRUE}, then the exposures will be normalized
#' to between 0 and 1 by dividing by the total number of counts for each sample.
#' Default \code{TRUE}.
#' @param annotation Sample annotation used to color the points. One used
#' when \code{color_by = "annotation"}. Default \code{NULL}.
#' @param point_size Scatter plot point size. Default \code{0.7}.
#' @param same_scale If \code{TRUE}, then all points will share the same color
#' scale in each signature subplot. If \code{FALSE}, then each signature subplot
#' will be colored by a different scale with different maximum values. Only
#' used when \code{color_by = "signature"}. Setting to \code{FALSE} is most
#' useful when the maximum value of various signatures are vastly different 
#' from one another. Default \code{TRUE}.
#' @param add_annotation_labels If \code{TRUE}, labels for each group in the
#' \code{annotation} variable will be displayed. Only used if
#' \code{color_by = "annotation"}. This not recommended if the annotation is
#' a continuous variable. The label is plotting using the centriod of each
#' group within the \code{annotation} variable. Default \code{FALSE}.
#' @param annotation_label_size Size of annotation labels. Only used if 
#' \code{color_by = "annotation"} and \code{add_annotation_labels = TRUE}.
#' Default \code{3}.
#' @param annotation_text_box Place a white box around the annotation labels
#' to improve readability. Only used if \code{color_by = "annotation"} and
#' \code{add_annotation_labels = TRUE}. Default \code{TRUE}.
#' @param plotly If \code{TRUE}, the the plot will be made interactive
#' using \code{\link[plotly]{plotly}}. Not used if \code{color_by = "signature"}
#' and \code{same_scale = FALSE}. Default \code{FALSE}.
#' @param clust Add cluster labels as annotation
#' @param legend Plot legend
#' @param strip_axes Remove axes labels for cleaner looking plots
#' @return Generates a ggplot or plotly object
#' @seealso See \link{create_umap} to generate a UMAP in a musica result.
#' @examples
#' data(res_annot)
#' create_umap(res_annot, "res_annot")
#' plot_umap(res_annot, "res_annot", color_by = "none")
#' @export
plot_umap <- function(musica, model_name, modality = "SBS96", 
                      result_name = "result", 
                      color_by = c("signatures", "annotation", 
                                           "cluster", "none"),
                      proportional = TRUE, annotation = NULL,
                      point_size = 0.7, same_scale = TRUE,
                      add_annotation_labels = FALSE,
                      annotation_label_size = 3,
                      annotation_text_box = TRUE, plotly = FALSE,
                      clust = NULL,
                      legend = TRUE,
                      strip_axes = FALSE) {
  
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
  
  umap <- umap(result)
  color_by <- match.arg(color_by)
  
  if (color_by %in% c("annotation", "cluster")) {
    if (is.null(annotation)) {
      stop("If the parameter or 'color_by' are is to 'annotation', ",
           "then the 'annotation' parameter must be supplied.")
    }
    
    # Add annotation to data.frame
    annot <- samp_annot(musica)
    
    # Manually override annotations with cluster labels
    if (color_by == "cluster") {
      annot <- clust
      annotation <- "cluster"
    }
    
    umap <- umap %>% as.data.frame %>%
      tibble::rownames_to_column(var = "sample")
    umap <- .add_annotation_to_df(musica, plot_dat = umap,
                                  annotation = annotation, 
                                  clust = clust)
    
    # Create base ggplot object
    if (plotly) {
      p <- ggplot(umap, aes_string(x = "UMAP_1", y = "UMAP_2",
                                   color = "annotation", text = "sample", 
                                   label = "annotation"))
    } else {
      p <- ggplot(umap, aes_string(x = "UMAP_1", y = "UMAP_2",
                                   color = "annotation"))
    }
    p <- p + geom_point(size = point_size) 
    
    # Add labels for discrete annotation groups
    if (isTRUE(add_annotation_labels) & !plotly) {
      centroid_list <- lapply(unique(umap$annotation), function(x) {
        df_sub <- umap[umap$annotation == x, ]
        median_1 <- stats::median(df_sub[, "UMAP_1"])
        median_2 <- stats::median(df_sub[, "UMAP_2"])
        cbind(median_1, median_2, x)
      })
      centroid <- do.call(rbind, centroid_list)
      centroid <- data.frame(
        UMAP_1 = as.numeric(centroid[, 1]),
        UMAP_2 = as.numeric(centroid[, 2]),
        annotation = centroid[, 3]
      )
      p <- p + ggplot2::geom_point(data = centroid, mapping =
          ggplot2::aes_string(x = "UMAP_1", y = "UMAP_1"),
          size = 0, alpha = 0) 
      
      if (isTRUE(annotation_text_box)) {
        p <- p + ggrepel::geom_label_repel(data = centroid, mapping =
                                          ggplot2::aes_string(label = "annotation"),
                                           size = annotation_label_size,
                                          show.legend = FALSE)
      } else {
        p <- p + ggrepel::geom_text_repel(data = centroid, mapping =
                                            ggplot2::aes_string(label = "annotation"),
                                          size = annotation_label_size,
                                          show.legend = FALSE)
      }
    }
    
    # Add theme and title 
    p <- .gg_default_theme(p)
    p <- p + ggplot2::scale_color_discrete(name = annotation)
    
  } else if (color_by == "signatures") {
    cols <- c("blue", "green", "yellow", "orange", "red")

    # Create df in long format with signatures and UMAP coords
    exposures <- exposures(result)
    if (isTRUE(proportional)) {
      color_lab <- "Fraction"
      exposures <- sweep(exposures, 2, colSums(exposures), FUN = "/")
    } else {
      color_lab <- "Counts"
    }

    cbind(umap, t(exposures)) %>%
      as.data.frame %>%
      tibble::rownames_to_column(var = "sample") %>%
      tidyr::pivot_longer(cols = rownames(exposures),
                          names_to = "signature",
                          values_to = "exposure",
                          names_repair = "minimal") -> df
    
    # Ensure that the signature order will be the same as in the exposures
    df$signature <- factor(df$signature, 
                           levels = gtools::mixedsort(rownames(exposures)))
    
    # Create base ggplot object for "signatures"
    if (isTRUE(same_scale)) {
      # Uses facet_grid to plot all signatures 
      breaks <- signif(seq(0, max(df$exposure), length.out = 5), 2)
      limits <- c(0, max(breaks))
      p <- ggplot(df, aes_string(x = "UMAP_1", y = "UMAP_2", 
                                 colour = "exposure")) +
        ggplot2::facet_wrap(~ signature, drop = TRUE, scales = "free") +
        geom_point() +
        ggplot2::scale_colour_gradientn(colors = cols, name = color_lab,
                                        breaks = breaks, limits = limits)
      p <- .gg_default_theme(p)
    } else {
      # Uses grid.arrange to plot many ggplot objects
      out <- by(data = df, INDICES = df$signature, FUN = function(m) {
        m <- droplevels(m)
        breaks <- signif(seq(0, max(m$exposure), length.out = 5), 2)
        limits <- c(0, max(breaks))
        m <- ggplot(m, aes_string(x = "UMAP_1",
                                  y = "UMAP_2",
                                  color = "exposure")) +
          ggplot2::geom_point() + 
          ggplot2::scale_colour_gradientn(colors = cols, name = color_lab,
                                          breaks = breaks, limits = limits) +
          ggplot2::ggtitle(unique(m$signature))
        m <- .gg_default_theme(m)
      })
      
      # Need to use invisible since this is a gtable object
      return(invisible(do.call(gridExtra::grid.arrange, out)))
    }
  } else {
    # Create base ggplot object without any color
    p <- ggplot(as.data.frame(umap), aes_string(x = "UMAP_1", y = "UMAP_2")) +
      geom_point(size = point_size)
    p <- .gg_default_theme(p)
  }
  
  if (!isTRUE(legend)) {
    p <- p + theme(legend.position = "none")
  }
  
  if (isTRUE(strip_axes)) {
    p <- p + theme(axis.title.x=element_blank(), 
                   axis.text.x=element_blank(), 
                   axis.ticks.x=element_blank(), 
                   axis.title.y=element_blank(), 
                   axis.text.y=element_blank(), 
                   axis.ticks.y=element_blank())
  }
  
  # Toggle plotly
  if (isTRUE(plotly)) {
    if (color_by == "annotation") {
      p <- plotly::ggplotly(p, tooltip = c("text", "x", "y", "type", "label"))
    } else {
      p <- plotly::ggplotly(p) 
    }
  }
  
  # Return all other scenarios other than
  # color_by = "signatures" & same_scale = TRUE
  return(p)
}
