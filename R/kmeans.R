#' @title Perform clustering analysis from a musica result object
#' @description Proportional sample exposures will be used as input to perform
#' clustering.
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' discovery or prediction.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing desired model.
#' Default \code{"result"}.
#' @param nclust Pre-defined number of clusters.
#' @param proportional Logical, indicating if proportional exposure (default)
#' will be used for clustering.
#' @param method Clustering algorithms. Options are "kmeans" (K-means, with
#' \link[stats]{kmeans}, seeded by \code{seed}), "fanny" (fuzzy clustering, with
#' \link[cluster]{fanny}), "hkmeans" (hybrid of hierarchical K-means),
#' "hclust" (hierarchical clustering), "pam" (PAM), and "clara" (Clara).
#' @param dis.method Methods to calculate dissimilarity matrix. Options are
#' "euclidean" (default), "manhattan", "jaccard", "cosine" (1 - cosine
#' similarity), and "canberra". "kmeans" accepts only "euclidean",
#' "hkmeans" accepts "euclidean", "manhattan" and "canberra", and "clara"
#' accepts "euclidean" and "manhattan" (\link[cluster]{clara}'s "jaccard" can
#' crash R).
#' @param hc.method Methods to perform hierarchical clustering. Options are
#' "ward.D" (default), "ward.D2", "single", "complete", "average", "mcquitty",
#' "median", and "centroid".
#' @param clara.samples Number of samples to be drawn from dataset. Only used
#' when "clara" is selected. Default is 5.
#' @param iter.max Maximum number of iterations, for "kmeans", "hkmeans" and
#' "fanny". Default \code{NULL} uses 10 for "kmeans" and "hkmeans" and 500
#' for "fanny".
#' @param tol Tolerance for convergence, for "fanny" only.
#' @param nstart Number of random starts, for "kmeans" only. Default
#' \code{10}.
#' @param seed Seed for the random starts of "kmeans", so the clusters are the
#' same each time; the caller's random number state is left unchanged.
#' \code{NULL} uses the caller's random number state. Default \code{1}.
#' @return A one-column data frame with sample IDs as row names and cluster
#' number for each sample.
#' @seealso \link[stats]{kmeans}
#' @examples
#' set.seed(123)
#' data(res_annot)
#' clust_out <- cluster_exposure(res_annot, model_name = "res_annot",
#' nclust = 2)
#' @export
# Exported argument names predate snake_case; renaming them would break calls.
# nolint start: object_name_linter.
cluster_exposure <- function(musica, model_name, modality = "SBS96",
                             result_name = "result", nclust,
                             proportional = TRUE, method = "kmeans",
                             dis.method = "euclidean", hc.method = "ward.D",
                             clara.samples = 5, iter.max = NULL, tol = 1e-15,
                             nstart = 10, seed = 1) {
  # nolint end
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

  method <- match.arg(method, .cluster_methods)
  dis.method <- match.arg(dis.method, # nolint: object_name_linter.
                          .cluster_dis_methods)
  hc.method <- match.arg(hc.method, # nolint: object_name_linter.
                         .cluster_hc_methods)
  .check_cluster_distance(method, dis.method)
  # read exposure data
  expos <- exposures(result)
  if (isTRUE(proportional)) {
    expos <- t(.proportional_exposures(expos))
  } else {
    expos <- t(expos)
  }
  clusters <- .cluster_rows(expos, nclust,
    method = method, dis_method = dis.method,
    hc_method = hc.method, clara_samples = clara.samples,
    iter_max = iter.max, tol = tol, nstart = nstart, seed = seed
  )
  clust_out <- data.frame(cluster = clusters)
  return(clust_out)
}

# Choices shared by cluster_exposure() and k_select()
.cluster_methods <- c("kmeans", "fanny", "hkmeans", "hclust", "pam", "clara")
.cluster_dis_methods <- c("euclidean", "manhattan", "jaccard", "cosine",
                          "canberra")
.cluster_hc_methods <- c("ward.D", "ward.D2", "single", "complete",
                         "average", "mcquitty", "median", "centroid")

# Stops if a clustering method can't use the chosen dissimilarity
.check_cluster_distance <- function(method, dis_method) {
  if (method == "kmeans" && dis_method != "euclidean") {
    stop("For kmeans clustering, dis.method must be \"euclidean\"; ",
         "use \"pam\", \"fanny\" or \"hclust\" for other distances.")
  }
  if (method == "hkmeans" &&
        !dis_method %in% c("euclidean", "manhattan", "canberra")) {
    stop("For hkmeans clustering, please choose the following methods to ",
         "calculate dissimilarity matrix: euclidean, manhattan, canberra")
  }
  # cluster::clara(metric = "jaccard") can crash R (cluster 2.1.8.2)
  if (method == "clara" && !dis_method %in% c("euclidean", "manhattan")) {
    stop("For clara clustering, please choose the following methods to ",
         "calculate dissimilarity matrix: euclidean, manhattan. ",
         "Use \"pam\", \"fanny\" or \"hclust\" for other distances.")
  }
}

# Dissimilarity between the rows (samples) of expos
.exposure_dissimilarity <- function(expos, dis_method) {
  diss <- philentropy::distance(x = expos, method = dis_method,
                                use.row.names = TRUE, as.dist.obj = TRUE)
  # philentropy's "cosine" is a similarity
  if (dis_method == "cosine") {
    diss <- 1 - diss
    # Rounding can leave identical directions slightly below zero
    diss[diss < 0] <- 0
  }
  diss
}

# Evaluates expr without philentropy::distance()'s "Metric: ..." messages,
# which it prints even with mute.message = TRUE (philentropy 0.10.0)
.muffle_philentropy <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (startsWith(conditionMessage(m), "Metric:")) {
      invokeRestart("muffleMessage")
    }
  })
}

# Clusters the rows (samples) of expos into k clusters. Returns the cluster
# numbers, named by sample. The arguments are those of cluster_exposure().
.cluster_rows <- function(expos, k, method, dis_method, hc_method,
                          clara_samples = 5, iter_max = NULL, tol = 1e-15,
                          nstart = 10, seed = NULL) {
  if (method %in% c("fanny", "hclust", "pam")) {
    diss <- .exposure_dissimilarity(expos, dis_method)
  }
  if (method == "kmeans") {
    run_kmeans <- function() {
      stats::kmeans(x = expos, centers = k,
                    iter.max = if (is.null(iter_max)) 10 else iter_max,
                    nstart = nstart)
    }
    res <- if (is.null(seed)) {
      run_kmeans()
    } else {
      withr::with_seed(seed, run_kmeans())
    }
    res$cluster
  } else if (method == "fanny") {
    cluster::fanny(
      x = diss, k = k, diss = TRUE,
      maxit = if (is.null(iter_max)) 500 else iter_max,
      tol = tol
    )$clustering
  } else if (method == "hkmeans") {
    factoextra::hkmeans(
      x = expos, k = k, hc.metric = dis_method,
      hc.method = hc_method,
      iter.max = if (is.null(iter_max)) 10 else iter_max
    )$cluster
  } else if (method == "hclust") {
    factoextra::hcut(x = diss, k = k, isdiss = TRUE,
                     hc_method = hc_method)$cluster
  } else if (method == "pam") {
    cluster::pam(x = diss, k = k, diss = TRUE, cluster.only = TRUE)
  } else {
    cluster::clara(x = expos, k = k, metric = dis_method,
                   samples = clara_samples)$clustering
  }
}

#' @title Visualize clustering results
#' @description The clustering results can be visualized on a UMAP panel.
#' Three different types of plots can be generated using this function:
#' cluster-by-signature plot, cluster-by-annotation plot, and a single UMAP
#' plot.
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' discovery or prediction. A two-dimensional UMAP has to be stored in this
#' object.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing desired model.
#' Default \code{"result"}.
#' @param clusters The result generated from cluster_exposure function.
#' @param group A single character string indicating the grouping factor.
#' Possible options are: "signature" (columns are signatures in a grid),
#' "annotation" (columns are sample annotation), and "none" (a single UMAP
#' plot). Default is "signature".
#' @param annotation Column name of annotation.
#' @param plotly If TRUE, the plot will be made interactive using plotly.
#' @return Generate a ggplot or plotly object.
#' @seealso \link{create_umap}
#' @examples
#' set.seed(123)
#' data(res_annot)
#' # Get clustering result
#' clust_out <- cluster_exposure(
#'   musica = res_annot, model_name = "res_annot",
#'   nclust = 2, iter.max = 15
#' )
#' # UMAP
#' create_umap(musica = res_annot, model_name = "res_annot")
#' # generate cluster X signature plot
#' plot_cluster(
#'   musica = res_annot, model_name = "res_annot",
#'   clusters = clust_out, group = "signature"
#' )
#' # generate cluster X annotation plot
#' plot_cluster(
#'   musica = res_annot, model_name = "res_annot",
#'   clusters = clust_out, group = "annotation",
#'   annotation = "Tumor_Subtypes"
#' )
#' # generate a single UMAP plot
#' plot_cluster(
#'   musica = res_annot, model_name = "res_annot",
#'   clusters = clust_out, group = "none"
#' )
#' @export

plot_cluster <- function(musica, model_name, modality = "SBS96",
                         result_name = "result", clusters, group = "signature",
                         annotation = NULL, plotly = TRUE) {
  # dummy variables
  UMAP_1 <- NULL # nolint: object_name_linter.
  UMAP_2 <- NULL # nolint: object_name_linter.
  exposure <- NULL

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

  group <- match.arg(group, c("signature", "annotation", "none"))
  if (length(umap(result)) == 0) {
    stop("UMAP not found in musica_result object. 
         Run create_umap(", deparse(substitute(result)), ") first.")
  }
  k_toplot <- cbind(result@umap, clusters)
  k_toplot$cluster <- factor(k_toplot$cluster)
  if (group == "signature") {
    expos <- t(.proportional_exposures(exposures(result)))
    clust_by_sigs <- cbind(k_toplot, expos) |>
      tibble::rownames_to_column(var = "sample") |>
      tidyr::pivot_longer(
        cols = colnames(expos),
        names_to = "signature",
        values_to = "exposure",
        names_repair = "minimal"
      )
    p <- ggplot2::ggplot(clust_by_sigs,
                         aes(x = UMAP_1, y = UMAP_2, colour = exposure)) +
      geom_point() +
      facet_grid(cluster ~ signature) +
      ggplot2::scale_colour_gradientn(colors =
                                        c("blue", "green", "yellow", "orange",
                                          "red"), name = "Fraction") +
      ggplot2::theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )
  } else if (group == "annotation") {
    if (is.null(annotation) || !annotation %in% colnames(samp_annot(musica))) {
      stop("Sample annotation not found or invalid annotation column name.")
    } else {
      # Match the chosen annotation to the UMAP rows by sample name
      annot_values <- .annot_for_samples(musica, annotation,
                                         rownames(k_toplot))
      clust_by_annot <- cbind(k_toplot,
                              annotation = factor(annot_values))
      p <- ggplot2::ggplot(clust_by_annot,
                           aes(x = .data$UMAP_1, y = .data$UMAP_2,
                               colour = .data$cluster)) +
        geom_point() +
        facet_grid(cluster ~ annotation) +
        ggplot2::theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        )
    }
  } else {
    p <- ggplot2::ggplot(k_toplot, aes(x = .data$UMAP_1, y = .data$UMAP_2,
                                       colour = .data$cluster)) +
      geom_point() +
      ggplot2::theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )
  }
  if (plotly) {
    p <- plotly::ggplotly(p)
    return(p)
  } else {
    return(p)
  }
}

#' @title Plots for helping decide number of clusters
#' @description To help decide the number of cluster, three different methods
#' are provided: total within cluster sum of squares, average silhouette
#' coefficient, and gap statistics.
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' discovery or prediction. A two-dimensional UMAP has to be stored in this
#' object.
#' @param model_name The name of the desired model.
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing desired model.
#' Default \code{"result"}.
#' @param method A single character string indicating which statistic to use for
#' plot. Options are "wss" (total within cluster sum of squares), "silhouette"
#' (average silhouette coefficient), and "gap_stat" (gap statistic). Default is
#' "wss".
#' @param clust.method A character string indicating clustering method. Options
#' are "kmeans" (default), "fanny", "hkmeans", "hclust" (hierarchical
#' clustering), "pam", and "clara". The samples are clustered as
#' \code{\link{cluster_exposure}} clusters them, with its default settings
#' for the other arguments.
#' @param n An integer indicating maximum number of clusters to test. Default is
#' 10. It is lowered, with a message, to the number of samples minus one, or
#' for "fanny" to half the number of samples minus one. At least three samples
#' (six for "fanny") are needed.
#' @param proportional Logical, indicating if proportional exposure (default)
#' will be used for clustering.
#' @param seed Seed for the random steps (k-means starts, and the bootstrap of
#' "gap_stat"), so the plot is the same each time; the caller's random number
#' state is left unchanged. \code{NULL} uses the caller's random number state.
#' Default \code{1}.
#' @param dis.method Method to calculate the dissimilarity between samples,
#' as in \code{\link{cluster_exposure}}. Options are "euclidean" (default),
#' "manhattan", "jaccard", "cosine" (1 - cosine similarity), and "canberra".
#' "kmeans" accepts only "euclidean", "hkmeans" accepts "euclidean",
#' "manhattan" and "canberra", and "clara" accepts "euclidean" and
#' "manhattan". "wss" and "silhouette" are computed from this
#' dissimilarity (for "hkmeans", with \link[stats]{dist}, as
#' \link[factoextra]{hkmeans} uses). "gap_stat" clusters with it, but \link[cluster]{clusGap}
#' measures the spread within clusters with Euclidean distance. With
#' "gap_stat", "fanny" also clusters the data \link[cluster]{clusGap}
#' simulates, and may warn that it has not converged on some of it.
#' @param hc.method Method of hierarchical clustering, for "hclust" and
#' "hkmeans", as in \code{\link{cluster_exposure}}. Options are "ward.D"
#' (default), "ward.D2", "single", "complete", "average", "mcquitty",
#' "median", and "centroid".
#' @return A ggplot object.
#' @seealso \link[factoextra]{fviz_nbclust}
#' @examples
#' data(res_annot)
#' set.seed(123)
#' # Make an elbow plot
#' k_select(res_annot, model_name = "res_annot", method = "wss", n = 6)
#' # Plot average silhouette coefficient against number of clusters
#' k_select(res_annot, model_name = "res_annot", method = "silhouette", n = 6)
#' # Plot gap statistics against number of clusters
#' k_select(res_annot, model_name = "res_annot", method = "gap_stat", n = 6)
#' # Choose k for hierarchical clustering on 1 - cosine similarity
#' k_select(res_annot, model_name = "res_annot", method = "silhouette",
#'   clust.method = "hclust", dis.method = "cosine", n = 6)
#' @export
# nolint start: object_name_linter.
k_select <- function(musica, model_name, modality = "SBS96",
                     result_name = "result", method = "wss",
                     clust.method = "kmeans", n = 10, proportional = TRUE,
                     seed = 1, dis.method = "euclidean", hc.method = "ward.D") {
  # nolint end
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

  method <- match.arg(method, c("wss", "silhouette", "gap_stat"))
  clust.method <- match.arg(clust.method, # nolint: object_name_linter.
                            .cluster_methods)
  dis.method <- match.arg(dis.method, # nolint: object_name_linter.
                          .cluster_dis_methods)
  hc.method <- match.arg(hc.method, # nolint: object_name_linter.
                         .cluster_hc_methods)
  .check_cluster_distance(clust.method, dis.method)
  expos <- exposures(result)
  if (isTRUE(proportional)) {
    expos <- t(.proportional_exposures(expos))
  } else {
    expos <- t(expos)
  }
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 2) {
    stop("'n' must be a single number of at least 2.")
  }
  n_samples <- nrow(expos)
  if (n_samples < 3) {
    stop("k_select() needs at least three samples; the model has ",
         n_samples, ".")
  }
  if (clust.method == "fanny") {
    max_k <- n_samples %/% 2 - 1
    if (max_k < 2) {
      stop("k_select() with clust.method = \"fanny\" needs at least six ",
           "samples; the model has ", n_samples, ".")
    }
    if (n > max_k) {
      message("Testing up to ", max_k, " clusters, the most fanny allows ",
              "for ", n_samples, " samples, instead of n = ", n, ".")
      n <- max_k
    }
  } else if (n > n_samples - 1) {
    message("Testing up to ", n_samples - 1, " clusters, one fewer than the ",
            n_samples, " samples, instead of n = ", n, ".")
    n <- n_samples - 1
  }
  # Cluster as cluster_exposure() does
  cluster_fun <- function(x, k) {
    clusters <- .muffle_philentropy(.cluster_rows(x, k,
      method = clust.method, dis_method = dis.method,
      hc_method = hc.method
    ))
    list(cluster = clusters)
  }
  run_nbclust <- function() {
    if (method == "gap_stat") {
      # clusGap() clusters simulated data, so diss can't be given
      return(factoextra::fviz_nbclust(expos, cluster_fun, method = method,
                                      k.max = n))
    }
    # Score with the distance the clusters come from; hkmeans builds its
    # tree with stats::dist(), whose canberra differs from philentropy's
    diss <- if (clust.method == "hkmeans") {
      stats::dist(expos, method = dis.method)
    } else {
      .muffle_philentropy(.exposure_dissimilarity(expos, dis.method))
    }
    factoextra::fviz_nbclust(expos, cluster_fun, method = method,
                             diss = diss, k.max = n)
  }
  if (is.null(seed)) {
    run_nbclust()
  } else {
    withr::with_seed(seed, run_nbclust())
  }
}
