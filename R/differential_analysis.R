#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom MASS glm.nb
#' @importFrom matrixTests row_kruskalwallis
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @title Compare exposures of annotated samples
#' @description \code{exposure_differential_analysis} is used to run
#' differential analysis on the signature exposures of annotated samples.
#' @param musica A \code{\linkS4class{musica}} object.
#' @param model_name The name of the model.
#' @param annotation Column in the sample_annotations table of the
#' \code{\linkS4class{musica}} object
#' @param modality The modality. Must be "SBS96", "DBS78", or "IND83". Default
#' \code{"SBS96"}.
#' @param result_name Name of the result list entry. Default \code{"result"}.
#' @param method Any method in \code{c("wilcox", "kruskal", "glm.nb")}
#' used to perform differential analysis on signature exposures. For
#' \code{"glm.nb"}, the first level of \code{annotation} (the first in sorted
#' order for a character column) is the reference group.
#' @param group1 character vector used in the Wilcox test. Elements in
#' \code{group1} are compared to elements in \code{group2}. This is
#' required for \code{annotation} with more than 2 levels.
#' @param group2 character vector used in the Wilcox test. Elements in
#' \code{group2} are compared to elements in \code{group1}. This is
#' required for \code{annotation} with more than 2 levels.
#' @param ... Additional arguments to be passed to the chosen method
#' @return A data frame whose columns depend on the method. P-values are
#' adjusted with the Benjamini-Hochberg method (\code{fdr} columns).
#' \describe{
#'   \item{\code{"wilcox"}}{One row per signature and comparison, with columns
#'   \code{Signature}, \code{Group1}, \code{Group2}, \code{p-value} and
#'   \code{fdr}. The FDR is adjusted within each comparison, across
#'   signatures.}
#'   \item{\code{"kruskal"}}{One row per signature, with columns
#'   \code{K-W chi-squared}, \code{df}, \code{p-value} and \code{fdr}.}
#'   \item{\code{"glm.nb"}}{One row per signature. \code{Intercept(coef)},
#'   \code{Intercept(Std. Error)}, \code{Intercept(z)} and
#'   \code{Intercept(Pr(>|z|))} describe the reference group. For each other
#'   group, \code{<group>(coef)} is its log fold change against the reference,
#'   with \code{<group>(Std. Error)}, \code{<group>(z)},
#'   \code{<group>(Pr(>|z|))} and \code{<group>(fdr)}, adjusted across
#'   signatures. \code{ANOVA(Pr(>Chi))} and \code{ANOVA(fdr)} test the
#'   annotation as a whole.}
#' }
#' @examples
#' data("res_annot")
#' exposure_differential_analysis(res_annot,
#'   model_name = "res_annot",
#'   annotation = "Tumor_Subtypes", method = "wilcox"
#' )
#' @export
exposure_differential_analysis <- function(musica, model_name, annotation,
                                           modality = "SBS96",
                                           result_name = "result",
                                           method = c("wilcox", "kruskal",
                                                      "glm.nb"),
                                           group1 = NULL, group2 = NULL, ...) {
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

  # dummy variables
  statistic <- NULL
  df <- NULL
  pvalue <- NULL

  method <- match.arg(method)
  model <- get_model(musica, result_name, modality, model_name)

  exposures <- exposures(model)
  annotations <- .annot_for_samples(musica, annotation, colnames(exposures))
  # The groups follow the factor levels the models use
  annotations <- factor(annotations)
  groups <- levels(annotations)
  if (length(groups) < 2) {
    stop("annotation must have at least 2 unique values")
  }
  if (!is.null(c(group1, group2))) {
    if (method != "wilcox" || length(groups) == 2) {
      message("'annotations' is of length 2. 'group1' and 'group2' were
                ignored.")
    }
  }
  diff.out <- 0
  if (method == "wilcox") {
    if (length(groups) > 2) {
      if (is.null(group1) || is.null(group2)) {
        stop("Parameters 'group1' and 'group2' are required for annotations
        containing more than 2 unique elements.")
      }
      if (!is.character(group1) || !is.character(group2)) {
        stop("Parameter 'group1' and 'group2' must be character vectors.")
      }
      if (length(group1) != length(group2)) {
        stop("'group1' and 'group2' must be the same length.")
      }
      if (any(group1 == group2)) {
        stop("All pairs of 'group1' and 'group2' must be unique")
      }
      if (!all(group1 %in% annotations)) {
        stop("'group1' does not exist in annotations.")
      }
      if (!all(group2 %in% annotations)) {
        stop("'group2' does not exist in annotations.")
      }
    } else {
      group1 <- groups[1]
      group2 <- groups[2]
    }
    diff.out <- lapply(seq_along(group1), FUN = function(i) {
      x <- exposures[, annotations == group1[i], drop = FALSE]
      y <- exposures[, annotations == group2[i], drop = FALSE]
      p <- matrixTests::row_wilcoxon_twosample(x, y, ...)$pvalue
      # BH within each comparison, across signatures
      data.frame(
        Signature = rownames(exposures),
        Group1 = as.character(group1[i]),
        Group2 = as.character(group2[i]),
        `p-value` = p,
        fdr = p.adjust(p, method = "BH"),
        check.names = FALSE
      )
    })
    diff.out <- do.call(rbind, diff.out)
  } else if (method == "kruskal") {
    header <- c("K-W chi-squared", "df", "p-value", "fdr")
    diff.out <- matrixTests::row_kruskalwallis(exposures, annotations, ...) %>%
      dplyr::select(statistic, df, pvalue)
    diff.out$fdr <- p.adjust(diff.out$pvalue, method = "BH")
    colnames(diff.out) <- header
    rownames(diff.out) <- rownames(exposures)
  } else if (method == "glm.nb") {
    reserved <- intersect(groups, c("Intercept", "ANOVA"))
    if (length(reserved) > 0) {
      stop(
        "For glm.nb, annotation values can't be named ",
        paste0("'", reserved, "'", collapse = " or "),
        ", which label the intercept and ANOVA columns."
      )
    }
    # Coefficients are read by name, so an aliased term gives NA
    terms <- c("(Intercept)", paste0("annotations", groups[-1]))
    labels <- c("Intercept", groups[-1])
    coef_cols <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    fits <- apply(exposures, 1, FUN = function(y) {
      fit <- MASS::glm.nb(round(y) ~ annotations, ...)
      coefs <- summary(fit)$coefficients
      coefs <- coefs[match(terms, rownames(coefs)), coef_cols, drop = FALSE]
      # Column-major: every term's estimate, then standard error, z and p
      c(as.vector(coefs), anova(fit)["annotations", "Pr(>Chi)"])
    })
    fits <- matrix(fits, nrow = nrow(exposures), byrow = TRUE)
    n_terms <- length(terms)
    group_p <- fits[, 3 * n_terms + seq_len(n_terms)[-1], drop = FALSE]
    group_fdr <- apply(group_p, 2, p.adjust, method = "BH")
    group_fdr <- matrix(group_fdr, nrow = nrow(group_p))
    anova_p <- fits[, ncol(fits)]
    diff.out <- data.frame(
      fits[, -ncol(fits), drop = FALSE], group_fdr, anova_p,
      p.adjust(anova_p, method = "BH")
    )
    colnames(diff.out) <- c(
      paste0(labels, "(coef)"), paste0(labels, "(Std. Error)"),
      paste0(labels, "(z)"), paste0(labels, "(Pr(>|z|))"),
      paste0(groups[-1], "(fdr)"), "ANOVA(Pr(>Chi))", "ANOVA(fdr)"
    )
    rownames(diff.out) <- rownames(exposures)
  }
  return(diff.out)
}


#' @title Plot a differential analysis of exposures
#' @description \code{plot_differential_analysis} plots the false discovery
#' rate of each signature from a table made by
#' \code{\link{exposure_differential_analysis}}, as \eqn{-log_{10}(FDR)} bars.
#' The method is read from the table's columns:
#' \describe{
#'   \item{\code{"wilcox"}}{one bar per comparison, labelled
#'   \code{"<Group1> vs <Group2>"};}
#'   \item{\code{"kruskal"}}{one bar per signature;}
#'   \item{\code{"glm.nb"}}{one bar per non-reference group, against the
#'   reference group, and one for the ANOVA.}
#' }
#' An FDR of zero is plotted at the smallest positive double.
#' @param analysis A table made by \code{exposure_differential_analysis}.
#' @param analysis_type Ignored, with a warning. The method is read from
#' \code{analysis}.
#' @param samp_num Ignored, with a warning.
#' @param threshold The FDR drawn as a dashed line. Default \code{0.05}.
#' @return A ggplot object.
#' @examples
#' data("res_annot")
#' analysis <- exposure_differential_analysis(res_annot,
#'   model_name = "res_annot",
#'   annotation = "Tumor_Subtypes", method = "wilcox"
#' )
#' plot_differential_analysis(analysis)
#' @export
plot_differential_analysis <- function(analysis, analysis_type = NULL,
                                       samp_num = NULL, threshold = 0.05) {
  if (!is.null(analysis_type) || !is.null(samp_num)) {
    warning(
      "'analysis_type' and 'samp_num' are ignored; the method is read ",
      "from the columns of 'analysis'."
    )
  }
  valid_threshold <- is.numeric(threshold) && length(threshold) == 1 &&
    !is.na(threshold) && threshold > 0 && threshold <= 1
  if (!valid_threshold) {
    stop("'threshold' must be a single number in (0, 1].")
  }
  dt <- .diff_plot_data(analysis)
  dt$value <- -log10(pmax(dt$fdr, .Machine$double.xmin))
  p <- ggplot2::ggplot(dt, aes(
    x = .data$signature, y = .data$value,
    fill = .data$comparison
  )) +
    geom_col(position = position_dodge()) +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    scale_x_discrete(limits = unique(dt$signature)) +
    labs(x = NULL, y = expression(-log[10](FDR)))
  p <- .gg_default_theme(p)
  p <- p + theme(legend.title = element_blank())
  return(p)
}

# One row per signature and comparison: signature, comparison, fdr
.diff_plot_data <- function(analysis) {
  cols <- colnames(analysis)
  if (!is.data.frame(analysis)) {
    cols <- NULL
  }
  if (all(c("Signature", "Group1", "Group2", "fdr") %in% cols)) {
    out <- data.frame(
      signature = as.character(analysis$Signature),
      comparison = paste(analysis$Group1, "vs", analysis$Group2),
      fdr = analysis$fdr
    )
  } else if (all(c("K-W chi-squared", "fdr") %in% cols)) {
    out <- data.frame(
      signature = rownames(analysis), comparison = "Kruskal-Wallis",
      fdr = analysis$fdr
    )
  } else if ("ANOVA(fdr)" %in% cols) {
    fdr_cols <- grep("\\(fdr\\)$", cols, value = TRUE)
    out <- data.frame(
      signature = rep(rownames(analysis), length(fdr_cols)),
      comparison = rep(sub("\\(fdr\\)$", "", fdr_cols),
                       each = nrow(analysis)),
      fdr = unlist(analysis[fdr_cols], use.names = FALSE)
    )
  } else {
    stop(
      "'analysis' must be a table returned by ",
      "exposure_differential_analysis()."
    )
  }
  out$comparison <- factor(out$comparison, levels = unique(out$comparison))
  out
}
