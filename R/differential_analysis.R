#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom MASS glm.nb
#' @importFrom matrixTests row_kruskalwallis
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @title Compare exposures of annotated samples
#' @description \code{exposure_differential_analysis} is used to run
#' differential analysis on the signature exposures of annotated samples within
#' the \code{\linkS4class{musica_result}} object.
#' @param musica_result A \code{\linkS4class{musica_result}} object
#' @param annotation Column in the sample_annotations table of the
#' \code{\linkS4class{musica_result}} object
#' @param method Any method in \code{c("wilcox", "kruskal", "glm.nb")}
#' used to perform differential analysis on signature exposures
#' @param group1 character vector used in the Wilcox test. Elements in
#' \code{group1} are compared to elements in \code{group2}. This is
#' required for \code{annotation} with more than 2 levels.
#' @param group2 character vector used in the Wilcox test. Elements in
#' \code{group2} are compared to elements in \code{group1}. This is
#' required for \code{annotation} with more than 2 levels.
#' @param ... Additional arguments to be passed to the chosen method
#' @return A matrix containing statistics summarizing the analysis dependent
#' on the chosen method
#' @examples
#' data("res_annot")
#' exposure_differential_analysis(res_annot, "Tumor_Subtypes", method="wilcox")
#' @export
exposure_differential_analysis <- function(musica_result, annotation,
                                  method = c("wilcox", "kruskal", "glm.nb"),
                                  group1 = NULL, group2 = NULL,
                                  ...) {
  method <- match.arg(method)
  if (!methods::is(musica_result, "musica_result")) {
    stop("Input to exposure_differential_analysis must be a musica_result
         object.")
  }
  annotations <- samp_annot(musica_result)[[annotation]]
  if (is.null(annotations)) {
    stop(annotation, " does not exist in musica_result.")
  }
  exposures <- exposures(musica_result)
  groups <- unique(annotations) %>% sort()
  annotations <- factor(annotations)
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
    header <- data.frame(x = group1, y = group2) %>%
      dplyr::mutate(p = paste0(.data$x, "-", .data$y, "(Pr(>|z|))"),
                    adj = paste0(.data$x, "-", .data$y, "(fdr)"))
    diff.out <- sapply(seq_len(length(group1)), FUN = function(i) {
      x <- exposures[, annotations == group1[i]] %>% as.matrix()
      y <- exposures[, annotations == group2[i]] %>% as.matrix()
      out <- matrixTests::row_wilcoxon_twosample(x, y, ...)$pvalue
      return(out)
      })
    p <- p.adjust(
      diff.out[, (ncol(diff.out) - length(group1) + 1):ncol(diff.out)],
      method = "BH") %>% matrix(ncol = length(group1), byrow = F)
    diff.out <- cbind(diff.out, p) %>% as.data.frame()
    colnames(diff.out) <- c(header$p, header$adj)
  } else if (method == "kruskal") {
    header <- c("K-W chi-squared", "df", "p-value", "fdr")
    diff.out <- matrixTests::row_kruskalwallis(exposures, annotations, ...) %>%
      dplyr::select(.data$statistic, .data$df, .data$pvalue)
    diff.out$fdr <- p.adjust(diff.out$pvalue, method = "BH")
    colnames(diff.out) <- header
  } else if (method == "glm.nb") {
    header <- data.frame(y = groups) %>%
      dplyr::mutate(coef = paste0(.data$y, "(coef)"),
             sd = paste0(.data$y, "(Std. Error)"),
             z = paste0(.data$y, "(z)"),
             p = paste0(.data$y, "(Pr(>|z|))"),
             adj = paste0(.data$y, "(fdr)"))
    diff.out <- apply(exposures, 1, FUN = function(y) {
      fit <- MASS::glm.nb(round(y) ~ annotations, ...)
      out <- c(summary(fit)$coefficients,
               anova(fit)["annotations", "Pr(>Chi)"])
      return(out)
    }) %>% t()
    anova.out <- data.frame(anova = diff.out[, ncol(diff.out)])
    anova.out$fdr <- p.adjust(diff.out[, ncol(diff.out)], method = "BH")
    p <- p.adjust(
      diff.out[, (ncol(diff.out) - length(groups) + 1):ncol(diff.out)],
                  method = "BH") %>% matrix(ncol = length(groups), byrow = F)
    diff.out <- cbind(diff.out[, -ncol(diff.out)], p)
    diff.out <- cbind(diff.out, anova.out)
    colnames(diff.out) <- c(header$coef, header$sd, header$z, header$p,
                            header$adj, "ANOVA(Pr(>Chi))", "ANOVA(fdr)")
  }
  rownames(diff.out) <- rownames(exposures)
  return(diff.out)
}

#' @title Compare exposures of annotated samples
#' @description \code{plot_differential_analysis} is used to plot
#' differential analysis created by \code{exposure_differential_analysis}.
#' @param analysis Analysis created by \code{exposure_differential_analysis}
#' @param analysis_type Currently only \code{"glm"} supported
#' @param samp_num Number of samples that went into the analysis
#' @examples
#' data("res_annot")
#' analysis <- exposure_differential_analysis(res_annot, "Tumor_Subtypes", 
#' method="wilcox")
#' plot_differential_analysis(analysis, "glm", 2)
#' @export
plot_differential_analysis <- function(analysis, analysis_type, samp_num) {
  if (analysis_type == "glm") {
    dt <- data.table::melt(data.table::setDT(analysis[, seq_len(samp_num)], 
                                             keep.rownames = TRUE), "rn")
    dt$signif <- ifelse(dt$value < 0.01, 1, 0)
    p <- ggplot2::ggplot(dt, aes_string(fill = "rn", y = "value", 
                                        x = "variable")) + 
      geom_bar(position="dodge", stat="identity")
    p <- .gg_default_theme(p)
    p <- p + theme(legend.title = element_blank())
    return(p)
  }
}