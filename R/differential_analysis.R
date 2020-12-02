#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom MASS glm.nb

#' @title Compare exposures of annotated samples
#' @description \code{differential_exposure} is used to run differential analysis on 
#' the signature exposures of annotated samples within the 
#' \code{\linkS4class{musica_result}} object.
#' @param musica_result A \code{\linkS4class{musica_result}} object 
#' @param annotation Column in the sample_annotations table of the 
#' \code{\linkS4class{musica_result}} object
#' @param method Any method in \code{c("wilcox", "kruskal", "glm.nb")} 
#' used to perform differential analysis on signature exposures
#' @param pair Paired annotations for use in Wilcox test if \code{annotation}
#' has more than 2 levels. This pair in \code{annotation}  will be compared.
#' @param ... Additional arguments to be passed to the chosen method
#' @return A matrix containing statistics summarizing the analysis dependent
#' on the chosen method
#' @examples 
#' data("res_annot")
#' differential_exposure(res_annot, "Tumor_Subtypes", method="wilcox")
#' @export
#' 
differential_exposure <- function(musica_result, annotation, 
                                  method=c("wilcox","kruskal", "glm.nb"),
                                  pair=NULL,...) {
  method <- match.arg(method)
  if (!methods::is(musica_result, "musica_result")) {
    stop("Input to differential_exposure must be a musica_result object.")
  }
  annotations <- samp_annot(musica_result)[[annotation]]
  if (is.null(annotations)) {
    stop(annotation," does not exist in musica_result.")
  }
  exposures <- exposures(musica_result)
  groups <- unique(annotations)
  # Set order of levels by the order specified in the user input.
  # This is done to preserve correct order during diff anal testing.
  # Check if already a factor
  annotations <- factor(annotations, levels = unique(groups))
  if (length(groups) < 2) {
    stop("annotation must have at least 2 unique values")
  }
  diff.out <- 0
  l <- length(exposures)
  if (method == "wilcox") {
    if (length(groups) > 2) {
      if (is.null(pair)) {
        stop("pair must be provided for annotations consisting of more than 2 unique elements")
      }
      if (length(factor(pair)) != 2) {
        stop("pair must be of length 2")
      }
    } else {
      pair <- groups
    }
    header <- c(paste0(pair[1],"(mean)"),paste0(pair[2],"(mean)"), 
                "p-value", "fdr")
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- wilcox.test(y[annotations == pair[1]],y[annotations == pair[2]],...)
      m1 <- mean(y[annotations==pair[1]])
      m2 <- mean(y[annotations==pair[2]])
      return(c(m1, m2, out$p.value))
    }) %>% t()
    p <- p.adjust(diff.out[,3], method="BH")
    diff.out <- cbind(diff.out, p)
    colnames(diff.out) <- header
  } else if (method=="kruskal") {
    header <- c("K-W chi-squared", "df", "p-value")
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- kruskal.test(y ~ annotations, ...)
      return (c(out$statistic, out$parameter, out$p.value))
    }) %>% t()
    colnames(diff.out) <- c(header)
  } else if (method=="glm.nb") {
    header <- data.frame(y=groups) %>%
      dplyr::mutate(coef=paste(.data$y,"(coef)", sep=""),
             sd=paste(.data$y,"(Std. Error)", sep=""),
             z=paste(.data$y,"(z)", sep=""),
             p=paste(.data$y, "(Pr(>|z|))", sep=""),
             adj=paste(.data$y, "(fdr)", sep=""))
    diff.out <- apply(exposures, 1, FUN=function(y) {
      fit <- MASS::glm.nb(round(y) ~ annotations, ...)
      out <- c(summary(fit)$coefficients, anova(fit)["annotations", "Pr(>Chi)"])
      return(out)
    }) %>% t()
    anova.out <- data.frame(anova=diff.out[, ncol(diff.out)])
    anova.out$fdr <- p.adjust(diff.out[, ncol(diff.out)], method="BH")
    p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
                  method="BH") %>% matrix(ncol=length(groups), byrow=F)
    diff.out <- cbind(diff.out[,-ncol(diff.out)], p)
    diff.out <- cbind(diff.out, anova.out)
    colnames(diff.out) <- c(header$coef, header$sd, header$z, header$p, 
                            header$adj, "ANOVA (Pr(>Chi))", "ANOVA (fdr)")
  } else {
    stop("Method is not supported. Please provide one of: 
         wilcox, kruskal, glm.nb")
  }
  
  return (diff.out)
}

