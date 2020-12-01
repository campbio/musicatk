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
#' @param groups Paired annotations for use in Wilcox test if \code{annotation}
#' has more than 2 levels.
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
                                  groups=NULL,...) {
  method <- match.arg(method)
  if (!methods::is(musica_result, "musica_result")) {
    stop("Input to differential_exposure must be a musica_result object.")
  }
  #samp_annot
  annotations <- 
    musica_result@musica@sample_annotations[[annotation]]
  if (is.null(annotations)) {
    stop(annotation," does not exist in musica_result.")
  }
  exposures <- exposures(musica_result)
  groups <- unique(annotations)

  # Set order of levels by the order specified in the user input.
  # This is done to preserve correct order during diff anal testing.
  # Check if already a factor
  annotations <- factor(annotations, levels = unique(groups))
  diff.out <- 0
  l <- length(exposures)
  #Change wilcox to only take 2 pairs
  if (method == "wilcox" || is.null(method)) {
    pairs <- combn(groups,2) %>% t()
    header <- #format_differential_exposures(paste(pairs[,1],pairs[,2],sep="-"), c("(mean)","(p-value)"))
    data.frame(y=pairs[,1], x=pairs[,2]) %>%
      dplyr::mutate(c=paste(.data$y,"-",.data$x,"(W)", sep=""),
             p=paste(.data$y,"-",.data$x,"(p-value)", sep=""),
             f=paste(.data$y,"-",.data$x,"(fdr)", sep=""))

    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- apply(pairs, 1, FUN=function(p) {
        # Add means 
        out <- wilcox.test(y[annotations==p[1]],y[annotations==p[2]],...)
        return (c(s=out$statistic, p=out$p.value))
      })
      return(c(out[1,], out[2,]))
    }) %>% t()
    browser()
    if (length(pairs)>2) {
      p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
                    method="BH") %>% matrix(ncol=length(groups), byrow=F)
      diff.out <- cbind(diff.out, p)
      colnames(diff.out) <- c(header$c, header$p, header$f)
    }
    else {colnames(diff.out) <- c(header$c, header$p)}
  } else if (method=="kruskal") {
    header <- format_differential_exposures(c(''), c("(K-W chi-squared)", "(df)", "(p-value)"))
      # data.frame(y=c('')) %>%
      # dplyr::mutate(c=paste(.data$y,"(K-W chi-squared)", sep=""),
      #        df=paste(.data$y,"(df)", sep=""),
      #        p=paste(.data$y,"(p-value)", sep=""))
    
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- kruskal.test(y ~ annotations, ...)
      return (c(out$statistic, out$parameter, out$p.value))
    }) %>% t()
    colnames(diff.out) <- c(header[-1])#header$c, header$df, header$p)
  } else if (method=="glm.nb") {
    header <- data.frame(y=groups) %>%
      dplyr::mutate(coef=paste(.data$y,"(coef)", sep=""),
             sd=paste(.data$y,"(Std. Error)", sep=""),
             z=paste(.data$y,"(z)", sep=""),
             p=paste(.data$y, "(Pr(>|z|))", sep=""),
             adj=paste(.data$y, "(fdr)", sep=""))
    diff.out <- apply(exposures, 1, FUN=function(y) {
      #Add anova p and fdr values
      out <- summary(MASS::glm.nb(round(y) ~ annotations, ...))$coefficients
      return(out)
    }) %>% t()
    p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
                  method="BH") %>% matrix(ncol=length(groups), byrow=F)
    diff.out <- cbind(diff.out, p)
    colnames(diff.out) <- c(header$coef, header$sd, header$z, header$p, header$adj)
  } else {
    stop("Method is not supported. Please provide one of: 
         wilcox, kruskal, glm.nb")
  }
  
  return (diff.out)
}

#' Formatting output of \code{differential_exposures} function
#' @param prefix A list of strings which will prefix the rest of the output columns.
#' @param columns A list of strings corresponding to the output of the 
#' \code{differential_exposures} function
#' @return A dataframe of strings which will be used in 
#' \code{differential_exposures} to format the header
#' @example 
#' x = c("a", "b", "c")
#' y = c("x", "y", "z")
#' format_output(paste(x,y,sep"-"), c("(coef)", "(p-value)"))
#' @keywords internal
format_differential_exposures <- function(prefix, columns) {
  header <- data.frame(c(prefix))
  for (column in columns){
    header[[column]] = paste0(prefix,column)
  }
  return (header)
}

TESTDATA <- "/Users/nathansahelijo/Documents/BostonUniversity/class/challenge/scripts/data/"
