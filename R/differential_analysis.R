#' Compare exposures from annotated samples
#' 
#' @importFrom magrittr %>%
#' @importFrom MASS glm
#' @importFrom MASS glm.nb
#' @param musica_result A \code{\linkS4class{musica_result}} object 
#' @param annotation Column in the sample_annotations table of the 
#' \code{\linkS4class{musica_result}} object
#' @param method Any method in c("wilcox", "kruskal", "glm.nb") used to perform differential analysis on signature 
#' exposures
#' @param ... Additional arguments to be passed to the chosen method
#' @return A matrix containing statistics summarizing the analysis dependent
#' on the chosen method
#' @examples 
#' musica <- readRDS(system.file("testdata", "res_annot.rds", package = "musicatk"))
#' compare_samples(musica, "Tumor_Subtypes", method="wilcox")
#' @export
#' 
compare_samples <- function(musica_result, annotation, method="wilcox",...) {
  annotations <- factor(
    musica_result@musica@sample_annotations[[annotation]])
  annotations <- factor(
    annotations[match(musica_result@musica@sample_annotations$Samples,
                      colnames(musica_result@exposures))])
  diff.out <- 0
  exposures <- musica_result@exposures
  l <- length(exposures)
  groups <- unique(annotations)
  if (method=="wilcox" || is.null(func)) {
    annotations <- as.integer(annotations)
    pairs <- combn(groups,2) %>% t()
    header <- data.frame(y=pairs[,1], x=pairs[,2]) %>%
      mutate(c=paste(y,"-",x,"(W)", sep=""),
             p=paste(y,"-",x,"(p-value)", sep=""),
             f=paste(y,"-",x,"(fdr)", sep=""))
    
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- apply(pairs, 1, FUN=function(p) {
        out <- wilcox.test(y[annotations==p[1]],y[annotations==p[2]],...)
        
        return (c(s=out$statistic, p=out$p.value))
      })
      return(c(out[1,], out[2,]))
    }) %>% t()
    
    if (length(pairs)>2) {
      p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
                    method="BH") %>% matrix(ncol=length(groups), byrow=F)
      diff.out <- cbind(diff.out, p)
      colnames(diff.out) <- c(header$c, header$p, header$f)
    }
    else {colnames(diff.out) <- c(header$c, header$p)}
  } else if (method=="kruskal") {
    header <- data.frame(y=c('')) %>%
      mutate(c=paste(y,"(K-W chi-squared)", sep=""),
             df=paste(y,"(df)", sep=""),
             p=paste(y,"(p-value)", sep=""))
    
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- kruskal.test(y ~ annotations, ...)
      return (c(out$statistic, out$parameter, out$p.value))
    }) %>% t()
    colnames(diff.out) <- c(header$c, header$df, header$p)
    
  } else if (method=="glm.nb") {
    header <- data.frame(y=groups) %>%
      mutate(coef=paste(y,"(coef)", sep=""),
             sd=paste(y,"(Std. Error)", sep=""),
             z=paste(y,"(z)", sep=""),
             p=paste(y, "(Pr(>|z|))", sep=""),
             adj=paste(y, "(p.adj)", sep=""))
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- summary(MASS::glm.nb(round(y) ~ annotations))$coefficients 
    }) %>% t()
    p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
                  method="BH") %>% matrix(ncol=length(groups), byrow=F)
    diff.out <- cbind(diff.out, p)
    colnames(diff.out) <- c(header$coef, header$sd, header$z, header$p, header$adj)
  } else {
    stop("Invalid method given.")
  }
  
  return (diff.out)
}
compare_samples(res_annot, "Tumor_Subtypes", method="glm.nb")
