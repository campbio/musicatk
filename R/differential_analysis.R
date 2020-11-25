#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom MASS glm.nb

#' @title Compare exposures of annotated samples
#' @description \code{compare_samples} is used to run differential analysis on 
#' the signature exposures of annotated samples within the 
#' \code{\linkS4class{musica_result}} object.
#' @param musica_result A \code{\linkS4class{musica_result}} object 
#' @param annotation Column in the sample_annotations table of the 
#' \code{\linkS4class{musica_result}} object
#' @param method Any method in \code{c("wilcox", "kruskal", "glm.nb")} 
#' used to perform differential analysis on signature exposures
#' @param ... Additional arguments to be passed to the chosen method
#' @return A matrix containing statistics summarizing the analysis dependent
#' on the chosen method
#' @examples 
#' data("res_annot")
#' compare_samples(res_annot, "Tumor_Subtypes", method="wilcox")
#' @export
compare_samples <- function(musica_result, annotation, method="wilcox",...) {
  if (!methods::is(musica_result, "musica_result")) {
    stop("Input to compare_samples must be a 'musica_result' object.")
  }
  annotations <- 
    musica_result@musica@sample_annotations[[annotation]]
  if (is.null(annotations)) {
    stop(paste('"',annotation,'" does not exist in musica_result.', sep=""))
  }
  exposures <- exposures(musica_result)
  groups <- unique(annotations)
  annotations <- 
    annotations[match(musica_result@musica@sample_annotations$Samples,
                      colnames(exposures))]
  annotations <- factor(annotations, levels=unique(groups))
  diff.out <- 0
  l <- length(exposures)
  browser()
  if (method=="wilcox" || is.null(method)) {
    pairs <- combn(groups,2) %>% t()
    header <- data.frame(y=pairs[,1], x=pairs[,2]) %>%
      dplyr::mutate(c=paste(.data$y,"-",.data$x,"(W)", sep=""),
             p=paste(.data$y,"-",.data$x,"(p-value)", sep=""),
             f=paste(.data$y,"-",.data$x,"(fdr)", sep=""))
    
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
      dplyr::mutate(c=paste(.data$y,"(K-W chi-squared)", sep=""),
             df=paste(.data$y,"(df)", sep=""),
             p=paste(.data$y,"(p-value)", sep=""))
    
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- kruskal.test(y ~ annotations, ...)
      return (c(out$statistic, out$parameter, out$p.value))
    }) %>% t()
    colnames(diff.out) <- c(header$c, header$df, header$p)
    
  } else if (method=="glm.nb") {
    header <- data.frame(y=groups) %>%
      dplyr::mutate(coef=paste(.data$y,"(coef)", sep=""),
             sd=paste(.data$y,"(Std. Error)", sep=""),
             z=paste(.data$y,"(z)", sep=""),
             p=paste(.data$y, "(Pr(>|z|))", sep=""),
             adj=paste(.data$y, "(p.adj)", sep=""))
    diff.out <- apply(exposures, 1, FUN=function(y) {
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

# lusc_maf <- system.file("extdata", "public_TCGA.LUSC.maf", package = "musicatk")
# lusc.variants <- extract_variants_from_maf_file(maf_file=lusc_maf) %>% 
#   mutate(.keep="all", Tumor_Type="meso")
# luad_vcf <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf", package = "musicatk") 
# luad.variants <- extract_variants_from_vcf_file(vcf_file = luad_vcf) %>% 
#   mutate(.keep="all", Tumor_Type="thca")
# 
# dmy1 <- lusc.variants %>% mutate(Tumor_Type="lusc_A",
#                             sample=paste("lusc_A-", sample, sep=""))
# dmy2 <- lusc.variants %>% mutate(Tumor_Type="lusc_B",
#                             sample=paste("lusc_B-", sample, sep=""))
# 
# mix.variants <- rbind(dmy1, luad.variants, dmy2)
# mix.musica <- create_musica(x=mix.variants, genome=g)
# 
# build_standard_table(mix.musica, g = g, table_name = "SBS96")
# samp_annot(mix.musica)
# samp_annot(mix.musica, "Tumor_Type") <- 
#   c("lusc_A", "lusc_A", "lusc_A", "lusc_B", "lusc_B", "lusc_B","luad")
# 
# res_diff <- discover_signatures(musica = mix.musica,
#                                   table_name = "SBS96",
#                                   num_signatures = 7,
#                                   method = "lda")
# compare_samples(res_diff, "Tumor_Type")
