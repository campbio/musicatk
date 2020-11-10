#' @importFrom magrittr %>%
#' @importFrom MASS glm
#' @importFrom MASS glm.nb
#' 
compare_samples <- function(musica_result, annotation,...) {
  annotations <- factor(
    musica_result@musica@sample_annotations[[annotation]])
  annotations <- 
    annotations[match(musica_result@musica@sample_annotations$Samples,
                      colnames(musica_result@exposures))]
  diff.out <- apply(musica_result@exposures, 1, FUN=function(y) {
    out <- summary(glm.nb(y~annotations,...))$coefficients %>% data.frame()
    #c(glm.score=out[2,"z value"], glm.pvalue=out[2,"Pr(>|z|)"])
    # LM
    # out <- summary(lm(y~annotations))$coefficients
    # c(lm.score=out[2,"t value"], lm.pvalue=out[2,"Pr(>|t|)"])
    #out[[4]] <- p.adjust(out[[4]], method="BH")
    return (out)
  })
  #browser()
  return (diff.out)
}
