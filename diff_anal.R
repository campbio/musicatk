library(TCGAbiolinks)
library(MASS)
library(dplyr)
g <- select_genome("hg38")

meso_maf <- GDCquery_Maf("MESO", pipelines = "mutect") %>% 
  mutate(.keep="all", Tumor_Type="meso")
thca_maf <- GDCquery_Maf("THCA", pipelines = "mutect") %>% 
  mutate(.keep="all", Tumor_Type="thca")
dmy1 <- meso_maf %>% mutate(Tumor_Type="A")
dmy2 <- meso_maf %>% mutate(Tumor_Type="B",
                            Tumor_Sample_Barcode=paste("DummyB-",
                                                       Tumor_Sample_Barcode, 
                                                       sep=""))
dmy_maf <- rbind(dmy1,dmy2)
mix_maf <- rbind(meso_maf, thca_maf, dmy2)

meso.variants <- extract_variants_from_matrix(meso_maf,
                                             chromosome_col="Chromosome", 
                                             end_col="End_Position", 
                                             start_col="Start_Position", 
                                             ref_col="Tumor_Seq_Allele1", 
                                             alt_col="Tumor_Seq_Allele2", 
                                             sample_col="Tumor_Sample_Barcode"
)
meso.musica <- create_musica(x=meso.variants, genome=g)
build_standard_table(meso.musica, g = g, table_name = "SBS96")
meso.result <- discover_signatures(musica = meso.musica, 
                                  table_name = "SBS96", 
                                  num_signatures = 7, 
                                  method = "lda")
init_sample_annotations(meso.musica)
add_sample_annotations(meso.musica, data.table::data.table(meso_maf),
                       sample_column="Tumor_Sample_Barcode",
                       columns_to_add="Tumor_Type")

meso.result <- discover_signatures(musica = meso.musica,
                                  table_name = "SBS96",
                                  num_signatures = 5,
                                  method = "lda")

# THCA
thca.variants <- extract_variants_from_matrix(thca_maf,
                                             chromosome_col="Chromosome", 
                                             end_col="End_Position", 
                                             start_col="Start_Position", 
                                             ref_col="Tumor_Seq_Allele1", 
                                             alt_col="Tumor_Seq_Allele2", 
                                             sample_col="Tumor_Sample_Barcode"
)
thca.musica <- create_musica(x=thca.variants, genome=g)
build_standard_table(thca.musica, g = g, table_name = "SBS96")
thca.result <- discover_signatures(musica = thca.musica, 
                                  table_name = "SBS96", 
                                  num_signatures = 7, 
                                  method = "lda")
init_sample_annotations(thca.musica)
add_sample_annotations(thca.musica, data.table::data.table(thca_maf),
                       sample_column="Tumor_Sample_Barcode",
                       columns_to_add="Tumor_Type")

thca.result <- discover_signatures(musica = thca.musica,
                                  table_name = "SBS96",
                                  num_signatures = 5,
                                  method = "lda")

# MIXED
mix.variants <- extract_variants_from_matrix(mix_maf,
                                    chromosome_col="Chromosome", 
                                    end_col="End_Position", 
                                    start_col="Start_Position", 
                                    ref_col="Tumor_Seq_Allele1", 
                                    alt_col="Tumor_Seq_Allele2", 
                                    sample_col="Tumor_Sample_Barcode"
                                    )
mix.musica <- create_musica(x=mix.variants, genome=g)
build_standard_table(mix.musica, g = g, table_name = "SBS96")
mix.result <- discover_signatures(musica = mix.musica, 
                                  table_name = "SBS96", 
                                  num_signatures = 7, 
                                  method = "lda")
init_sample_annotations(mix.musica)
add_sample_annotations(mix.musica, data.table::data.table(mix_maf),
                       sample_column="Tumor_Sample_Barcode",
                       columns_to_add="Tumor_Type")

mix.result <- discover_signatures(musica = mix.musica,
                                  table_name = "SBS96",
                                  num_signatures = 5,
                                  method = "lda")

# DUMMY
dmy.variants <- extract_variants_from_matrix(dmy_maf,
                                             chromosome_col="Chromosome", 
                                             end_col="End_Position", 
                                             start_col="Start_Position", 
                                             ref_col="Tumor_Seq_Allele1", 
                                             alt_col="Tumor_Seq_Allele2", 
                                             sample_col="Tumor_Sample_Barcode"
)
dmy.musica <- create_musica(x=dmy.variants, genome=g)
build_standard_table(dmy.musica, g = g, table_name = "SBS96")
mix.result <- discover_signatures(musica = dmy.musica, 
                                  table_name = "SBS96", 
                                  num_signatures = 7, 
                                  method = "lda")
init_sample_annotations(dmy.musica)
add_sample_annotations(dmy.musica, data.table::data.table(dmy_maf),
                       sample_column="Tumor_Sample_Barcode",
                       columns_to_add="dummy")

dmy.result <- discover_signatures(musica = dmy.musica,
                                  table_name = "SBS96",
                                  num_signatures = 5,
                                  method = "lda")

# MESO
meso.variants <- extract_variants_from_matrix(meso_maf,
                                             chromosome_col="Chromosome", 
                                             end_col="End_Position", 
                                             start_col="Start_Position", 
                                             ref_col="Tumor_Seq_Allele1", 
                                             alt_col="Tumor_Seq_Allele2", 
                                             sample_col="Tumor_Sample_Barcode"
)
meso.musica <- create_musica(x=meso.variants, genome=g)
build_standard_table(meso.musica, g = g, table_name = "SBS96")
meso.result <- discover_signatures(musica = meso.musica, 
                                  table_name = "SBS96", 
                                  num_signatures = 7, 
                                  method = "lda")

init_sample_annotations(meso.musica)

#Must run init_sample_anonotations and use data.table. 
#DO NOT USE set_sample_annotations
add_sample_annotations(meso.musica, data.table::data.table(meso_maf),
                       sample_column="Tumor_Sample_Barcode",
                       columns_to_add="Tumor_Type")
meso.result <- discover_signatures(musica = meso.musica,
                                  table_name = "SBS96",
                                  num_signatures = 5,
                                  method = "lda")



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
    p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
             method="BH") %>% matrix(ncol=length(groups), byrow=F)
    diff.out <- cbind(diff.out, p)
    colnames(diff.out) <- c(header$c, header$p, header$f)
  } else if (method=="kruskal") {
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- kruskal.test(y ~ annotations, ...)
      browser()
    })
  } else if (method=="glm") {
    header <- list(coef=c(), sd=c(), z=c(), p=c(), adj=c())
    for (a in groups) {
      header$coef <- append(header$coef,paste(a,"(coef)", sep=""))
      header$sd <- append(header$sd, paste(a,"Std. Error", sep=""))
      header$z <- append(header$z, paste(a, "(z)", sep=""))
      header$p <- append(header$p, paste(a, "(Pr(>|z|))", sep=""))
      header$adj <- append(header$adj, paste(a, "(p.adj)", sep=""))
    }
    header <- unlist(header)
    diff.out <- apply(exposures, 1, FUN=function(y) {
      out <- summary(MASS::glm.nb(round(y) ~ annotations))$coefficients 
      }) %>% t()
    p <- p.adjust(diff.out[,(ncol(diff.out)-length(groups)+1):ncol(diff.out)], 
             method="BH") %>% matrix(ncol=length(groups), byrow=F)
    diff.out <- cbind(diff.out, p)
    colnames(diff.out) <- header
  } else {
    stop("Invalid method given.")
  }
  return (diff.out)
}

compare_samples(mix.result, "Tumor_Type", method="wilcox")

plot_exposures(meso.result, proportional = F)
plot_exposures(thca.result, proportional = F)


compare_samples(dmy.result, "dummy")
plot(dmy.result@exposures["Signature1",
                          dmy.result@musica@sample_annotations$dummy=="A"],
     dmy.result@exposures["Signature1",
                          dmy.result@musica@sample_annotations$dummy=="B"])



B <- dmy.musica@sample_annotations %>% filter(dummy=="B") 
B <- B$Samples
testB <- dmy.result@exposures[, colnames(dmy.result@exposures) %in% B] * 2
test.result <- dmy.result
test.result@exposures <- cbind(
  dmy.result@exposures[,!(colnames(dmy.result@exposures) %in% B)], 
  testB)
compare_samples(test.result, "dummy")




