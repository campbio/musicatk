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



compare_samples <- function(musica_result, annotation,method="wilcox",...) {
  annotations <- factor(
    musica_result@musica@sample_annotations[[annotation]])
  annotations <- 
    annotations[match(musica_result@musica@sample_annotations$Samples,
                      colnames(musica_result@exposures))]
  diff.out <- 0

  if (method=="wilcox" || is.null(func)) {
    diff.out <- apply(musica_result@exposures, 1, FUN=function(y) {
      wilcox.test(y~annotations)
    })
    #   apply(
    # func <- function(y){
    #   return(sapply(annotations, FUN=function(a){
    #     browser()
    #   return()
    # }))}
  } else if (method=="kruskal") {
    func <- kruskal.test
  } else if (method=="glm") {
    # func <- glm.nb
    # table_header <- unique(annotations) %>% data.frame() %>%
    #   mutate(coef=paste(.,"(coef)",sep=""),
    #          stat=paste(.,"(z)", sep=""),
    #          p=paste(.,"(Pr(>|z|))", sep="")) %>% as.vector()
    # #browser()
    # diff.out <- apply(musica_result@exposures, 1, FUN=function(y) {
    #   out <- summary(func(y~annotations))$coefficients %>% as.numeric() %>%
    #     as.data.frame() %>% t()
    #   #browser()
    #   colnames(out) <- table_header
    #   print(out)
    #   return (out)
    #   })
  } else {
    stop("Invalid function given")
  }
  # diff.out <- apply(musica_result@exposures, 1, FUN=function(y) {
  #   out <- func(y)
    #out <- c(w.stat=out$statistic[[1]], p.value=out[["p.value"]])
    #out <- summary(glm.nb(y~annotations,...))$coefficients %>% data.frame()
    #c(glm.score=out[2,"z value"], glm.pvalue=out[2,"Pr(>|z|)"])
    # LM
    # out <- summary(lm(y~annotations))$coefficients
    # c(lm.score=out[2,"t value"], lm.pvalue=out[2,"Pr(>|t|)"])
    #out[[4]] <- p.adjust(out[[4]], method="BH")
    #browser()
    #return (out)
  #})
  #browser()
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




