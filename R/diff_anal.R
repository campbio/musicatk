library(TCGAbiolinks)
library(MASS)
meso_maf <- GDCquery_Maf("MESO", pipelines = "mutect")
thca_maf <- GDCquery_Maf("THCA", pipelines = "mutect")

g <- select_genome("hg38")
meso.variants <- extract_variants_from_matrix(meso_maf, 
                                              chromosome_col="Chromosome", 
                                              end_col="End_Position", 
                                              start_col="Start_Position", 
                                              ref_col="Tumor_Seq_Allele1", 
                                              alt_col="Tumor_Seq_Allele2", 
                                              sample_col="Tumor_Sample_Barcode")
meso.musica <- create_musica(x = meso.variants, genome = g)
build_standard_table(meso.musica, g = g, table_name = "SBS96")
meso.result <- discover_signatures(musica = meso.musica, 
                                   table_name = "SBS96", 
                                   num_signatures = 5, 
                                   method = "lda", 
                                   nstart = 10)
meso.samples <- get_sample_names(meso.musica)
init_sample_annotations(meso.musica)
thca.variants <- extract_variants_from_matrix(thca_maf, 
                                              chromosome_col="Chromosome", 
                                              end_col="End_Position", 
                                              start_col="Start_Position", 
                                              ref_col="Tumor_Seq_Allele1", 
                                              alt_col="Tumor_Seq_Allele2", 
                                              sample_col="Tumor_Sample_Barcode")

# add_sample_annotations(meso.musica, meso_maf, 
#                        sample_column="Tumor_Sample_Barcode",
#                        columns_to_add="BIOTYPE")
# set_sample_annotations(meso.musica, 
#                        data.table::data.table(meso_maf)$BIOTYPE)
# protein_coding.musica <- subset_musica_by_annotation(meso.musica, 
#                                                      "V1", 
#                                                      "protein_coding")
# lincRNA.musica <- subset_musica_by_annotation(meso.musica, 
#                                               "V1", 
#                                               "lincRNA")
# Compare signatures in musica by annotation
# pc.result <- discover_signatures(musica = protein_coding.musica, 
#                                  table_name = "SBS96", 
#                                  num_signatures = 5, 
#                                  method = "lda", 
#                                  nstart = 5)
# lincRNA.result <- discover_signatures(musica = lincRNA.musica, 
#                                       table_name = "SBS96", 
#                                       num_signatures = 5, 
#                                       method = "lda", 
#                                       nstart = 5)

thca.result <- discover_signatures(musica = thca.musica,
                                 table_name = "SBS96",
                                 num_signatures = 4,
                                 method = "lda")
meso.result <- discover_signatures(musica = meso.musica,
                                      table_name = "SBS96",
                                      num_signatures = 5,
                                      method = "lda")
# sig <- data.frame(x=meso.result@signatures[,1], y=thca.result@signatures[,1])
# lm(y~x, data=sig)
glm.nb(x~y, data=sig)
for (x in colnames(meso.result@signatures)) {
  exposures <- cbind(x=meso.result@exposures[x,], 
                     y=thca.result@exposures[x,]) %>% data.frame()
  #exposures <- Seurat::NormalizeData(exposures) %>% Seurat::ScaleData()
  #print(lm(y~x, data=data.frame(exposures))$coefficients)
  print(MASS::glm.nb(y~x, data=exposures)$coefficients)
}


compare_samples(meso.result, thca.result)
#Performs differential analysis on 2 sets of signatures
compare_samples <- function(a.result, b.result) {
  # sapply(colnames(a.result@signatures), FUN = function(x) {
  #   #print(data.frame(a.result@signatures) )
  #   length(a.result@signatures[,x]) == length(b.result@signatures[,x]) 
  #   #out <- glm.nb(x~b.result@signatures[,x], data.frame(a.result@signatures) )%>% summary()$coefficients
  #   #other <- b.result
  #   #glm.nb()
  # })
  for (x in colnames(a.result@signatures)) {
    exposures <- cbind(x=a.result@exposures[x,], 
                       y=b.result@exposures[x,]) %>% data.frame()
    #exposures <- Seurat::NormalizeData(exposures) %>% Seurat::ScaleData()
    #print(lm(y~x, data=data.frame(exposures))$coefficients)
    coef <- summary(MASS::glm.nb(y~x, data=exposures))$coefficients
    print(coef)
  }
}





# https://github.com/compbiomed/singleCellTK
# https://satijalab.org/seurat/
# https://github.com/satijalab/seurat
# 
# https://www.sctk.science/index.html