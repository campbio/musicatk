library(TCGAbiolinks)
meso_maf <- GDCquery_Maf("MESO", pipelines = "mutect")

#data.frame("Missense Variant" = 
#             str_detect(meso_maf$all_effects, "missense_variant")) 
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
meso.result <- discover_signatures(musica = meso.musica, table_name = "SBS96", 
                                   num_signatures = 10, method = "lda", 
                                   nstart = 10)
meso.samples <- get_sample_names(meso.musica)
set_sample_annotations(meso.musica, data.table::data.table(meso_maf$BIOTYPE))

#https://github.com/compbiomed/singleCellTK
https://satijalab.org/seurat/
https://github.com/satijalab/seurat

https://www.sctk.science/index.html