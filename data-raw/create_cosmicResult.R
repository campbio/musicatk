library(musicatk)
cosmic=read.table("cosmic_signatures.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
cosmic_mat=cosmic[,4:ncol(cosmic)]
rownames(cosmic_mat)=paste(cosmic$Substitution.Type, cosmic$Trinucleotide, sep="_")
colnames(cosmic_mat)=paste("Signature", seq_len(30), sep="")
cosmic_v2_sigs <- new("musica_result", signatures = as.matrix(cosmic_mat), type = "Cosmic")
usethis::use_data(cosmic_v2_sigs, internal = FALSE, overwrite = TRUE)

v3_snv_tab_exome <- read.table("sigProfiler_exome_SBS_signatures.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
exome_mat=v3_snv_tab_exome[,3:ncol(v3_snv_tab_exome)]
rownames(exome_mat)=paste(v3_snv_tab_exome$Type, v3_snv_tab_exome$SubType, sep="_")
cosmic_v3_sbs_sigs_exome <- new("musica_result", signatures = as.matrix(exome_mat), type = "Cosmic_v3_SNV_Exome")
usethis::use_data(cosmic_v3_sbs_sigs_exome, internal = FALSE, overwrite = TRUE)

v3_snv_tab <- read.table("sigProfiler_SBS_signatures_2019_05_22.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
snv_mat=v3_snv_tab[,3:ncol(v3_snv_tab)]
rownames(snv_mat)=paste(v3_snv_tab$Type, v3_snv_tab$SubType, sep="_")
cosmic_v3_sbs_sigs <- new("musica_result", signatures = as.matrix(snv_mat), type = "Cosmic_v3_SNV_Genome")
usethis::use_data(cosmic_v3_sbs_sigs, internal = FALSE, overwrite = TRUE)

v3_dbs_tab <- read.table("sigProfiler_DBS_signatures.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
row.names(v3_dbs_tab) = gsub('>', '>NN_', row.names(v3_dbs_tab))
cosmic_v3_dbs_sigs <- new("musica_result", signatures = as.matrix(v3_dbs_tab), type = "Cosmic_v3_DBS_Genome")
usethis::use_data(cosmic_v3_dbs_sigs, internal = FALSE, overwrite = TRUE)

v3_indel_tab <- read.table("sigProfiler_ID_signatures.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)
cosmic_v3_indel_sigs <- new("musica_result", signatures = as.matrix(v3_indel_tab), type = "Cosmic_v3_Indel_Genome")
usethis::use_data(cosmic_v3_indel_sigs, internal = FALSE, overwrite = TRUE)

usethis::use_data(cosmic_v2_sigs, cosmic_v3_sbs_sigs_exome, cosmic_v3_sbs_sigs, cosmic_v3_dbs_sigs, cosmic_v3_indel_sigs, internal = TRUE, overwrite = TRUE)
