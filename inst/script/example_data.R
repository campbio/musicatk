# The example files in inst/extdata. The 12 TCGA samples (hg38) were
# chosen by data-raw/choose_example_samples.R from the GDC open-access
# "Aliquot Ensemble Somatic Variant Merging and Masking" MAFs of TCGA-LUAD,
# TCGA-LUSC and TCGA-SKCM: four lung tumors each with a strong smoking
# (COSMIC SBS4) or APOBEC (SBS2 and SBS13) signature, and four melanomas
# with a strong UV (SBS7a and SBS7b) signature. data-raw/example_files.R
# downloads their MAFs by GDC file ID and writes the files below.

#' Lung squamous cell carcinoma VCF from TCGA
#'
#' @name raw_public_LUSC_TCGA-22-5472-01A-01D-1632-08
#'
#' @section public_LUSC_TCGA-22-5472-01A-01D-1632-08.vcf
#'
#' A lung tumor with a strong smoking signature, written as a VCF by
#' data-raw/example_files.R. This data is used in the vignette and in
#' several examples
NULL

#' Melanoma VCFs from TCGA
#'
#' @name raw_public_SKCM
#'
#' @section public_SKCM_TCGA-BF-A1PX-01A-12D-A19A-08.vcf,
#' public_SKCM_TCGA-EB-A3XD-01A-22D-A23B-08.vcf,
#' public_SKCM_TCGA-EB-A41B-01A-11D-A24R-08.vcf,
#' public_SKCM_TCGA-ER-A3PL-06A-11D-A23B-08.vcf
#'
#' Four melanomas with a strong UV signature, one VCF each, written by
#' data-raw/example_files.R. This data is used in the vignette and in
#' several examples
NULL

#' Lung tumor MAF from TCGA
#'
#' @name raw_public_TCGA_lung
#'
#' @section public_TCGA_lung.maf
#'
#' Seven lung tumors, three with a strong smoking and four with a strong
#' APOBEC signature, with the standard TCGA MAF columns and read counts,
#' written by data-raw/example_files.R. This data is used in the vignette
#' and in several examples
NULL

#' Example sample annotations that maps sample names to tumor subtypes
#'
#' @name raw_sample_annotations
#'
#' @section sample_annotations.txt
#'
#' The tumor type (Lung or Skin) of each of the 12 example samples, written
#' by data-raw/example_files.R. This data is used in the vignette and in
#' several examples
NULL

#' Example count table with LUAD data from TCGA
#'
#' @name luad_tcga_count_table
#'
#' @section luad_tcga_count_table.csv
#'
#' This data is used in the vignette
NULL
