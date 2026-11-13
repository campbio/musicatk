#' musica
#'
#' Variants of 12 TCGA tumors from the example files in
#' \code{inst/extdata}, with an SBS96 table: four lung tumors with a strong
#' smoking signature (COSMIC SBS4), four lung tumors with a strong APOBEC
#' signature (SBS2 and SBS13) and four melanomas with a strong UV signature
#' (SBS7a and SBS7b).
#'
#' @docType data
#'
#' @usage data(musica)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{create_musica_from_variants}}.
#'
#' @source TCGA, hg38: the example files in \code{inst/extdata}
#' (\code{public_LUSC_TCGA-22-5472-01A-01D-1632-08.vcf},
#' \code{public_TCGA_lung.maf} and four \code{public_SKCM_*.vcf} files),
#' from GDC open-access aliquot ensemble masked MAFs. The samples were chosen
#' by \code{data-raw/choose_example_samples.R}, the files written by
#' \code{data-raw/example_files.R}, and the objects built by
#' \code{data-raw/example_musica.R}.
#'
#' @keywords datasets
#'
"musica"

#' musica_annot
#'
#' \code{\link{musica}} with the \code{Tumor_Subtypes} sample annotation
#' (\code{"Lung"} or \code{"Skin"}) from
#' \code{inst/extdata/sample_annotations.txt}.
#'
#' @docType data
#'
#' @usage data(musica_annot)
#'
#' @format An object of class \code{musica}. See \code{\link{samp_annot}}.
#'
#' @inherit musica source
#'
#' @keywords datasets
#'
"musica_annot"

#' musica_sbs96
#'
#' The single base substitutions (SBS) of \code{\link{musica}}, with an
#' SBS96 table.
#'
#' @docType data
#'
#' @usage data(musica_sbs96)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{build_standard_table}}.
#'
#' @inherit musica source
#'
#' @keywords datasets
#'
"musica_sbs96"

#' musica_sbs96_tiny
#'
#' A very small musica for testing: the first 40 SBS by position of two
#' lung tumors in \code{\link{musica}}, one with a smoking and one with an
#' APOBEC signature, with an SBS96 table.
#'
#' @docType data
#'
#' @usage data(musica_sbs96_tiny)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{build_standard_table}}.
#'
#' @inherit musica source
#'
#' @keywords datasets
#'
"musica_sbs96_tiny"

#' dbs_musica
#'
#' All variants of five TCGA lung adenocarcinomas, the open-access GDC
#' TCGA-LUAD MAFs with the most doublet base substitutions (DBS), with
#' SBS96 and DBS78 tables.
#'
#' @docType data
#'
#' @usage data(dbs_musica)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{create_musica_from_variants}} and
#' \code{\link{build_standard_table}}.
#'
#' @source GDC open-access "Aliquot Ensemble Somatic Variant Merging and
#' Masking" MAFs (TCGA, WXS, hg38), GDC Data Release 46.0. Built by
#' \code{data-raw/gdc_example_musica.R}, which lists the GDC file IDs.
#'
#' @keywords datasets
#'
"dbs_musica"

#' indel_musica
#'
#' The insertions and deletions of five TCGA uterine corpus endometrial
#' carcinomas, with no count tables.
#'
#' @docType data
#'
#' @usage data(indel_musica)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{create_musica_from_variants}} and
#' \code{\link{build_standard_table}}.
#'
#' @inherit dbs_musica source
#'
#' @keywords datasets
#'
"indel_musica"

#' res
#'
#' \code{\link{musica}} with a three-signature NMF model on SBS96
#' (\code{model_id = "res"}). Signature1 is an APOBEC signature (closest to
#' COSMIC SBS2 and SBS13) and Signature2 a smoking signature (SBS4), found in
#' the lung tumors; Signature3 is a UV signature (SBS7a), found in the
#' melanomas.
#'
#' @docType data
#'
#' @usage data(res)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{discover_signatures}}.
#'
#' @inherit musica source
#'
#' @keywords datasets
#'
"res"

#' res_annot
#'
#' \code{\link{musica_annot}} with a three-signature NMF model on SBS96
#' (\code{model_id = "res_annot"}), as in \code{\link{res}}.
#'
#' @docType data
#'
#' @usage data(res_annot)
#'
#' @format An object of class \code{musica}. See
#' \code{\link{discover_signatures}}.
#'
#' @inherit musica source
#'
#' @keywords datasets
#'
"res_annot"

#' synthetic_breast_counts
#'
#' A data.frame containing the SBS96 mutation counts for a synthetic breast
#' cancer dataset with 214 samples.
#'
#' @docType data
#'
#' @usage data(synthetic_breast_counts)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
"synthetic_breast_counts"

#' synthetic_breast_true_exposures
#'
#' A data.frame containing true signature exposure levels for a synthetic
#' breast cancer dataset with 214 samples and 8 signatures.
#'
#' @docType data
#'
#' @usage data(synthetic_breast_true_exposures)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
"synthetic_breast_true_exposures"

#' example_predicted_sigs
#'
#' A matrix containing the predicted signatures for a synthetic breast cancer
#' dataset with 214 samples, generated using NMF and k=8 signatures.
#'
#' @docType data
#'
#' @usage data(example_predicted_sigs)
#'
#' @format An object of class \code{matrix}
#'
#' @keywords datasets
#'
"example_predicted_sigs"

#' example_predicted_exp
#'
#' A matrix containing the predicted signature exposure levels for a synthetic
#' breast cancer dataset with 214 samples, generated using NMF and k=8
#' signatures.
#'
#' @docType data
#'
#' @usage data(example_predicted_exp)
#'
#' @format An object of class \code{matrix}
#'
#' @keywords datasets
#'
"example_predicted_exp"

#' full_benchmark_example
#'
#' A \code{full_benchmark} object containing ground truth signatures and
#' exposures
#' from a synthetic breast cancer dataset with 8 known signatures as well as
#' a disocvery prediction that has been benchmarked.
#'
#' @docType data
#'
#' @usage data(full_benchmark_example)
#'
#' @format An object of class \code{full_benchmark}
#'
#' @keywords datasets
#'
"full_benchmark_example"
