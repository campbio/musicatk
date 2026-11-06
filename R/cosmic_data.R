#' COSMIC v2 SBS96 Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v2_sigs)
#'
#' @format An object of class \code{result_model}
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L., Nik-Zainal, S., Wedge, D. et al. (2013)
#' Signatures of mutational processes in human cancer. Nature 500, 415–421
#' ([Nature](https://www.ncbi.nlm.nih.gov/pubmed/23945592))
#'
#' @source COSMIC v2, <https://cancer.sanger.ac.uk/cosmic/signatures_v2>
"cosmic_v2_sigs"

#' COSMIC v3 SBS96 Exome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_sbs_sigs_exome)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94–101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3, <https://cancer.sanger.ac.uk/cosmic/signatures>
"cosmic_v3_sbs_sigs_exome"

#' COSMIC v3 SBS96 Genome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_sbs_sigs)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94–101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3, <https://cancer.sanger.ac.uk/cosmic/signatures>
"cosmic_v3_sbs_sigs"

#' COSMIC v3 DBS Genome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_dbs_sigs)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94–101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3, <https://cancer.sanger.ac.uk/cosmic/signatures>
"cosmic_v3_dbs_sigs"

#' COSMIC v3 Indel Genome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_indel_sigs)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94–101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3, <https://cancer.sanger.ac.uk/cosmic/signatures>
"cosmic_v3_indel_sigs"

#' Replication strand of the hg19 genome as a GRanges object
#'
#' Leading and lagging replication strand in 10 kb bins on the hg19
#' (GRCh37) genome build, from the replication fork directionality (RFD)
#' that GLOE-seq measured in HCT116 colon cancer cells with ligase 1
#' depleted. Use it only with variants called on hg19; for hg38 use
#' [rep_range_hg38].
#'
#' Each bin is labelled by the sign of its RFD, (REV - FWD) / (REV + FWD).
#' In these cells the breaks GLOE-seq captures lie on the newly made
#' lagging strand. A positive RFD therefore means the fork moves toward
#' higher coordinates, and the plus (reference) strand is copied as the
#' lagging strand. So `"lagging"` means the plus strand is the
#' lagging-strand template, and `"leading"` means it is the leading-strand
#' template. Bins where the RFD is 0 carry no direction and are left out.
#'
#' `data-raw/rep_range.R` in the package source builds the object.
#'
#' @docType data
#'
#' @usage data(rep_range)
#'
#' @format A \code{"GRanges"} object with 280,770 ranges, tagged with
#' genome `"hg19"`. Its one metadata column, `X`, is `"leading"` or
#' `"lagging"`. See [annotate_replication_strand()].
#'
#' @keywords datasets
#'
#' @references Sriramachandran, A. M. et al. (2020)
#' Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns,
#' Single-Strand Breaks, and Lesions by GLOE-Seq. ([Molecular Cell]
#' (doi:10.1016/j.molcel.2020.03.027))
#'
#' @source GEO sample GSM4305465 (`WT_siLIG1_rep1`, 10 kb RFD), series
#' <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225>
"rep_range"

#' Replication strand of the hg38 genome as a GRanges object
#'
#' [rep_range] lifted over from hg19 to hg38 (GRCh38) with the UCSC
#' hg19ToHg38 chain. Use it with variants called on hg38. The labels mean
#' the same as in [rep_range]: `"lagging"` means the plus strand is the
#' lagging-strand template.
#'
#' In regions that are inverted between hg19 and hg38, the plus strand of
#' hg38 is the minus strand of hg19, so the label is swapped. Positions
#' where bins with different labels land on each other are left out, as
#' are the 95 bins that don't map. Neighbouring ranges with the same label
#' are merged, so the ranges are segments rather than 10 kb bins.
#'
#' `data-raw/rep_range.R` in the package source builds the object.
#'
#' @docType data
#'
#' @usage data(rep_range_hg38)
#'
#' @format A \code{"GRanges"} object with 89,121 ranges, tagged with genome
#' `"hg38"`. Its one metadata column, `X`, is `"leading"` or `"lagging"`.
#' See [annotate_replication_strand()].
#'
#' @keywords datasets
#'
#' @references Sriramachandran, A. M. et al. (2020)
#' Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns,
#' Single-Strand Breaks, and Lesions by GLOE-Seq. ([Molecular Cell]
#' (doi:10.1016/j.molcel.2020.03.027))
#'
#' @source Lifted over from [rep_range]; GEO sample GSM4305465, series
#' <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225>
"rep_range_hg38"
