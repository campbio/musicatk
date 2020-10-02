#' COSMIC v2 SBS96 Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts.
#'
#' @docType data
#'
#' @usage data(cosmic_v2_sigs)
#'
#' @format An object of class \code{musica_result}
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
#' @format An object of class \code{musica_result}.
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
#' @format An object of class \code{musica_result}.
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
#' @format An object of class \code{musica_result}.
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
#' @format An object of class \code{musica_result}.
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

#' Replication Timing Data as GRanges Object
#'
#' Supplementary data converted from bigWig to bedgraph to GRanges, with low
#' RFD indicating the leading strand and high RFD indicating lagging strand and
#' removing uninformative zero RFD intervals. Timing data is 10kb bins from
#' a colon cancer sample.
#'
#' @docType data
#'
#' @usage data(rep_range)
#'
#' @format An object of class \code{"GRanges"}; see
#' [annotate_replication_strand()].
#'
#' @keywords datasets
#'
#' @references Sriramachandran, A. M. et al. (2020)
#' Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns,
#' Single-Strand Breaks, and Lesions by GLOE-Seq. ([Molecular Cell]
#' (doi:10.1016/j.molcel.2020.03.027))
#'
#' @source GEO, <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225>
"rep_range"
