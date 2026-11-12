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

#' COSMIC v3.0 SBS96 Exome Signatures Result Object (GRCh37)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 65 COSMIC v3.0 SBS signatures normalized to the
#' GRCh37 exome.
#'
#' @section Version:
#' This object is COSMIC v3.0 and will stay v3.0. COSMIC stopped publishing
#' exome-normalized catalogues after v3.0, so it is also the catalogue
#' \code{compare_cosmic_v3(sample_type = "exome")} uses by default. Signatures
#' added after v3.0 can only be matched against the genome catalogues, such as
#' \code{\link{cosmic_v3_6_sbs_sigs_grch37}}. For exome data,
#' \code{cosmic_version = "latest"} and \code{"3.0"} both use this object;
#' \code{"3.6"} is an error.
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
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.0 (May 2019), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_sbs_sigs_exome"

#' COSMIC v3.0 SBS96 Genome Signatures Result Object (GRCh37)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 67 COSMIC v3.0 SBS signatures normalized to the
#' GRCh37 genome.
#'
#' @section Version:
#' This object is COSMIC v3.0 and will stay v3.0, so code that uses it gives the
#' same results across musicatk releases. Newer COSMIC versions ship as separate
#' objects with the version in their names, such as
#' \code{\link{cosmic_v3_6_sbs_sigs_grch37}}. \code{compare_cosmic_v3()} uses
#' the newest version by default; set \code{cosmic_version = "3.0"} to compare
#' against this object.
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
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.0 (May 2019), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_sbs_sigs"

#' COSMIC v3.0 DBS Genome Signatures Result Object (GRCh37)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 11 COSMIC v3.0 DBS signatures normalized to the
#' GRCh37 genome.
#'
#' @section Version:
#' This object is COSMIC v3.0 and will stay v3.0, so code that uses it gives the
#' same results across musicatk releases. Newer COSMIC versions ship as separate
#' objects with the version in their names, such as
#' \code{\link{cosmic_v3_6_dbs_sigs_grch37}}. \code{compare_cosmic_v3()} uses
#' the newest version by default; set \code{cosmic_version = "3.0"} to compare
#' against this object.
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
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.0 (May 2019), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_dbs_sigs"

#' COSMIC v3.0 Indel Genome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 17 COSMIC v3.0 Indel signatures. COSMIC publishes
#' one Indel catalogue per version, used for any build.
#'
#' @section Version:
#' This object is COSMIC v3.0 and will stay v3.0, so code that uses it gives the
#' same results across musicatk releases. Newer COSMIC versions ship as separate
#' objects with the version in their names, such as
#' \code{\link{cosmic_v3_6_indel_sigs}}. \code{compare_cosmic_v3()} uses the
#' newest version by default; set \code{cosmic_version = "3.0"} to compare
#' against this object.
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
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.0 (May 2019), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_indel_sigs"

#' COSMIC v3.6 SBS96 Genome Signatures Result Object (GRCh37)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 101 COSMIC v3.6 SBS signatures normalized to the
#' GRCh37 trinucleotide frequencies. \code{\link{cosmic_v3_6_sbs_sigs_grch38}}
#' has the same signatures normalized to GRCh38. Use the object that matches the
#' build the variants were called against. \code{compare_cosmic_v3()} picks it
#' from \code{genome_build}.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_6_sbs_sigs_grch37)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.6 (May 2026), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_6_sbs_sigs_grch37"

#' COSMIC v3.6 SBS96 Genome Signatures Result Object (GRCh38)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 101 COSMIC v3.6 SBS signatures normalized to the
#' GRCh38 trinucleotide frequencies. \code{\link{cosmic_v3_6_sbs_sigs_grch37}}
#' has the same signatures normalized to GRCh37. Use the object that matches the
#' build the variants were called against. \code{compare_cosmic_v3()} picks it
#' from \code{genome_build}.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_6_sbs_sigs_grch38)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.6 (May 2026), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_6_sbs_sigs_grch38"

#' COSMIC v3.6 DBS Genome Signatures Result Object (GRCh37)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 22 COSMIC v3.6 DBS signatures normalized to the
#' GRCh37 dinucleotide frequencies. \code{\link{cosmic_v3_6_dbs_sigs_grch38}}
#' has the same signatures normalized to GRCh38. Use the object that matches the
#' build the variants were called against. \code{compare_cosmic_v3()} picks it
#' from \code{genome_build}.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_6_dbs_sigs_grch37)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.6 (May 2026), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_6_dbs_sigs_grch37"

#' COSMIC v3.6 DBS Genome Signatures Result Object (GRCh38)
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 22 COSMIC v3.6 DBS signatures normalized to the
#' GRCh38 dinucleotide frequencies. \code{\link{cosmic_v3_6_dbs_sigs_grch37}}
#' has the same signatures normalized to GRCh37. Use the object that matches the
#' build the variants were called against. \code{compare_cosmic_v3()} picks it
#' from \code{genome_build}.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_6_dbs_sigs_grch38)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.6 (May 2026), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_6_dbs_sigs_grch38"

#' COSMIC v3.6 Indel Genome Signatures Result Object
#'
#' Data from COSMIC formatted to be used for prediction with individual tumors
#' and cohorts. Contains the 25 COSMIC v3.6 Indel signatures. COSMIC publishes
#' one Indel catalogue per version, used for any build.
#'
#' @docType data
#'
#' @usage data(cosmic_v3_6_indel_sigs)
#'
#' @format An object of class \code{result_model}.
#' See [predict_exposure()].
#'
#' @keywords datasets
#'
#' @references Alexandrov, L.B., Kim, J., Haradhvala, N.J. et al. (2020)
#' The repertoire of mutational signatures in human cancer. Nature 578, 94-101
#' ([Nature](https://doi.org/10.1038/s41586-020-1943-3))
#'
#' @source COSMIC v3.6 (May 2026), <https://cancer.sanger.ac.uk/signatures/>
"cosmic_v3_6_indel_sigs"

#' Replication strand of the hg19 genome as a GRanges object
#'
#' Leading and lagging replication strand on the hg19 (GRCh37) genome
#' build, from the replication fork directionality (RFD) that GLOE-seq
#' measured in 10 kb bins in HCT116 colon cancer cells with ligase 1
#' depleted. Neighbouring bins with the same RFD are stored as one range,
#' so a few ranges are 20 to 50 kb. Use it only with variants called on
#' hg19; for hg38 use \code{\link{rep_range_hg38}}.
#'
#' Each bin is labelled by the sign of its RFD, (REV - FWD) / (REV + FWD).
#' In these cells the breaks GLOE-seq captures lie on the newly made
#' lagging strand. A positive RFD therefore means the fork moves toward
#' higher coordinates, and the plus (reference) strand is copied as the
#' lagging strand. So \code{"lagging"} means the plus strand is the
#' lagging-strand template, and \code{"leading"} means it is the
#' leading-strand template. Bins where the RFD is 0 carry no direction and
#' are left out.
#'
#' \code{data-raw/rep_range.R} in the package source builds the object.
#'
#' @docType data
#'
#' @usage data(rep_range)
#'
#' @format A \code{"GRanges"} object with 280,770 ranges, tagged with
#' genome \code{"hg19"}. Its one metadata column, \code{X}, is
#' \code{"leading"} or \code{"lagging"}. See
#' \code{\link{annotate_replication_strand}}.
#'
#' @keywords datasets
#'
#' @references Sriramachandran, A. M. et al. (2020)
#' Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns,
#' Single-Strand Breaks, and Lesions by GLOE-Seq. ([Molecular Cell]
#' (doi:10.1016/j.molcel.2020.03.027))
#'
#' @source GEO sample GSM4305465 (\code{WT_siLIG1_rep1}, 10 kb RFD), series
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225}
"rep_range"

#' Replication strand of the hg38 genome as a GRanges object
#'
#' \code{\link{rep_range}} lifted over from hg19 to hg38 (GRCh38) with the
#' UCSC hg19ToHg38 chain. Use it with variants called on hg38. The labels
#' mean the same as in \code{\link{rep_range}}: \code{"lagging"} means the
#' plus strand is the lagging-strand template.
#'
#' In regions that are inverted between hg19 and hg38, the plus strand of
#' hg38 is the minus strand of hg19, so the label is swapped. Positions
#' where bins with different labels land on each other are left out, as
#' are the 95 bins that don't map and hg38 bases with no hg19 counterpart.
#' Neighbouring ranges with the same label are merged, so the ranges are
#' segments rather than 10 kb bins.
#'
#' \code{data-raw/rep_range.R} in the package source builds the object.
#'
#' @docType data
#'
#' @usage data(rep_range_hg38)
#'
#' @format A \code{"GRanges"} object with 89,121 ranges, tagged with genome
#' \code{"hg38"}. Its one metadata column, \code{X}, is \code{"leading"} or
#' \code{"lagging"}. See \code{\link{annotate_replication_strand}}.
#'
#' @keywords datasets
#'
#' @references Sriramachandran, A. M. et al. (2020)
#' Genome-wide Nucleotide-Resolution Mapping of DNA Replication Patterns,
#' Single-Strand Breaks, and Lesions by GLOE-Seq. ([Molecular Cell]
#' (doi:10.1016/j.molcel.2020.03.027))
#'
#' @source Lifted over from \code{\link{rep_range}}; GEO sample
#' GSM4305465, series
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134225}
"rep_range_hg38"
