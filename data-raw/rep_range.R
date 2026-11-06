# Build the replication strand objects rep_range (hg19) and rep_range_hg38
# (#220). Run from the package root: Rscript data-raw/rep_range.R
#
# Source: GEO GSM4305465 (GSE134225), GLOE-seq of HCT116 cells with
# ligase 1 depleted (Sriramachandran et al. 2020, Mol Cell), mapped to
# hg19. The 10 kb file holds the replication fork directionality,
# RFD = (REV - FWD) / (REV + FWD), where FWD and REV count breaks on the
# plus and minus strands. Without ligase 1 the breaks are the unsealed
# nicks between Okazaki fragments on the newly made lagging strand. So
# RFD > 0 means the fork moves toward higher coordinates, and the plus
# (reference) strand is the lagging-strand template. Each bin is labelled
# "lagging" (RFD > 0) or "leading" (RFD < 0), for the plus strand. Bins
# where RFD = 0 carry no direction and are dropped.
#
# The original script for rep_range isn't in the repository. This script
# rebuilds it and stops unless the result matches the committed object,
# whose starts were 0-based before this rebuild.
#
# rep_range_hg38 is lifted over with the UCSC hg19ToHg38 chain:
# - Bins are lifted as "+" ranges. A piece that lands on "-" is in a region
#   inverted between the builds, where the fork direction reverses relative
#   to the coordinates, so its label is swapped.
# - Abutting pieces of the same bin are merged.
# - Where pieces of different bins overlap, the overlap is kept only if
#   their labels agree. Bins that don't map are dropped.
# - Neighbouring ranges with the same label are merged, so rep_range_hg38
#   holds segments rather than 10 kb bins.
#
# Needs rtracklayer, which isn't a package dependency because only this
# script uses it.

bw_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4305nnn/GSM4305465/suppl/",
  "GSM4305465_WT_siLIG1_rep1_rfd_10kb.bw"
)
chain_url <- paste0(
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
  "hg19ToHg38.over.chain.gz"
)
keep_chr <- paste0("chr", c(1:22, "X", "Y", "M"))

swap_label <- function(x) {
  ifelse(x == "leading", "lagging", "leading")
}

set_build <- function(gr, g) {
  gr <- GenomeInfoDb::keepSeqlevels(gr,
    intersect(GenomeInfoDb::seqlevels(g), GenomeInfoDb::seqlevels(gr)),
    pruning.mode = "coarse"
  )
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(g)[
    GenomeInfoDb::seqlevels(gr)
  ]
  sort(gr)
}

# hg19 ----------------------------------------------------------------------

bw_file <- file.path(tempdir(), basename(bw_url))
utils::download.file(bw_url, bw_file, mode = "wb")
rfd <- rtracklayer::import(bw_file)
keep <- as.character(GenomeInfoDb::seqnames(rfd)) %in% keep_chr &
  !is.na(rfd$score) & rfd$score != 0
rfd <- rfd[keep]

hg19 <- GenomicRanges::granges(rfd)
S4Vectors::mcols(hg19)$X <- ifelse(rfd$score > 0, "lagging", "leading")
g19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
hg19 <- set_build(hg19, g19)

# Check against the committed object
old_env <- new.env()
load("data/rep_range.rda", envir = old_env)
old <- old_env$rep_range
if (anyNA(GenomeInfoDb::genome(old))) {
  # The original object had 0-based starts
  BiocGenerics::start(old) <- BiocGenerics::start(old) + 1L
}
# The original object merged 10 pairs of neighbouring bins with the same
# label into one 20 kb range, so compare the positions each label covers
# rather than the ranges themselves
GenomeInfoDb::seqlevels(old) <- GenomeInfoDb::seqlevels(hg19)
GenomeInfoDb::seqinfo(old) <- GenomeInfoDb::seqinfo(hg19)
for (lab in c("leading", "lagging")) {
  old_cov <- GenomicRanges::reduce(old[old$X == lab])
  new_cov <- GenomicRanges::reduce(hg19[hg19$X == lab])
  stopifnot(identical(old_cov, new_cov))
}
rep_range <- hg19

# hg38 ----------------------------------------------------------------------

chain_gz <- file.path(tempdir(), basename(chain_url))
utils::download.file(chain_url, chain_gz, mode = "wb")
chain_file <- sub("\\.gz$", "", chain_gz)
chain_con <- gzfile(chain_gz)
writeLines(readLines(chain_con), chain_file)
close(chain_con)
chain <- rtracklayer::import.chain(chain_file)

src <- hg19
BiocGenerics::strand(src) <- "+"
GenomeInfoDb::seqinfo(src) <- GenomeInfoDb::Seqinfo(
  GenomeInfoDb::seqlevels(src)
)
lifted <- rtracklayer::liftOver(src, chain)

# Merge abutting pieces of each bin, keeping strand and chromosome apart
merged <- GenomicRanges::reduce(lifted, min.gapwidth = 1L)
pieces <- unlist(merged)
bin <- rep(seq_along(merged), lengths(merged))
label <- hg19$X[bin]
inverted <- as.character(BiocGenerics::strand(pieces)) == "-"
label[inverted] <- swap_label(label[inverted])
BiocGenerics::strand(pieces) <- "*"

# Split where pieces of different bins overlap; keep agreeing labels
dis <- GenomicRanges::disjoin(pieces, with.revmap = TRUE)
dis_labels <- lapply(dis$revmap, function(i) unique(label[i]))
agree <- lengths(dis_labels) == 1
kept <- GenomicRanges::granges(dis[agree])
kept_labels <- unlist(dis_labels[agree])

# Merge neighbouring ranges with the same label. Every position keeps its
# label, and the object is about a third the size.
hg38 <- unlist(GenomicRanges::reduce(S4Vectors::split(kept, kept_labels)))
S4Vectors::mcols(hg38)$X <- names(hg38)
names(hg38) <- NULL
g38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
hg38 <- set_build(hg38, g38)
rep_range_hg38 <- hg38

message(
  "rep_range: ", length(rep_range), " bins; rep_range_hg38: ",
  length(rep_range_hg38), " ranges; ", sum(lengths(lifted) == 0),
  " bins unmapped; ", sum(inverted), " inverted pieces; ", sum(!agree),
  " overlapping pieces dropped for conflicting labels"
)

save(rep_range, file = "data/rep_range.rda", compress = "xz")
save(rep_range_hg38, file = "data/rep_range_hg38.rda", compress = "xz")
