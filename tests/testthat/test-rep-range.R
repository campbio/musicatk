library(musicatk)

# Replication strand of SBS at the given positions, using rep_obj. Each
# position is its own sample, so neighbouring positions don't become a DBS.
rep_strand_at <- function(chr, pos, genome, rep_obj) {
  tab <- data.frame(
    chr = chr, start = pos, end = pos, ref = "C", alt = "A",
    sample = paste0("S", seq_along(pos))
  )
  musica <- create_musica_from_variants(tab, select_genome(genome),
    check_ref_bases = FALSE, verbose = FALSE
  )
  annotate_replication_strand(musica, rep_obj, build_table = FALSE)
  variants(musica)$Replication_Strand
}

test_that("Replication strand objects record their genome build", {
  # Issue #220
  data_items <- utils::data(package = "musicatk")$results[, "Item"]
  expect_true("rep_range_hg38" %in% data_items)

  expect_identical(unique(GenomeInfoDb::genome(rep_range)), "hg19")
  expect_identical(unique(GenomeInfoDb::genome(rep_range_hg38)), "hg38")
  expect_false(anyNA(GenomeInfoDb::seqlengths(rep_range)))
  expect_false(anyNA(GenomeInfoDb::seqlengths(rep_range_hg38)))
})

test_that("Replication strand bins don't overlap", {
  # Issue #220: rep_range had 0-based starts, so bins shared their edges
  expect_true(GenomicRanges::isDisjoint(rep_range))
  expect_true(GenomicRanges::isDisjoint(rep_range_hg38))
})

test_that("rep_range gives a position on a bin edge one label", {
  # hg19 chr8:220001-230000 is leading and chr8:230001-240000 lagging
  expect_identical(
    rep_strand_at("chr8", c(225000L, 230000L, 230001L), "hg19", rep_range),
    c("leading", "leading", "lagging")
  )
})

test_that("rep_range_hg38 keeps labels where the liftover is direct", {
  # hg19 chr8:80001-90000 (lagging) lifts to hg38 chr8:130001-140000
  expect_identical(
    rep_strand_at("chr8", 85000L, "hg19", rep_range),
    "lagging"
  )
  expect_identical(
    rep_strand_at("chr8", 135000L, "hg38", rep_range_hg38),
    "lagging"
  )
})

test_that("rep_range_hg38 swaps labels in inverted regions", {
  # hg19 chr2:90410001-90420000 (lagging) lifts to hg38
  # chr2:90354037-90364036 on the minus strand
  expect_identical(
    rep_strand_at("chr2", 90415000L, "hg19", rep_range),
    "lagging"
  )
  expect_identical(
    rep_strand_at("chr2", 90359036L, "hg38", rep_range_hg38),
    "leading"
  )
})

test_that("annotate_replication_strand accepts a build name", {
  # "hg19" and "hg38" pick the bundled object for that build
  expect_identical(
    rep_strand_at("chr8", c(85000L, 230000L), "hg19", "hg19"),
    rep_strand_at("chr8", c(85000L, 230000L), "hg19", rep_range)
  )
  expect_identical(
    rep_strand_at("chr2", c(90359036L, 135000L), "hg38", "hg38"),
    rep_strand_at("chr2", c(90359036L, 135000L), "hg38", rep_range_hg38)
  )
  expect_identical(
    rep_strand_at("chr8", 135000L, "hg38", "38"),
    "lagging"
  )
})

test_that("annotate_replication_strand rejects other rep_range values", {
  msg <- "rep_range must be a GRanges object, or \"hg19\" or \"hg38\""
  data(musica)
  expect_error(annotate_replication_strand(musica, "mm10"), msg,
    fixed = TRUE
  )
  expect_error(annotate_replication_strand(musica, c("hg19", "hg38")), msg,
    fixed = TRUE
  )
  expect_error(annotate_replication_strand(musica, data.frame()), msg,
    fixed = TRUE
  )
})
