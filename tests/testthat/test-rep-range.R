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
  for (build in c("38", "HG38", "GRCh38")) {
    expect_identical(rep_strand_at("chr8", 135000L, "hg38", build), "lagging")
  }
  expect_identical(rep_strand_at("chr8", 85000L, "hg19", "GRCh37"), "lagging")
})

test_that("annotate_replication_strand rejects other rep_range values", {
  msg <- "Replication strand data is bundled for hg19 and hg38 only"
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

test_that("A>* and G>* SBS take the label of the opposite strand", {
  # Issue #222: rep_range labels the plus strand, but SBS192 counts an
  # A>* or G>* SBS on the minus strand, where its pyrimidine is.
  # hg19 chr8:80001-90000 is lagging, so the minus strand there is leading.
  tab <- data.frame(
    chr = "chr8", start = 85000L, end = 85000L,
    ref = c("C", "T", "G", "A"), alt = c("A", "C", "T", "G"),
    sample = paste0("S", 1:4)
  )
  musica <- create_musica_from_variants(tab, select_genome("hg19"),
    check_ref_bases = FALSE, verbose = FALSE
  )
  annotate_replication_strand(musica, rep_range, build_table = FALSE)
  expect_identical(
    variants(musica)$Replication_Strand,
    c("lagging", "lagging", "leading", "leading")
  )
})

test_that("SBS192_Rep counts a G>T in a lagging bin as C>A leading", {
  # Issue #222
  g <- select_genome("hg19")
  tab <- data.frame(
    chr = "chr8", start = 85000L, end = 85000L, ref = "G", alt = "T",
    sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  annotate_replication_strand(musica, rep_range, build_table = FALSE)
  build_standard_table(musica, g, "SBS192", "Replication_Strand")

  counts <- get_count_table(tables(musica)$SBS192_Rep)
  counted <- rownames(counts)[counts[, "S1"] > 0]
  expect_identical(length(counted), 1L)
  expect_match(counted, "^C>A_.C._leading$")
})
