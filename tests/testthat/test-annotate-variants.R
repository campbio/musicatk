library(musicatk)

# On hg38, chr1:43430028-43430032 reads CAGTC on the plus strand
flank_variants <- function() {
  data.frame(
    chr = "chr1", start = c(43430030L, 43430032L),
    end = c(43430030L, 43430032L), ref = c("G", "C"), alt = c("T", "A"),
    sample = c("S1", "S2")
  )
}

test_that("add_flank_to_variants reads flanks on the pyrimidine strand", {
  # Issue #104
  g <- select_genome("hg38")
  musica <- create_musica_from_variants(flank_variants(), g, verbose = FALSE)
  add_flank_to_variants(musica, g, -2, -1, build_table = FALSE)
  add_flank_to_variants(musica, g, 1, 2, build_table = FALSE)
  dat <- variants(musica)
  purine <- dat$ref == "G"

  # The C strand of CAGTC is GACTG: GA lies 5' of the mutated C, TG 3' of it
  expect_identical(dat$lflank_2_to_1[purine], "GA")
  expect_identical(dat$rflank_1_to_2[purine], "TG")
  # Pyrimidine reference: plus-strand bases 5' of the C
  expect_identical(dat$lflank_2_to_1[!purine], "GT")
})

test_that("Intergenic SBS have no transcript strand", {
  # Issue #105
  data(musica)
  annotate_transcript_strand(musica, "hg38")
  dat <- variants(musica)
  sbs_strand <- dat$Transcript_Strand[dat$Variant_Type == "SBS"]

  expect_true(any(is.na(sbs_strand)))
  expect_true(all(sbs_strand %in% c("T", "U", NA)))
  counts <- get_count_table(tables(musica)$Transcript_Strand)
  expect_equal(sum(counts), sum(!is.na(sbs_strand)))
})

test_that("annotate_transcript_strand accepts a TxDb object", {
  # Issue #117
  data(musica)
  by_name <- musica
  by_txdb <- musica
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  annotate_transcript_strand(by_name, "hg38", build_table = FALSE)
  annotate_transcript_strand(by_txdb, txdb, build_table = FALSE)

  expect_identical(
    variants(by_txdb)$Transcript_Strand,
    variants(by_name)$Transcript_Strand
  )
})

test_that("annotate_replication_strand works when the last rows aren't SBS", {
  # Issue #118
  g <- select_genome("hg38")
  data(rep_range)
  # rep_range calls chr8:80000-90000 lagging
  tab <- data.frame(
    chr = "chr8", start = c(85000L, 85010L), end = c(85000L, 85010L),
    ref = c("C", "-"), alt = c("A", "A"), sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_identical(variants(musica)$Variant_Type, c("SBS", "INS"))

  annotate_replication_strand(musica, rep_range, build_table = FALSE)
  expect_identical(variants(musica)$Replication_Strand, c("lagging", NA))
})

test_that("annotate_variant_length gives signed indel lengths", {
  # Issue #119
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1", start = c(1000L, 2000L, 3000L, 4000L),
    end = c(1000L, 2001L, 3000L, 4002L),
    ref = c("C", "AC", "-", "ATG"), alt = c("A", "GT", "AT", "-"),
    sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  annotate_variant_length(musica)
  dat <- variants(musica)
  lengths <- stats::setNames(dat$Variant_Length, dat$Variant_Type)

  expect_equal(lengths[c("SBS", "DBS", "INS", "DEL")],
    c(SBS = 1, DBS = 2, INS = 2, DEL = -3)
  )
})
