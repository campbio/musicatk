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
