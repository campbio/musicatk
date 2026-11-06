library(musicatk)

# A variant table in the column layout create_musica_from_variants() expects
variant_table <- function(chr = "chr1", start, end = start, ref, alt,
                          sample = "S1", ...) {
  data.frame(
    chr = chr, start = start, end = end, ref = ref, alt = alt,
    sample = sample, ...
  )
}

test_that("create_musica_from_variants accepts matrices and factors", {
  # Issue #114
  g <- select_genome("hg38")
  tab <- variant_table(
    start = c(5L, 2000L), ref = c("C", "G"), alt = c("A", "T"),
    sample = c("S1", "S2")
  )
  from_matrix <- create_musica_from_variants(as.matrix(tab), g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_identical(variants(from_matrix)$start, c(5L, 2000L))
  expect_identical(variants(from_matrix)$Variant_Type, c("SBS", "SBS"))

  factors <- variant_table(
    start = c(5L, 2000L), ref = c("C", "G"), alt = c("A", "T"),
    sample = c("S1", "S2"), stringsAsFactors = TRUE
  )
  from_factors <- create_musica_from_variants(factors, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_identical(variants(from_factors)$ref, c("C", "G"))
  expect_identical(variants(from_factors)$Variant_Type, c("SBS", "SBS"))
})

test_that("Lowercase bases are read as uppercase, not as compound", {
  # Issue #123
  g <- select_genome("hg38")
  tab <- variant_table(
    start = c(1000L, 2000L), ref = c("c", "g"), alt = c("t", "gtt")
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_identical(variants(musica)$ref, c("C", "-"))
  expect_identical(variants(musica)$alt, c("T", "TT"))
})

test_that("VCF-style deletions get the coordinates of the deleted bases", {
  # Issue #108: the same deletion and insertion written VCF- and MAF-style
  g <- select_genome("hg38")
  vcf_style <- variant_table(
    start = c(100L, 200L), end = c(101L, 200L),
    ref = c("CT", "G"), alt = c("C", "GTT")
  )
  maf_style <- variant_table(
    start = c(101L, 200L), end = c(101L, 201L),
    ref = c("T", "-"), alt = c("-", "TT")
  )
  from_vcf <- variants(create_musica_from_variants(vcf_style, g,
    check_ref_bases = FALSE, verbose = FALSE
  ))
  from_maf <- variants(create_musica_from_variants(maf_style, g,
    check_ref_bases = FALSE, verbose = FALSE
  ))
  expect_identical(from_vcf[1, ], from_maf[1, ])
  # Insertions keep 'start' as the base before the insertion
  expect_identical(from_vcf$start[2], 200L)
  expect_identical(from_vcf$ref[2], "-")
  expect_identical(from_vcf$alt[2], "TT")
})

test_that("Indels whose first base doesn't match are removed as complex", {
  # Issue #109
  g <- select_genome("hg38")
  tab <- variant_table(
    start = c(100L, 200L, 300L, 400L),
    end = c(100L, 202L, 300L, 401L),
    ref = c("A", "AGT", "G", "CA"),
    alt = c("GGT", "C", "GTT", "C")
  )
  expect_message(
    musica <- create_musica_from_variants(tab, g,
      check_ref_bases = FALSE, verbose = FALSE
    ),
    regexp = "2 complex"
  )
  expect_identical(variants(musica)$start, c(300L, 401L))
  expect_identical(variants(musica)$ref, c("-", "A"))
  expect_identical(variants(musica)$alt, c("TT", "-"))
})

test_that("Adjacent SBS are merged within a sample only", {
  # Issue #107
  g <- select_genome("hg38")
  tab <- variant_table(
    start = c(1000L, 1001L, 2001L, 2000L, 3000L, 3001L, 3002L),
    ref = c("C", "G", "T", "A", "A", "C", "G"),
    alt = c("A", "T", "G", "C", "C", "G", "T"),
    sample = c("S1", "S2", "S1", "S1", "S1", "S1", "S1")
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  v <- variants(musica)
  # Adjacent positions in different samples stay SBS
  separate <- v$start %in% c(1000L, 1001L)
  expect_identical(v$Variant_Type[separate], c("SBS", "SBS"))
  # Unsorted adjacent SBS in one sample become one DBS
  dbs <- v[v$Variant_Type == "DBS", ]
  expect_identical(nrow(dbs), 1L)
  expect_identical(c(dbs$start, dbs$end), c(2000L, 2001L))
  expect_identical(c(dbs$ref, dbs$alt), c("AT", "CG"))
  # A run of three is one MNV with all three bases
  mnv <- v[v$Variant_Type == "MNV", ]
  expect_identical(nrow(mnv), 1L)
  expect_identical(c(mnv$start, mnv$end), c(3000L, 3002L))
  expect_identical(c(mnv$ref, mnv$alt), c("ACG", "CGT"))
  expect_identical(nrow(v), 4L)
})

test_that("Variants on contigs missing from the genome are dropped", {
  # Issue #115
  g <- select_genome("hg38")
  base <- as.character(BSgenome::getSeq(g, "chr1", 1000000, 1000000))
  tab <- variant_table(
    chr = c("chr1", "chrFoo"), start = c(1000000L, 500L),
    ref = c(base, "A"), alt = c(setdiff(c("A", "C"), base)[1], "C")
  )
  expect_warning(
    musica <- create_musica_from_variants(tab, g, verbose = FALSE),
    regexp = "chrFoo"
  )
  expect_identical(variants(musica)$chr, "chr1")
})

test_that("create_musica_from_counts matches rows to motifs by name", {
  # Issue #116
  data(musica)
  counts <- get_count_table(extract_count_tables(musica)$SBS96)
  reversed <- counts[rev(seq_len(nrow(counts))), ]
  from_reversed <- create_musica_from_counts(reversed, "SBS96")
  expect_identical(get_count_table(tables(from_reversed)$SBS96), counts)

  misnamed <- counts
  rownames(misnamed)[1] <- "A[C>A]A"
  expect_error(
    create_musica_from_counts(misnamed, "SBS96"),
    regexp = "row names"
  )

  unnamed <- unname(counts)
  colnames(unnamed) <- colnames(counts)
  expect_message(
    from_unnamed <- create_musica_from_counts(unnamed, "SBS96"),
    regexp = "no row names"
  )
  expect_identical(get_count_table(tables(from_unnamed)$SBS96), counts)

  data(dbs_musica)
  dbs <- get_count_table(extract_count_tables(dbs_musica)$DBS78)
  from_dbs <- create_musica_from_counts(dbs[rev(seq_len(nrow(dbs))), ], "DBS")
  expect_identical(get_count_table(tables(from_dbs)$DBS78), dbs)
})
