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
