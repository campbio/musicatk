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
  # A hand-made table with 'end' equal to 'start' gets the same coordinates
  one_base_end <- vcf_style
  one_base_end$end[1] <- 100L
  from_one_base_end <- variants(create_musica_from_variants(one_base_end, g,
    check_ref_bases = FALSE, verbose = FALSE
  ))
  expect_identical(from_one_base_end[1, ], from_maf[1, ])
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
  # add_variant_type() gives the same types again
  expect_identical(musicatk:::add_variant_type(v)$Variant_Type, v$Variant_Type)

  # Two SBS at one position are not merged with their neighbour
  same_position <- variant_table(
    start = c(4000L, 4000L, 4001L), ref = c("A", "A", "C"),
    alt = c("C", "G", "T")
  )
  v <- variants(create_musica_from_variants(same_position, g,
    check_ref_bases = FALSE, verbose = FALSE
  ))
  expect_identical(v$Variant_Type, c("SBS", "SBS", "SBS"))
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

  # Also when the chromosome style check is turned off
  expect_warning(
    musica <- create_musica_from_variants(tab, g,
      check_ref_chromosomes = FALSE, verbose = FALSE
    ),
    regexp = "chrFoo"
  )
  expect_identical(variants(musica)$chr, "chr1")

  # A clear error when no variant is left
  unmapped <- variant_table(chr = "1", start = 100L, ref = "A", alt = "C")
  expect_error(
    suppressWarnings(create_musica_from_variants(unmapped, g,
      check_ref_chromosomes = FALSE, verbose = FALSE
    )),
    regexp = "No variants are left"
  )
})

test_that("create_musica_from_counts matches rows to motifs by name", {
  # Issue #116
  data(musica)
  counts <- get_count_table(extract_count_tables(musica)$SBS96)
  reversed <- counts[rev(seq_len(nrow(counts))), ]
  from_reversed <- create_musica_from_counts(reversed, "SBS96")
  expect_identical(get_count_table(tables(from_reversed)$SBS96), counts)

  misnamed <- counts
  rownames(misnamed)[1] <- "not_a_motif"
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

  duplicated_row <- counts
  rownames(duplicated_row)[2] <- rownames(counts)[1]
  expect_error(
    create_musica_from_counts(duplicated_row, "SBS96"),
    regexp = "Duplicated: C>A_ACA.*Missing: C>A_ACC"
  )

  # COSMIC/SigProfiler names (A[C>A]A, AC>CA) are converted, in any order
  bracket <- counts
  motif <- rownames(counts)
  rownames(bracket) <- paste0(
    substr(motif, 5, 5), "[", substr(motif, 1, 3), "]", substr(motif, 7, 7)
  )
  bracket <- bracket[order(rownames(bracket)), ]
  from_bracket <- create_musica_from_counts(bracket, "SBS96")
  expect_identical(get_count_table(tables(from_bracket)$SBS96), counts)

  data(dbs_musica)
  dbs <- get_count_table(extract_count_tables(dbs_musica)$DBS78)
  from_dbs <- create_musica_from_counts(dbs[rev(seq_len(nrow(dbs))), ], "DBS")
  expect_identical(get_count_table(tables(from_dbs)$DBS78), dbs)

  dbs_bracket <- dbs
  rownames(dbs_bracket) <- sub("NN_", "", rownames(dbs))
  dbs_bracket <- dbs_bracket[order(rownames(dbs_bracket)), ]
  from_dbs_bracket <- create_musica_from_counts(dbs_bracket, "DBS")
  expect_identical(get_count_table(tables(from_dbs_bracket)$DBS78), dbs)
})

test_that("MAFs without records or with only silent variants are read", {
  # Issue #241
  maf_file <- system.file("extdata", "public_TCGA_lung.maf",
    package = "musicatk"
  )
  lines <- readLines(maf_file)
  header <- startsWith(lines, "#") | startsWith(lines, "Hugo_Symbol")
  records <- lines[!header]
  classes <- vapply(strsplit(records, "\t", fixed = TRUE), "[",
    character(1), 9
  )

  empty <- tempfile(fileext = ".maf")
  writeLines(lines[header], empty)
  expect_warning(
    expect_null(extract_variants_from_maf_file(empty)),
    regexp = "no variant records"
  )

  silent <- tempfile(fileext = ".maf")
  writeLines(c(lines[header], records[classes == "Silent"][1:3]), silent)
  dt <- extract_variants_from_maf_file(silent)
  expect_identical(nrow(dt), 3L)

  # Also when the classes aren't ones read.maf() knows
  fields <- strsplit(records[classes == "Silent"][1:3], "\t", fixed = TRUE)
  renamed <- vapply(fields, function(x) {
    x[9] <- "Splice_Site_SNP"
    paste(x, collapse = "\t")
  }, character(1))
  renamed_maf <- tempfile(fileext = ".maf")
  writeLines(c(lines[header], renamed), renamed_maf)
  expect_identical(nrow(extract_variants_from_maf_file(renamed_maf)), 3L)

  # One empty MAF doesn't stop the others from being read
  lung <- extract_variants_from_maf_file(maf_file)
  expect_warning(
    both <- extract_variants(c(empty, maf_file), verbose = FALSE),
    regexp = "no variant records"
  )
  expect_identical(both, lung)

  # Gzipped MAFs are checked for records too
  skip_if_not_installed("R.utils")
  gz <- tempfile(fileext = ".maf.gz")
  con <- gzfile(gz, "w")
  writeLines(c(lines[header], records[1:2]), con)
  close(con)
  expect_identical(nrow(extract_variants_from_maf_file(gz)), 2L)
  gz_empty <- tempfile(fileext = ".maf.gz")
  con <- gzfile(gz_empty, "w")
  writeLines(lines[header], con)
  close(con)
  expect_warning(
    expect_null(extract_variants_from_maf_file(gz_empty)),
    regexp = "no variant records"
  )
})
