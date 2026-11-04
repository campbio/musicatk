library(musicatk)

test_that("sample_names() reads the Samples column exactly, as character", {
  variants <- data.table::data.table(
    chr = "chr1", start = 10, end = 10, ref = "C", alt = "A",
    sample = factor(c("s2", "s1")), Variant_Type = "SBS"
  )
  annot <- data.frame(
    Samples = factor(c("s2", "s1"), levels = c("s2", "s1")),
    Sample = c("user_a", "user_b")
  )
  musica <- new("musica", variants = variants, sample_annotations = annot)
  expect_identical(sample_names(musica), c("s2", "s1"))
})

test_that("create_musica_from_counts() sets sample names for annotation", {
  data(musica_sbs96)
  counts <- get_count_table(tables(musica_sbs96)$SBS96)
  musica <- create_musica_from_counts(counts, "SBS96")
  expect_identical(sample_names(musica), colnames(counts))

  samp_annot(musica, "group") <- rep(c("a", "b"), c(3, 4))
  expect_identical(samp_annot(musica)$group, rep(c("a", "b"), c(3, 4)))

  unnamed <- unname(counts)
  expect_error(
    create_musica_from_counts(unnamed, "SBS96"),
    regexp = "'x' must have column names"
  )
  duplicated_names <- counts
  colnames(duplicated_names)[2] <- colnames(duplicated_names)[1]
  expect_error(
    create_musica_from_counts(duplicated_names, "SBS96"),
    regexp = "must be unique"
  )
})

test_that("sample_names() without annotations keeps the variant level order", {
  variants <- data.table::data.table(
    chr = "chr1", start = 10, end = 10, ref = "C", alt = "A",
    sample = factor(c("b", "a", "b"), levels = c("stale", "b", "a")),
    Variant_Type = "SBS"
  )
  musica <- new("musica", variants = variants)
  expect_identical(sample_names(musica), c("b", "a"))
})

test_that("Subsetting a musica object also subsets its results", {
  data(res_annot)
  lung <- subset_musica_by_annotation(res_annot, "Tumor_Subtypes", "Lung")
  lung_exposures <- exposures(lung, "result", "SBS96", "res_annot")
  expect_identical(colnames(lung_exposures), sample_names(lung))
  expect_length(sample_names(lung), 4)
  expect_identical(levels(variants(lung)$sample), sample_names(lung))

  counts <- get_count_table(tables(res_annot)$SBS96)
  keep <- colnames(counts)[colSums(counts) >= 100]
  high <- subset_musica_by_counts(res_annot, "SBS96", 100)
  expect_setequal(sample_names(high), keep)
  expect_identical(
    colnames(exposures(high, "result", "SBS96", "res_annot")),
    sample_names(high)
  )
})
