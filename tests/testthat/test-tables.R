library(musicatk)

test_that("IND83 columns follow sample_names(), matched by name", {
  data(indel_musica)
  g <- select_genome("hg38")
  samples <- as.character(unique(variants(indel_musica)$sample))
  forward <- new("musica",
    variants = variants(indel_musica),
    sample_annotations = data.frame(Samples = samples)
  )
  reverse <- new("musica",
    variants = variants(indel_musica),
    sample_annotations = data.frame(Samples = rev(samples))
  )
  build_standard_table(forward, g, "IND83")
  build_standard_table(reverse, g, "IND83")
  forward_counts <- get_count_table(tables(forward)$IND83)
  reverse_counts <- get_count_table(tables(reverse)$IND83)

  expect_identical(colnames(forward_counts), samples)
  expect_identical(colnames(reverse_counts), rev(samples))
  expect_equal(reverse_counts[, samples], forward_counts)
})

test_that("Count tables have no columns for stale sample factor levels", {
  g <- select_genome("hg38")
  data(musica_sbs96_tiny)
  tiny <- musica_sbs96_tiny
  expect_warning(
    build_standard_table(tiny, g, "SBS96", overwrite = TRUE),
    regexp = "Overwriting counts table"
  )
  expect_identical(
    colnames(get_count_table(tables(tiny)$SBS96)),
    sample_names(tiny)
  )

  data(dbs_musica)
  dbs <- dbs_musica
  expect_warning(
    build_standard_table(dbs, g, "DBS78", overwrite = TRUE),
    regexp = "Overwriting counts table"
  )
  expect_identical(
    colnames(get_count_table(tables(dbs)$DBS78)),
    sample_names(dbs)
  )
})

test_that("combine_count_tables() and custom tables align samples by name", {
  data(musica_sbs96)
  sbs96 <- tables(musica_sbs96)$SBS96
  counts <- get_count_table(sbs96)
  samples <- colnames(counts)

  musica <- new("musica",
    variants = variants(musica_sbs96),
    count_tables = list(SBS96 = sbs96),
    sample_annotations = data.frame(Samples = samples)
  )

  # Columns in another order are put back in sample order by name
  reversed_counts <- counts[, rev(samples)]
  expect_identical(
    .align_count_columns(musica, reversed_counts, "reversed"),
    counts
  )

  # A table in another order, as saved by older versions, is not bound by
  # position: combining stops and asks for the table to be rebuilt
  reversed <- sbs96
  reversed@count_table <- reversed_counts
  legacy <- musica
  legacy@count_tables$reversed <- reversed
  expect_error(
    combine_count_tables(legacy, c("SBS96", "reversed"), "combo"),
    regexp = "count table 'reversed'"
  )

  # Custom tables use sample_names() order
  custom <- new("musica",
    variants = variants(musica_sbs96),
    sample_annotations = data.frame(Samples = rev(samples))
  )
  build_custom_table(custom, "Variant_Type", "type_table")
  type_counts <- get_count_table(tables(custom)$type_table)
  expect_identical(colnames(type_counts), rev(samples))
  per_sample <- table(as.character(variants(musica_sbs96)$sample))
  expect_equal(colSums(type_counts), c(per_sample[rev(samples)]),
    ignore_attr = TRUE
  )

  # Tables whose samples differ from sample_names() are refused
  extra <- sbs96
  extra@count_table <- counts[, -1]
  wrong <- musica
  wrong@count_tables$extra <- extra
  expect_error(
    combine_count_tables(wrong, c("SBS96", "extra"), "combo"),
    regexp = "do not match the sample names"
  )
})
