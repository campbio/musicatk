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
  # Variants whose sample factor keeps levels for samples that are gone,
  # as left behind by subsetting in older versions
  add_stale_levels <- function(musica) {
    variants <- variants(musica)
    present <- sample_names(musica)
    variants$sample <- factor(as.character(variants$sample),
      levels = c("stale_1", present, "stale_2")
    )
    new("musica",
      variants = variants,
      sample_annotations = data.frame(Samples = present)
    )
  }

  data(musica_sbs96_tiny)
  tiny <- add_stale_levels(musica_sbs96_tiny)
  build_standard_table(tiny, g, "SBS96")
  expect_identical(
    colnames(get_count_table(tables(tiny)$SBS96)),
    sample_names(tiny)
  )

  data(dbs_musica)
  dbs <- add_stale_levels(dbs_musica)
  build_standard_table(dbs, g, "DBS78")
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

test_that(".subset_count_tables() selects samples by name, not factor code", {
  data(musica_sbs96)
  samples <- sample_names(musica_sbs96)
  # A factor whose integer codes point at other columns
  wanted <- factor(samples[c(6, 7)], levels = rev(samples))
  subset <- .subset_count_tables(musica_sbs96, wanted)
  expect_identical(
    colnames(get_count_table(subset$SBS96)),
    samples[c(6, 7)]
  )
})

test_that("SBS tables warn about variants whose motif isn't in the schema", {
  # Issue #135
  g <- select_genome("hg38")
  # hg38 chr1:43430030 is G, so the first ref doesn't match the genome;
  # chr1:5000 lies in the N gap at the start of chr1, so its flanks are N
  tab <- data.frame(
    chr = "chr1", start = c(43430030L, 43430032L, 5000L),
    end = c(43430030L, 43430032L, 5000L), ref = "C", alt = "A",
    sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_warning(
    build_standard_table(musica, g, "SBS96"),
    regexp = "2 of 3 variants have motifs not in the SBS96 schema"
  )
  expect_equal(sum(get_count_table(tables(musica)$SBS96)), 1)

  dat <- variants(musica)
  dat$Transcript_Strand <- "T"
  variants(musica) <- dat
  expect_warning(
    build_standard_table(musica, g, "SBS192", "Transcript_Strand"),
    regexp = "1 of 3 variants have motifs not in the SBS192_Trans schema"
  )
})
