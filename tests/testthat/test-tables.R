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
