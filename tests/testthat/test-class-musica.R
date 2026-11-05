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

test_that("Shipped musica objects satisfy the sample alignment invariant", {
  datasets <- c(
    "musica", "musica_annot", "musica_sbs96", "musica_sbs96_tiny",
    "dbs_musica", "indel_musica", "res", "res_annot"
  )
  for (dataset in datasets) {
    env <- new.env()
    data(list = dataset, envir = env)
    expect_true(validObject(get(dataset, envir = env)), label = dataset)
  }
  data(full_benchmark_example)
  expect_true(validObject(ground_truth(full_benchmark_example)))
})

test_that("musica validity rejects misaligned samples", {
  data(res_annot)
  samples <- sample_names(res_annot)
  sbs96 <- tables(res_annot)$SBS96
  reversed <- sbs96
  reversed@count_table <- get_count_table(sbs96)[, rev(samples)]

  expect_error(
    new("musica",
      variants = variants(res_annot),
      count_tables = list(SBS96 = reversed),
      sample_annotations = samp_annot(res_annot)
    ),
    regexp = "columns of count table 'SBS96'"
  )
  expect_error(
    tables(res_annot) <- list(SBS96 = reversed),
    regexp = "columns of count table 'SBS96'"
  )
  expect_error(
    new("musica",
      variants = variants(res_annot),
      sample_annotations = data.frame(Samples = samples[-1])
    ),
    regexp = "variant samples are missing"
  )
  expect_error(
    new("musica",
      sample_annotations = data.frame(Samples = c("a", "a"))
    ),
    regexp = "must be unique"
  )

  model <- get_model(res_annot, "result", "SBS96", "res_annot")
  model@exposures <- exposures(model)[, rev(samples)]
  results <- result_list(res_annot)
  results[["result"]]@modality[["SBS96"]][["res_annot"]] <- model
  expect_error(
    result_list(res_annot) <- results,
    regexp = "exposures of model 'res_annot'"
  )
})

test_that("result_model and result_collection validity", {
  data(res)
  model <- get_model(res, "result", "SBS96", "res")
  wrong_names <- exposures(model)
  rownames(wrong_names) <- c("x", "y", "z")
  expect_error(
    new("result_model",
      signatures = signatures(model), exposures = wrong_names
    ),
    regexp = "rownames of 'exposures' must match"
  )
  expect_error(
    new("result_collection",
      modality = S4Vectors::SimpleList(SBS96 = list(res = "not a model"))
    ),
    regexp = "result_model"
  )
})

test_that("Exposure setters and add_result() check sample alignment", {
  data(res_annot)
  samples <- sample_names(res_annot)
  model <- get_model(res_annot, "result", "SBS96", "res_annot")
  reversed <- exposures(model)[, rev(samples)]

  expect_error(
    exposures(res_annot, "result", "SBS96", "res_annot") <- reversed,
    regexp = "exposures of model 'res_annot'"
  )
  # add_result() matches columns by name and puts them in sample order
  added <- add_result(signatures(model), reversed, res_annot, "result",
    "SBS96", "reversed",
    make_copy = TRUE
  )
  expect_identical(
    exposures(added, "result", "SBS96", "reversed"),
    exposures(model)
  )
  renamed <- exposures(model)
  colnames(renamed)[1] <- "not_a_sample"
  expect_error(
    add_result(signatures(model), renamed, res_annot, "result", "SBS96",
      "renamed",
      make_copy = TRUE
    ),
    regexp = "exposures of model 'renamed'"
  )
  expect_error(
    add_result(signatures(model), unname(exposures(model)), res_annot,
      "result", "SBS96", "unnamed",
      make_copy = TRUE
    ),
    regexp = "exposures of model 'unnamed'"
  )
})

test_that("Count-only objects from older versions still have sample names", {
  data(musica_sbs96)
  sbs96 <- tables(musica_sbs96)$SBS96
  # Shape of a create_musica_from_counts() object saved before 2.7.2:
  # a count table, no variants and no sample annotations
  legacy <- new("musica")
  legacy@count_tables <- list(SBS96 = sbs96)
  expect_identical(sample_names(legacy), colnames(get_count_table(sbs96)))
  expect_true(validObject(legacy))
  samp_annot(legacy, "group") <- rep(c("a", "b"), c(3, 4))
  expect_identical(samp_annot(legacy)$group, rep(c("a", "b"), c(3, 4)))
})

test_that("subset_variants_by_samples() keeps every requested sample", {
  data(musica)
  samples <- sample_names(musica)[1:2]
  subset <- subset_variants_by_samples(musica, samples)
  expect_equal(nrow(subset), sum(variants(musica)$sample %in% samples))
  expect_setequal(as.character(unique(subset$sample)), samples)
})

test_that("subset_musica_by_annotation() reports values that are missing", {
  data(res_annot)
  expect_error(
    subset_musica_by_annotation(res_annot, "Tumor_Subtypes", c("Nope", "No")),
    regexp = "Nope, No not present in Tumor_Subtypes column",
    fixed = TRUE
  )
})

test_that("samp_annot<- accepts dates and variants<- accepts a data.frame", {
  data(musica)
  n <- length(sample_names(musica))
  dates <- rep(as.Date("2026-01-01"), n)
  samp_annot(musica, "date") <- dates
  expect_identical(samp_annot(musica)$date, dates)
  times <- rep(as.POSIXct("2026-01-01 10:00:00", tz = "UTC"), n)
  samp_annot(musica, "time") <- times
  expect_identical(samp_annot(musica)$time, times)
  expect_error(
    samp_annot(musica, "bad") <- matrix(1, n, 2),
    regexp = "must be a vector"
  )
  expect_error(
    samp_annot(musica, "bad") <- as.list(seq_len(n)),
    regexp = "must be a vector"
  )

  variants_df <- as.data.frame(variants(musica))
  variants(musica) <- variants_df
  expect_s3_class(variants(musica), "data.table")
  expect_equal(nrow(variants(musica)), nrow(variants_df))
})

test_that("name_signatures() returns the musica object with new names", {
  data(res)
  new_names <- c("smoking", "apobec", "unknown")
  renamed <- name_signatures(res, "res", new_names)
  expect_s4_class(renamed, "musica")
  expect_identical(
    colnames(signatures(renamed, "result", "SBS96", "res")), new_names
  )
  expect_identical(
    rownames(exposures(renamed, "result", "SBS96", "res")), new_names
  )

  # The caller's object is still renamed in place (see #205)
  data(res)
  name_signatures(res, "res", new_names)
  expect_identical(colnames(signatures(res, "result", "SBS96", "res")),
                   new_names)
})
