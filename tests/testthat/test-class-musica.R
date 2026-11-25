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

  samp_annot(musica, "group") <- rep(c("a", "b"), c(6, 6))
  expect_identical(samp_annot(musica)$group, rep(c("a", "b"), c(6, 6)))

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
  expect_length(sample_names(lung), 8)
  expect_identical(levels(variants(lung)$sample), sample_names(lung))

  counts <- get_count_table(tables(res_annot)$SBS96)
  keep <- colnames(counts)[colSums(counts) >= 300]
  high <- subset_musica_by_counts(res_annot, "SBS96", 300)
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
  samp_annot(legacy, "group") <- rep(c("a", "b"), c(6, 6))
  expect_identical(samp_annot(legacy)$group, rep(c("a", "b"), c(6, 6)))
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

test_that("samp_annot<- accepts Date and date-time vectors", {
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
  samp_annot(musica, "time_lt") <- as.POSIXlt(times)
  expect_identical(samp_annot(musica)$time_lt, times)
})

test_that("variants<- accepts a data.frame", {
  data(musica)
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

  # The first argument does not have to be a variable
  renamed <- name_signatures(identity(res), "res", new_names)
  expect_identical(
    colnames(signatures(renamed, "result", "SBS96", "res")), new_names
  )

  # The caller's object is still renamed in place (see #205)
  data(res)
  name_signatures(res, "res", new_names)
  expect_identical(colnames(signatures(res, "result", "SBS96", "res")),
                   new_names)
})

test_that("add_samp_annot() matches annotation rows to samples by name", {
  data(musica_sbs96)
  samples <- sample_names(musica_sbs96)
  annot <- read.table(system.file("extdata", "sample_annotations.txt",
                                  package = "musicatk"),
                      sep = "\t", header = TRUE)
  # Rows in a different order from the samples, and one extra sample
  annot <- annot[rev(seq_len(nrow(annot))), ]
  annot <- rbind(annot, data.frame(Sample_Names = "Not_A_Sample",
                                   Tumor_Subtypes = "Skin"))
  annot$Age <- seq_len(nrow(annot))
  expect_false(identical(annot$Sample_Names[seq_along(samples)], samples))

  result <- add_samp_annot(musica_sbs96, annot, "Sample_Names")
  expected <- annot[match(samples, annot$Sample_Names), ]
  expect_identical(sample_names(result), samples)
  expect_identical(colnames(samp_annot(result)),
                   c("Samples", "Tumor_Subtypes", "Age"))
  expect_identical(samp_annot(result)$Tumor_Subtypes,
                   expected$Tumor_Subtypes)
  expect_identical(samp_annot(result)$Age, expected$Age)
  expect_true(methods::validObject(result))

  # An existing annotation of the same name is replaced
  annot$Tumor_Subtypes <- "Replaced"
  result <- add_samp_annot(result, annot, "Sample_Names")
  expect_identical(colnames(samp_annot(result)),
                   c("Samples", "Tumor_Subtypes", "Age"))
  expect_true(all(samp_annot(result)$Tumor_Subtypes == "Replaced"))
})

test_that("add_samp_annot() gives NA, with a warning, to missing samples", {
  data(musica_sbs96)
  samples <- sample_names(musica_sbs96)
  annot <- data.frame(id = factor(samples[3:1]), group = c("c", "b", "a"))
  expect_warning(
    result <- add_samp_annot(musica_sbs96, annot, "id"),
    regexp = paste0("9 samples have no row in 'annotations'.*", samples[4])
  )
  expect_identical(samp_annot(result)$group,
                   c("a", "b", "c", rep(NA, length(samples) - 3)))
})

test_that("add_samp_annot(convert = TRUE) keeps number-like IDs (#287)", {
  data(musica_sbs96)
  counts <- get_count_table(tables(musica_sbs96)$SBS96)[, 1:3]
  colnames(counts) <- c("001", "002", "100000")
  musica <- create_musica_from_counts(counts, "SBS96")
  annot_file <- withr::local_tempfile(fileext = ".txt")
  writeLines(c("id\tgroup\tscore\tflag",
               "100000\tA\t3\tTRUE",
               "001\tB\t1.5\tFALSE",
               "002\t\t\tTRUE"), annot_file)
  # Every column read as character, so the sample IDs keep their form
  annot <- utils::read.delim(annot_file, colClasses = "character")

  result <- add_samp_annot(musica, annot, "id", convert = TRUE)
  expect_identical(samp_annot(result)$Samples, c("001", "002", "100000"))
  expect_identical(samp_annot(result)$group, c("B", "", "A"))
  expect_identical(samp_annot(result)$score, c(1.5, NA, 3))
  expect_identical(samp_annot(result)$flag, c(FALSE, TRUE, TRUE))

  # The default leaves the columns as they are
  result <- add_samp_annot(musica, annot, "id")
  expect_identical(samp_annot(result)$score, c("1.5", "", "3"))
  expect_identical(samp_annot(result)$flag, c("FALSE", "TRUE", "TRUE"))
})

test_that("add_samp_annot() checks its inputs", {
  data(musica_sbs96)
  samples <- sample_names(musica_sbs96)
  annot <- data.frame(id = samples, group = "a")
  expect_error(add_samp_annot(musica_sbs96, as.matrix(annot), "id"),
               regexp = "'annotations' must be a data.frame")
  for (bad_column in list(c("id", "group"), NA_character_, 1)) {
    expect_error(add_samp_annot(musica_sbs96, annot, bad_column),
                 regexp = "'sample_column' must be a single column name")
  }
  expect_error(add_samp_annot(musica_sbs96, annot, "sample"),
               regexp = "'sample_column' \\('sample'\\) is not a column")
  expect_error(add_samp_annot(musica_sbs96, annot, "group"),
               regexp = "None of the sample names in column 'group'")
  duplicated_annot <- rbind(annot, annot[2, ])
  expect_error(add_samp_annot(musica_sbs96, duplicated_annot, "id"),
               regexp = paste0("more than once.*", samples[2]))
  samples_annot <- data.frame(id = samples, Samples = "a")
  expect_error(add_samp_annot(musica_sbs96, samples_annot, "id"),
               regexp = "'Samples' column")
  expect_error(add_samp_annot(musica_sbs96, annot["id"], "id"),
               regexp = "no columns other than 'sample_column'")
  for (bad_convert in list(NA, "yes", c(TRUE, FALSE))) {
    expect_error(add_samp_annot(musica_sbs96, annot, "id",
                                convert = bad_convert),
                 regexp = "'convert' must be TRUE or FALSE")
  }
  # A column that fails gives only the error, not the missing-sample warning
  list_annot <- data.frame(id = samples[-1])
  list_annot$group <- I(as.list(seq_len(nrow(list_annot))))
  expect_no_warning(
    expect_error(add_samp_annot(musica_sbs96, list_annot, "id"),
                 regexp = "must be a vector")
  )
})
