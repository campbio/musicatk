# Repair sample alignment in the shipped musica datasets (PR A of #208).
#
# The original scripts and sources for these objects are not in the
# repository, so this script repairs the existing objects instead of
# rebuilding them. No count or exposure value changes; only the order and
# labelling of samples do.
#
# For each musica object it:
# - keeps the current sample_names() order as the canonical sample order
#   (for musica, musica_annot and res_annot this matches
#   inst/extdata/sample_annotations.txt);
# - fills in the Samples annotation column where it is missing;
# - reorders count-table and exposure columns (and UMAP rows) to match;
# - drops count-table columns for samples with no variants, after checking
#   they are all zero;
# - drops unused levels from the variant sample factor.
#
# full_benchmark_example gets the same treatment for its ground truth, whose
# true exposures are transposed to signatures x samples.
#
# Run from the package root with the package source loaded:
#   Rscript data-raw/repair_sample_alignment.R
# No random numbers are used, so no seed is needed.

pkgload::load_all(".", quiet = TRUE)

# Repair one musica object; uses slots directly because this is a one-off
# migration of objects that do not yet satisfy the class invariant
repair_musica <- function(musica) {
  samples <- sample_names(musica)
  for (name in names(tables(musica))) {
    counts <- get_count_table(tables(musica)[[name]])
    missing <- setdiff(samples, colnames(counts))
    if (length(missing) > 0) {
      stop("Table ", name, " lacks samples: ", paste(missing, collapse = ", "))
    }
    extra <- setdiff(colnames(counts), samples)
    if (any(counts[, extra] != 0)) {
      stop("Table ", name, " has counts for samples with no variants.")
    }
  }
  if (nrow(samp_annot(musica)) == 0) {
    musica@sample_annotations <- data.frame(
      Samples = factor(samples, levels = samples)
    )
  }
  musica <- musicatk:::.subset_musica_samples(musica, samples)
  return(musica)
}

# Check that values are unchanged when matched by sample name
check_same_values <- function(before, after) {
  samples <- sample_names(after)
  for (name in names(tables(after))) {
    old <- get_count_table(tables(before)[[name]])
    new <- get_count_table(tables(after)[[name]])
    stopifnot(identical(colnames(new), samples))
    stopifnot(identical(old[, samples, drop = FALSE], new))
  }
  for (result in names(result_list(after))) {
    for (modality in names(get_result_list_entry(after, result)@modality)) {
      models <- get_modality(after, result, modality)
      for (model_id in names(models)) {
        old <- exposures(before, result, modality, model_id)
        new <- exposures(after, result, modality, model_id)
        stopifnot(identical(old[, colnames(new), drop = FALSE], new))
      }
    }
  }
  stopifnot(identical(
    as.character(variants(before)$sample),
    as.character(variants(after)$sample)
  ))
}

musica_datasets <- c(
  "musica", "musica_annot", "musica_sbs96", "musica_sbs96_tiny",
  "dbs_musica", "indel_musica", "res", "res_annot"
)
for (dataset in musica_datasets) {
  env <- new.env()
  utils::data(list = dataset, envir = env)
  before <- get(dataset, envir = env)
  after <- repair_musica(before)
  check_same_values(before, after)
  assign(dataset, after, envir = env)
  save(list = dataset, envir = env, compress = "gzip",
    file = file.path("data", paste0(dataset, ".rda"))
  )
  message(dataset, ": ", length(sample_names(after)), " samples")
}

# full_benchmark_example: transpose the true exposures and repair the truth
env <- new.env()
utils::data("full_benchmark_example", envir = env)
full_benchmark_example <- env$full_benchmark_example
truth <- ground_truth(full_benchmark_example)
model <- get_model(truth, "ground_truth", "SBS96", "ground_truth")
old_exposures <- exposures(model)
count_samples <- colnames(get_count_table(tables(truth)$SBS96))
true_signatures <- colnames(signatures(model))
if (identical(rownames(old_exposures), true_signatures)) {
  stop("full_benchmark_example is already signatures x samples.")
}
# The truth was stored samples x signatures in count-table column order;
# its row names were dropped, so check the shape and any names present
stopifnot(
  nrow(old_exposures) == length(count_samples),
  identical(colnames(old_exposures), true_signatures),
  is.null(rownames(old_exposures)) ||
    identical(rownames(old_exposures), count_samples)
)
new_exposures <- t(old_exposures)
colnames(new_exposures) <- count_samples
model@exposures <- new_exposures
truth@result_list[["ground_truth"]]@modality[["SBS96"]][["ground_truth"]] <-
  model
truth@sample_annotations <- data.frame(Samples = count_samples)
ground_truth(full_benchmark_example) <- truth
save(full_benchmark_example, compress = "bzip2",
  file = file.path("data", "full_benchmark_example.rda")
)
message("full_benchmark_example: true exposures now signatures x samples")
