# Build data/full_benchmark_example.rda (PR I of #208).
#
# The ground truth is the synthetic breast cancer data set: 214 samples
# with 8 COSMIC v2 signatures (synthetic_breast_counts and
# synthetic_breast_true_exposures). One prediction is benchmarked against
# it: the NMF k = 8 result in example_predicted_sigs and
# example_predicted_exp, as method "example_k8".
#
# The summaries are data.frames since musicatk 2.7.18 (ADR 0004), so the
# object is rebuilt rather than repaired.
#
# Run from the package root with the package source loaded:
#   Rscript data-raw/make_full_benchmark_example.R
# No random numbers are used, so no seed is needed.

pkgload::load_all(".", quiet = TRUE)

env <- new.env()
utils::data(
  list = c(
    "synthetic_breast_counts", "synthetic_breast_true_exposures",
    "example_predicted_sigs", "example_predicted_exp", "cosmic_v2_sigs"
  ),
  envir = env
)

# ground truth
true_sigs <- c(
  "SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26"
)
true_signatures <- signatures(env$cosmic_v2_sigs)[, true_sigs]
full_benchmark_example <- create_benchmark(
  true_signatures, env$synthetic_breast_true_exposures,
  env$synthetic_breast_counts, "SBS96"
)

# the prediction to benchmark
prediction <- create_musica_from_counts(env$synthetic_breast_counts, "SBS96")
add_result(
  env$example_predicted_sigs, env$example_predicted_exp, prediction,
  "result", "SBS96", "example_k8"
)

full_benchmark_example <- benchmark(
  full_benchmark_example, prediction, "example_k8", "SBS96",
  method_id = "example_k8", plot = FALSE, make_copy = TRUE
)

save(full_benchmark_example, compress = "bzip2",
  file = file.path("data", "full_benchmark_example.rda")
)
message("Saved data/full_benchmark_example.rda")
