# Spec: sample alignment invariant (PR A of #208)

Approved by the maintainer on 2026-11-04.

## Problem

Samples are paired by position, never by name, and nothing keeps the
positions in agreement. Annotation rows, count-table columns and exposure
columns can be in different orders, so positional consumers silently give
wrong statistics and mislabelled plots. The shipped `musica`,
`musica_annot` and `res_annot` objects are already out of sync, and
`dbs_musica` / `musica_sbs96_tiny` carry all-zero columns for stale sample
factor levels. The vignette assigns `sample_annotations.txt` by position
to an object whose samples are in a different order, mislabelling 6 of 7
samples.

Issues: #153, #155, #158, #159, #160, #172, #144, #131, #133, #134, #194.

## Invariant

`sample_names(musica)` is the source of truth: `samp_annot()$Samples` as
character, exact column match (no partial matching), falling back to the
sorted unique variant samples only for legacy objects without annotations.

- Every count table's columns equal `sample_names()` in the same order.
- Every model's exposure columns are an ordered subset of `sample_names()`
  (per-group models from the grid functions hold a subset).
- Every variant sample is in `sample_names()`.

Validity methods enforce this on `musica`, and check internal consistency
of `result_model` (signature names match exposure rows) and
`result_collection` (every entry is a named `result_model`). The `musica`
replacement methods for variants, tables, annotations, results and
exposures call `validObject()`, as does `add_result()`, which also puts
exposure columns in sample order by name.

Deviation found in review: `num_signatures` is not checked, because the
package itself builds models whose `num_signatures` is the full count
while they hold a subset of signatures (`compare_results.R`,
`benchmarking.R`); that belongs with PR K's model-storage cleanup.

## Changes

- `R/class_musica.R`: validity method; `sample_names()` exact match and
  character; `samp_annot<-` initialises from `sample_names()`; setters
  validate; subset functions also subset exposures and UMAP rows, and drop
  unused factor levels. Internal `.annot_for_samples()` returns an
  annotation column aligned by name to a vector of samples.
- `R/class_result_model.R`, `R/class_result_collection.R`: validity
  methods.
- `R/load_data.R`: `create_musica_from_counts()` fills in `Samples`; it
  stops when the matrix has no column names or duplicated ones.
- `R/standard_tables.R`, `R/table_utils.R`: every builder counts with
  `factor(sample, levels = sample_names(musica))`; IND83 names columns from
  `sample_names()`; `.create_count_table()` reorders columns to
  `sample_names()` or stops if the sample sets differ;
  `combine_count_tables()` aligns columns by name before binding.
- `R/differential_analysis.R`, `R/plot_heatmap.R`, `R/kmeans.R`: pair
  annotations with exposures by name; `plot_cluster(group = "annotation")`
  uses the `annotation` argument.
- `R/discovery_prediction.R`: `combine_predict_grid()` orders exposure
  columns by sample and sets `num_signatures` and `model_id`.
- `R/benchmarking.R`: `create_benchmark()` stores true exposures as
  signatures x samples, matched to the count table by sample name
  (accepting either orientation); all truth-exposure sites index by
  signature row and match samples by name.

## Data

`data-raw/repair_sample_alignment.R` repairs the shipped objects rather
than rebuilding them (their original sources are not in the repo). For
each `musica` object it keeps the current `sample_names()` order, fills in
`Samples` where missing, reorders count-table and exposure columns to
match, and drops count-table columns for samples with no variants after
asserting they are all zero. No count values change. It re-saves
`full_benchmark_example` with transposed true exposures.

## Docs

- Class docs describe the invariant.
- Vignette and examples assign annotations by matching on sample names.
- NEWS.md: 2.7.2 entry. DESCRIPTION: 2.7.1 -> 2.7.2.

## Tests

One failing-first regression test per issue, one commit per issue.
Results change: differential-analysis p-values, heatmap annotation bars,
cluster facets and the vignette's annotated plots. The Scientific
correctness section of the PR is for a person to fill in.

## Out of scope (follow-up issues)

- `updateObject()` method for objects saved with older versions (needs
  BiocGenerics in Imports, so an ADR).
- Shiny `server.R` annotation upload merges, re-sorts, then assigns by
  position; fixing it needs a new exported function.
