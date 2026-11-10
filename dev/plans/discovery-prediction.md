# Spec: discovery and prediction correctness (PR F of #208)

Decisions by the maintainer on 2026-11-09: keep zero-count samples with
zero exposures; `predict_exposure()` gets `seed = 1`; keep the grid's
reconstruction error and document it; devel only, no release port.

Changed during review (maintainer, 2026-11-09):

- The seed is applied with `withr::with_seed()`, with withr moved to
  Imports (ADR 0002), because BiocCheck warns on `set.seed()` in package
  code.
- Zero-count samples are left out of the exposures with a message, not
  kept with zero exposures: zero columns made `create_umap()`,
  `cluster_exposure()` and proportional plots fail (0/0). The grid and
  `compare_k_vals()` match samples by name instead.
- From `/code-review`: per-group fits use an object without stored
  results; a grid error names its group; `combine_predict_grid()` keeps
  the per-group models' motif order.

Issues: #140, #141, #142, #143, #145, #146, #147, #148, #149, #152, #162.
All reproduced on `devel` at 80350b9.

## Problem

- `generate_result_grid()` reads `Tumor_Type` whatever `annotation` says,
  so other columns store nothing, and every group writes the same
  `lda<k>` model (#140). A re-run fetches the old `lda2` after discovery
  stored the new fit as `lda2.1`, so the new fit is lost (#141).
  `auto_predict_grid()` has the same fetch-by-guessed-name pattern for its
  per-group models.
- `predict_exposure()` never compares signature motifs with the count
  table: permuted rows predict silently, DBS78 signatures on an SBS96
  table only warn about recycling (#142).
- Samples with no mutations are dropped silently, so exposures have fewer
  columns than the table; the grid and `compare_k_vals()` then crash
  (#143).
- `auto_predict_grid()` fails for `result_name` other than `"result"`
  (#145) and for a factor annotation column (#146);
  `auto_subset_sigs()` errors opaquely when no signature passes (#147).
- Model-id de-duplication strips two characters, so the 12th copy is
  `lda2..11`; the block is copied three times (#148).
- `compare_k_vals()` with an invalid `error_type` runs every fit, then
  returns NA (#149).
- LDA prediction has no seed; one sample gives 100 warnings; the grid
  error metric is undocumented (#152).
- `.cosine()` drops to a vector for one signature, so
  `compare_cosmic_v2()` on a one-signature model errors (#162).

## Changes

Shared helpers, in `R/discovery_prediction.R`:

- `.unique_model_id(musica, result_name, modality, model_id)` returns
  `model_id`, or `<model_id>.<n>` with the smallest `n` not in use, with
  the existing message. Replaces the three copied blocks in
  `discover_signatures()`, `predict_exposure()` and
  `combine_predict_grid()` (#148).
- `.pad_zero_exposures(exposures, samples)`: callers fit on samples with
  at least one mutation, then this adds all-zero exposure columns for the
  others, so exposure columns equal the table's columns in order. A
  message names the zero-count samples (#143).

`discover_signatures()`: stops when fewer than two samples have mutations
(was: fewer than two samples); zero-count samples get zero exposures
(#143).

`predict_exposure()`:

- Motifs matched by name: signatures are reordered to the count table's
  rows; it stops, naming the modality and a few differing motifs, when
  the motif sets differ or either has no motif names (#142).
- `algorithm = "decompTumor2Sig"` stops unless the table has the 96 SBS
  motifs, since its motif-name conversion only handles SBS96 (#142).
- New `seed = 1` argument (added after `make_copy`, so positional calls
  keep working). For `"lda"` it sets the seed and restores the caller's
  RNG state on exit (#152).
- Zero-count samples get zero exposures (#143).
- `lda_posterior()` skips the Dirichlet update of `theta` when there is
  only one sample, which can't be fit and gave 100 warnings (#152).

`generate_result_grid()`:

- Uses the `annotation` column; stops, listing available columns, when it
  doesn't exist. Groups are the column's non-missing values as
  character, subset with `subset_musica_by_annotation()`; samples with a
  missing value are left out with a message (#140).
- Each model id is chosen up front with `.unique_model_id()` against the
  object being filled: `<algorithm><k>`, or `<algorithm><k>_<group>` with
  an annotation. Discovery is called with that id and the model fetched
  by it, so a re-run stores the new fit as `lda2.1` with a message
  (#140, #141).
- Reconstruction error: same formula, now matching samples by name, and
  documented in `@details`: for each sample, the mean squared difference
  between observed and reconstructed counts over motifs; the stored value
  is the mean of the squares of those per-sample values. Noted as not
  comparable with `compare_k_vals()` (#152).
- `verbose` prints with `message()` instead of `cat()` (lines touched).

`auto_predict_grid()` and internal `auto_subset_sigs()`:

- The first, all-signature pass computes exposures without storing a
  model, through an internal `.predict_exposures(counts_table, signature,
  algorithm, seed, verbose)` that `predict_exposure()` also uses (#145).
- Group values are converted to character and missing values skipped
  (#146). Per-group ids are chosen up front with `.unique_model_id()`,
  passed to `auto_subset_sigs()`, fetched by that id, and the actual ids
  are passed to `combine_predict_grid()`.
- `auto_subset_sigs()` stops with a message naming the three thresholds
  (and the group, if any) when no signature passes (#147).
- New `seed = 1` argument, passed to both passes (#152).

`R/k_val_assit.R`:

- `compare_k_vals()` validates `error_type` before any fit (same accepted
  spellings) and stops otherwise; `.get_reconstruction_error()` stops
  instead of messaging (#149).
- `.resample_count_table()` keeps zero-count samples at zero instead of
  passing NaN probabilities to `rmultinom()`; the proportional error
  treats a zero-count sample's proportions as 0, so it adds 0 instead of
  NaN (#143).

`R/misc_functions.R`: `.cosine()` uses `drop = FALSE` (#162).

Not changed here (other groups): the second `.get_reconstruction_error()`
in `benchmarking.R` (BEN-12, PR I); raw `@` model insertion and other
copied validation (PR K); mutate-in-place API (#205).

## Tests

New `tests/testthat/test-discovery-prediction.R`, one regression test per
issue, each failing on its assertion on `devel` first. Small fixtures:
`musica_sbs96`/`musica_annot` (12 samples), `nstart = 1`, `k = 2`, a few
COSMIC signatures, `reps = 1`.

- #140 grid with `annotation = "Tumor_Subtypes"` stores `lda2_Lung`,
  `lda2_Skin`, each with only its group's samples; unknown column errors.
- #141 grid re-run keeps the old `lda2` and stores the new fit as
  `lda2.1`, which differs from it.
- #142 permuted signature rows give the same exposures as unpermuted;
  DBS78 signatures on SBS96 error; decompTumor2Sig on DBS78 errors.
- #143 a zero-count sample gives a zero exposure column in discovery
  (LDA, NMF), both prediction algorithms, the grid, and `compare_k_vals()`
  runs; fewer than two samples with mutations errors.
- #145 `auto_predict_grid(result_name = "pred")` stores its model there.
- #146 factor annotation works and gives the same result as character.
- #147 impossible thresholds error with the new message.
- #148 12 repeats give `.1` … `.11`.
- #149 invalid `error_type` errors before any fit.
- #152 two `predict_exposure()` calls give identical exposures and leave
  `.Random.seed` unchanged; one sample gives no warning.
- #162 `.cosine()` with one column returns a 1-row matrix and
  `compare_cosmic_v2()` works on a one-signature model.

## Docs, NEWS, version

- Roxygen: `seed` for `predict_exposure()` and `auto_predict_grid()`;
  motif matching and zero-count behaviour in `predict_exposure()` and
  `discover_signatures()`; `annotation`, model ids and the metric in
  `generate_result_grid()`. `make docs`.
- Vignette "Note on reproducibility": say `discover_signatures()` and
  `predict_exposure()` take `seed`, and fix the `set_seed` typo. Render
  with `make article FILTER=musicatk`.
- NEWS.md: 2.7.12 entries, results-changing ones first. DESCRIPTION
  2.7.11 → 2.7.12.

## Results change

Yes. Grid results with an annotation (were missing or overwritten),
grid re-runs, prediction with signatures whose motif order differs from
the table, samples with no mutations (now present with zero exposures),
and LDA prediction, now seeded by default. Shipped `data/*.rda` objects
are unaffected: none store predictions or grid results. The PR's
Scientific correctness section needs a person.

## Checks

`make test`, `make coverage` (must not drop), `make lint` on changed
files, `make check-full` and `make bioccheck` before the PR.
