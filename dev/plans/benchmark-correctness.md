# Spec: Benchmarking correctness (PR I of #208)

Changes results: yes. `benchmark()` summaries change for any prediction
with duplicates or composites. The summary tables become data.frames
(ADR 0004, accepted 2026-11-13), and `full_benchmark_example` is rebuilt. Devel only; no
release port proposed (a large change to a niche workflow, and most of
the bugs are crashes rather than silently wrong numbers).

Fixes #188, #189, #190, #191, #192, #193, #195, #196, #197, #198, #199,
#200, #201 and #94. It also fixes `benchmark()`'s own copy of the
`method_id` de-duplication bug from #148.

## Problem (reproduced on devel e6b0c07 unless noted)

- **#188.** A predicted signature that matches two true signatures, each
  of which also has another match, is merged into both duplicate groups.
  In a 5-signature simulation, total exposure grew from 63,759 (Initial)
  to 72,118 (Intermediate), and "Num Direct Matched" was 0 although P0
  matched SBS1 directly.
- **#189.** A composite of three true signatures stops with "undefined
  columns selected": the per-sample NNLS results go into a data.frame
  with only `factor1` and `factor2`.
- **#190.** Composite correction reads `musica@count_tables$SBS96` and
  builds 96-row matrices. A DBS composite stops with "no applicable
  method for `@` applied to NULL".
- **#191.** Composite correction replaces a composite with copies of the
  true signatures it matched (each NNLS-scaled component is renormalised
  back to the true signature), so the Final cosine is 1.000 by
  construction. A 3-signature prediction with a SBS3+SBS8 composite gets
  Mean CS 0.997 and Final RE 3.22. The correct 4-signature prediction
  gets 0.987 and 3.91.
- **#192.** Already fixed on devel by 8321ad2 (PR A, #210): `drop =
  FALSE` when one signature is left unchanged. The duplicate path has a
  test; the composite path has none.
- **#193.** `benchmark_get_comparison(..., "final")` passes validation
  and then stops with "object 'comparison' not found": dispatch uses
  exact `==` on the three canonical names.
- **#195.** `create_benchmark(..., "SBS")` stores the truth under
  modality `"SBS"`, so `benchmark(..., "SBS96")` stops with "SBS96 is not
  a valid modality". `create_benchmark(..., "DBS78")` is rejected by
  `create_musica_from_counts()` ("Only SBS and DBS classes are
  supported"), and `"DBS"` stores modality `"DBS"`, which never matches a
  prediction made from the `DBS78` table. So DBS benchmarks can't be run
  at all.
- **#196.** `benchmark(plot = FALSE)` draws three pages, because it calls
  `compare_results()`, which always plots.
- **#197.** (From the code.) `predict_and_benchmark()` benchmarks
  `paste0(algorithm, k)`. When that `model_id` already exists,
  `discover_signatures()` stores the new fit as `<id>.1`, so the old model
  is benchmarked.
- **#198.** (From the code.) "Num Found" counts rows of the Final
  comparison, so a true signature matched twice counts twice. "Num
  Missing" comes from the Initial comparison.
- **#199.** `.get_reconstruction_error()` is defined in
  `benchmarking.R` (2 arguments) and in `k_val_assit.R` (3 arguments).
  `.generate_summary()` calls it with 3 arguments and works only because
  `k_val_assit.R` is collated later.
- **#200.** A prediction with no match above `threshold` stops in
  `sig_compare()` with "No matches found, try lowering threshold."
- **#201.** `method_view_summary` and `sig_view_summary` are character
  matrices, so Mean CS etc. are strings. `.generate_summary()` uses
  `exists("duplicated_signature_components")` without `inherits =
  FALSE`, so a global variable of that name changes the counts.
  `.signature_view_summary()` calls `min(NULL)`/`max(NULL)` (warnings),
  and its `is.infinite` clean-up never fires.
- **#94.** Nine `class(x)[1] != "..."` checks (BiocCheck WARNING); a
  subclass of `full_benchmark` is rejected.
- **#148 (benchmark copy).** `benchmark()` de-duplicates `method_id` by
  stripping two characters, so the 11th repeat becomes `x..10`.

## Approach

### Matching (#196, #200)

- New internal `.signature_matches(sig1, sig2, threshold, metric)` in
  `R/compare_results.R` returns the comparison data.frame, possibly with
  0 rows. `sig_compare()` calls it and keeps its "No matches found"
  error, so `compare_results()` is unchanged.
- `benchmark()` compares with `.signature_matches()` followed by
  `.benchmark_comp_adj()`, with no plotting. A prediction with no match
  gives an empty comparison and a summary with Found 0, Missing = number
  of true signatures, Spurious = number of predicted signatures, and NA
  cosines. With `plot = TRUE`, the comparison and exposure plots are
  skipped for an empty stage, with a message. `benchmark_plot_comparison()`
  and `benchmark_plot_exposures()` stop on an empty stage, naming it.

### Duplicates (#188)

- New `.duplicate_groups(comparison)`: each predicted signature joins
  the group of the true signature it matches best (highest cosine; ties go
  to the first row, which is the comparison's order). A group with two or
  more predicted signatures is a duplicate. `.correct_duplicates()`,
  `.generate_summary()` and `benchmark_plot_duplicate_exposures()` all
  use it, so they agree.
- The existing `!grepl("like", x_sig_name)` filter is dropped. It
  silently excluded any user signature whose name contains "like" (e.g.
  "APOBEC-like"), and Initial predictions never contain the `_like`
  pieces it was meant for.
- The merge itself is unchanged: the exposure-weighted mean of the
  signatures, and the sum of their exposures. The ggplot objects built
  and thrown away inside `.correct_duplicates()` are removed with the
  rewrite (the rest of #202 stays in group K).

### Composites (#189, #190, #191)

For each composite C (a predicted signature matching two or more true
signatures T1..Tn in the Intermediate comparison):

1. `a <- lsei::nnls(T, C)$x`, and the fitted parts `F_k = a_k T_k`.
2. Split C motif by motif in proportion to the fitted parts:
   `D_k[m] = C[m] * F_k[m] / sum_j F_j[m]`. Where every `F_j[m]` is 0,
   `C[m]` is split by `a_k / sum(a)`. The pieces add up to C exactly, so
   the noise in C stays in the pieces.
3. Components with `a_k = 0` get no piece.
4. Piece k's signature is `D_k / sum(D_k)`, and its exposure is C's
   exposure times `sum(D_k) / sum(C)` (the maintainer's choice of a fixed
   fraction, 2026-11-13). The reconstruction `sigs %*% exposures` is
   unchanged by the split, so Final RE gets no credit from the truth.
5. Pieces are named `<true>_like`, made unique if two composites share a
   true signature. One helper, `.composite_pieces(comparison, ...)`,
   gives the composite, true signature and piece name for both
   `.correct_composites()` and `benchmark_plot_composite_exposures()`.

This works for any number of components and any modality; the count
table and the `musica` argument are no longer used. The discarded ggplot
objects in `.correct_composites()` go with the rewrite.

### Summary (#191, #198, #199, #201)

- **Final stage for the headline numbers** (maintainer's decision,
  2026-11-13): for each true signature, its best cosine in the Final
  comparison. Num Found = true signatures with a match there; Num Missing
  = the rest (they add up to the number of true signatures). Mean, Min
  and Max CS are over the found true signatures, one value each, and NA
  when none is found. Sigs Found lists them once each, in true-signature
  order.
- Num Spurious, Duplicates, Composites, Dup/Comp and Direct Matched still
  describe the Initial prediction, now counted with `.duplicate_groups()`.
  The documentation says which stage each row comes from.
- `sig_view_summary`: per run, a signature is found if it is in the
  Final comparison, with its best cosine for that run; Mean/Min/Max over
  runs are NA when it was never found, with no `min(NULL)` warnings.
- **ADR 0004, data.frame summaries** (maintainer's decision,
  2026-11-13): `method_view_summary` is a data.frame with one row per
  run (row names = method IDs), and `sig_view_summary` one row per true
  signature. Column names are the current row labels ("Num Found", "Mean
  CS", ...). Counts are integer, RE and CS numeric, "Sigs Found" and
  "Methods Found" character. `single_summary` (in `single_benchmark`) is
  the one-row data.frame. The slot classes, prototypes and replacement
  methods change from `matrix` to `data.frame`. `benchmark()` stops with
  "made with an older version of musicatk; recreate it with
  create_benchmark()" when given a `full_benchmark` with matrix summaries.
- `exists()` is replaced by an ordinary variable; dead `is.infinite`
  removed.
- **#199.** The 2-argument `.get_reconstruction_error()` in
  `benchmarking.R` is deleted. `re_method` goes through the
  `k_val_assit.R` version, which also leaves out zero-count samples
  (#143).

### Smaller fixes

- **#193.** `.match_benchmark_stage(prediction)` maps the accepted
  aliases, case-insensitively, to "Initial", "Intermediate" or "Final",
  and stops otherwise. It is used by `benchmark_get_prediction()`,
  `benchmark_get_comparison()` and the plot functions.
- **#195.** `create_benchmark()` turns `variant_class` into the table
  name: the SBS aliases into "SBS96", and the DBS aliases plus "DBS78"
  into "DBS78". Anything else stops, naming the two. The truth is stored
  under that modality. `create_musica_from_counts()` accepts "DBS78" as a
  DBS alias (a one-word change in `R/load_data.R`). `benchmark()` stops
  early if `prediction_modality` isn't the ground truth's modality,
  naming it. The docs say "SBS96" or "DBS78" (IND83 was listed but has
  never worked).
- **#197.** `predict_and_benchmark()` picks the new `model_id` with
  `.unique_model_id()` before discovery, passes it to
  `discover_signatures(model_id =)`, and benchmarks that model. Its
  `@return` text (which mentions a nonexistent `make_copy`) is corrected.
- **#148 copy.** A small `.unique_id(id, existing, label)` holds the
  `.<n>` logic. `.unique_model_id()` and `benchmark()`'s `method_id` both
  use it.
- **#94.** `methods::is()` for all nine checks.

Out of scope (group K, #202): the validation block copied into each
export, `print()` of the summary, `benchmark_plot_comparison()` not
returning its plot, and the `assign(parent.frame())` pattern.

## Tests

In `tests/testthat/test-benchmarking.R`, on small simulations (20
samples, COSMIC v2 or v3 DBS signatures, fixed seed). One regression
test per issue, each failing on its assertion first:

- #188: total exposure is the same in Initial and Intermediate, and the
  shared signature is in one group only.
- #189: a SBS3+SBS5+SBS8 composite gives three `_like` pieces.
- #190: a DBS78 composite benchmarks (`create_benchmark(..., "DBS78")`).
- #191: with a noisy composite, the pieces' Final cosines are below 1.
  `sigs %*% exposures` is equal in Intermediate and Final, and Mean CS is
  the mean of the best Final cosine per true signature.
- #192: a composite that leaves one signature unchanged keeps its name.
  This passes on devel (fixed in 8321ad2); it guards the composite path.
- #193: "final", "Init" and "inter" work in `benchmark_get_comparison()`.
- #195: "SBS" stores "SBS96", and benchmarking with "SBS96" works. A
  mismatched `prediction_modality` errors, naming the modality.
- #196: `plot = FALSE` with `.plot_compare_result_signatures` mocked to
  fail runs without error.
- #197: with `lda2` already present, `predict_and_benchmark(k_min = 2,
  k_max = 2, algorithm = "lda", nstart = 1)` benchmarks `lda2.1`.
- #198: in a run where one true signature is matched twice in Final,
  Found + Missing equals the number of true signatures.
- #199: Initial RE matches the prop and raw formulas. This can't fail
  first: behaviour is already right and the bug is the Collate
  dependence.
- #200: a prediction with no match gives Found 0, Missing 3, NA cosines,
  and no error.
- #201: summary columns are integer/numeric/character; a global
  `duplicated_signature_components` doesn't change the counts; no
  warnings from `sig_view_summary`.
- #94: a subclass of `full_benchmark` is accepted, and a non-benchmark
  still errors.
- #148 copy: the 11th repeat of a `method_id` is `x.10`.
- Existing tests are kept; any that check matrix summaries are updated
  to the data.frame layout. No test is weakened.

## Data, docs and other

- `data-raw/make_full_benchmark_example.R` (new, seed fixed) rebuilds
  `data/full_benchmark_example.rda` from the bundled synthetic breast
  data and `example_predicted_sigs`/`_exp` (method ID `example_k8`). Its
  numbers don't change: it has no duplicates or composites. Only the
  summary storage changes. The `R/test_data.R` entry points to the script.
- `dev/adr/0004-benchmark-summaries-as-data-frames.md`, status Proposed
  until the maintainer accepts it. Index row in `dev/adr/README.md`.
- Roxygen: `create_benchmark` `variant_class`, `benchmark`
  `prediction_modality` and a description of the summary rows and their
  stages, `predict_and_benchmark` `@return`, the class slot docs. Then
  `make docs`.
- The article `vignettes/articles/benchmarking_tutorial.Rmd` prints both
  summaries, which are now transposed. Update its text and preview it
  with `make article FILTER=benchmarking_tutorial`. The vignette doesn't
  use benchmarking.
- NEWS under 2.7.18; version 2.7.17 -> 2.7.18.
- Lint: new and edited code is lint-clean. `R/benchmarking.R` as a whole
  still has about 960 lints, so the changed-files lint CI job will fail
  on this PR. #257 plans that file's clean-up after these fixes.
- Checks: `make test`, `make coverage` (must not drop), `make lint`;
  `make check-full` and `make bioccheck` before the PR. The #94 BiocCheck
  WARNING should be gone.
