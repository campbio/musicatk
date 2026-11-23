# Spec: `overwrite` for model IDs; Shiny duplicate IDs and fit errors (#282, #283, part of #122)

The maintainer approved this design on 2026-11-23, in two sections, in
chat. A clash on `model_id` now stops with an error by default, and
`overwrite = TRUE` replaces the model. `add_result()` gets the same
argument. There is no exported function for deleting a model, and the
Shiny app has no replace option: it asks for another ID.

## Problem

- `discover_signatures()` and `predict_exposure()` rename a clashing
  `model_id` to `x.1`, `x.2`, and so on, with a message. A user who wants
  to replace a model with that ID can't, and there is no exported way to
  delete one.
- `add_result()` does the opposite: it silently replaces an existing
  model with that ID (#122).
- The Shiny app tries to ask "This model ID already exists" before
  Discover and Predict, but checks `%in% names(vals$res_names)`.
  `res_names` is an unnamed vector, so the prompt never shows. The
  package renames the model to `x.1`, while the app records `x` again.
  The new model can't be reached in the app (#282).
- In the app, `disc_sigs()` and `get_predict()` call the package with no
  `tryCatch()`, so any error ends the session. Picking DBS signatures for
  an SBS96 table is one example (#283). The success notification and the
  `res_names` update run even when nothing was stored.

## Package (`R/`)

### Internal helpers (`R/discovery_prediction.R`, next to `.unique_model_id()`)

- `.check_model_id(musica, result_name, modality, model_id, overwrite)`
  - Stops when `model_id` is already a model in `result_name`/`modality`
    and `overwrite` is FALSE. Message: "A model named '<id>' already
    exists in result '<result_name>', modality '<modality>'. Use another
    'model_id', or 'overwrite = TRUE' to replace it."
  - Also stops unless `overwrite` is a single TRUE or FALSE.
  - Returns nothing.
- `.store_model(musica, result_name, modality, model_id, model)`, in
  `R/class_result_collection.R`, since it uses `@`
  - Creates the `result_collection` entry and the modality list when
    they're missing.
  - Puts `model` at `model_id`, replacing any model there.
  - Calls `methods::validObject()`, and returns the `musica`.
  - Replaces the copied insertion blocks in `discover_signatures()`,
    `predict_exposure()` and `add_result()`. The grid functions keep
    theirs; that's progress on #165/#151, not their fix.

### Exported functions

- `discover_signatures()`, `predict_exposure()` and `add_result()` gain
  `overwrite = FALSE` as their last argument, so positional calls keep
  working.
- `discover_signatures()` and `predict_exposure()` call
  `.check_model_id()` right after the default `model_id` is filled in, in
  place of `.unique_model_id()`. That is before fitting, so a clash
  fails fast. In `predict_exposure()` that means moving the default-ID
  code above `.predict_exposures()`.
- `add_result()` calls `.check_model_id()` before building the model.
- When `overwrite = TRUE` replaces a model, it is replaced in place: its
  position in the modality list is kept.
- Roxygen:
  - `@param overwrite` for all three functions.
  - In `discover_signatures()`, the sentence "If `model_id` is already
    used, `.1`, `.2`, and so on is added to it." becomes the new rule.
  - `predict_exposure()` and `add_result()` get the same sentence.
- `auto_subset_sigs()` is internal, and only `auto_predict_grid()` calls
  it. It fills in the default `model_id` and makes it unique with
  `.unique_model_id()` before calling `predict_exposure()`, so
  `auto_predict_grid()` keeps adding `.1` with or without
  `sample_annotation`. (Changed during implementation: the spec first
  said it was exported and gave it `overwrite`; review found that
  `auto_predict_grid()` without an annotation then errored on a rerun.)
- Unchanged:
  - `generate_result_grid()`, `auto_predict_grid()`,
    `combine_predict_grid()` and `predict_and_benchmark()`. They make their IDs
    unique with `.unique_model_id()` before calling the functions above,
    so they never clash. Their docs, and the tests at
    `test-discovery-prediction.R:178-196` and `test-benchmarking.R:562`,
    still describe them correctly.
  - `compare_k_vals()`. Each fit starts from the user's unmodified
    object, so its fixed IDs `"actual"` and `"boot"` don't clash.

**Breaking:**
- Rerunning `discover_signatures()` or `predict_exposure()` with an ID
  already in that result and modality now errors. That includes the
  default ID, for example `"lda3"`.
- `add_result()` errors instead of replacing.

The NEWS entry says so in bold.

### Tests (`tests/testthat/`)

For each of `discover_signatures()`, `predict_exposure()` and
`add_result()`:
- A clash errors with the message above, and the stored models are
  unchanged.
- `overwrite = TRUE` replaces the model: same IDs, in the same order, and
  the new model's exposures.
- A new ID is still added.

Also:
- `discover_signatures()` checks the ID before fitting. With a clashing
  ID and `num_signatures = 1000`, it stops with the clash message, not a
  fitting error.
- A non-logical `overwrite` errors.

Use small fixtures: `musica_sbs96_tiny`, `nstart = 1`, `lda`. The
discovery-prediction test file is slow, so new tests go in a new file,
`test-model-id.R`.

## Shiny app (`inst/shiny/server.R`)

### Duplicate IDs (#282)

- The Discover and Predict guards also fail when
  `paste0(result, "/", table, "/", model) %in% vals$res_names`. Their
  validation output shows: "A model with this ID already exists in this
  result list entry. Choose another model ID."
- Remove the overwrite modals and the `confirm_result_overwrite` and
  `confirm_predict_overwrite` observers.
- The app calls the package with the default `overwrite = FALSE`, so the
  package error backs up the app's check.

### Fit errors (#283)

- `disc_sigs()` and `get_predict()` each run their package call inside
  `tryCatch()`.
- On error: `shinyalert::shinyalert(title = "Error", text = cond$message)`.
  `vals` is unchanged.
- On success: update `vals$musica`, `vals$res_names` and
  `vals$res_names_full`, then show the notification. The notifications
  move from the observers into these functions.

### Left out

Filtering the Predict "Result to Predict" choices to the chosen table's
signature type. With the `tryCatch()`, a mismatched set shows the
package's motif-mismatch message. A new enhancement issue tracks the
filter.

### Tests (`tests/testthat/test-shiny-server.R`)

- Discover with an ID already in `vals$res_names`: the duplicate
  message shows, and the stored models are unchanged. The same for
  Predict.
- Predict with `vals$p_res <- cosmic_v3_dbs_sigs` on the SBS96 table:
  - one alert, with the motif-mismatch message;
  - the session stays open;
  - `vals$res_names` is unchanged;
  - no success notification.
- A small successful Discover (`musica_sbs96_tiny`, `lda`, 2 signatures,
  `nstart = 1`):
  - `vals$res_names` gains the ID;
  - the model is in `vals$musica`;
  - the notification shows.

Each new test must fail on devel, except the successful Discover, which
guards the moved notification and record code and passes on devel too.

## Docs, NEWS, version, checks

- `make docs`; keep the roxygen version churn out (see the
  `RoxygenNote` memory).
- NEWS: version 2.7.28, with a breaking-change entry for the package and
  an entry for the app.
- Run `make test`, `make lint` on the changed files, `make check-full`,
  `make bioccheck` and `make coverage`.
- Check the vignette and `vignettes/articles/` for reruns that reuse a
  model ID, which would now error. Fix any by giving each run its own
  `model_id`.
- After merging:
  - close #282 and #283 by hand;
  - comment on #122 that the silent overwrite is fixed. Its
    deparse/assign, dimension-check and doc points stay open.
  - comment on #165 and #151 about `.store_model()`;
  - file the Predict-filter enhancement issue.
