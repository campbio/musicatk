# Agent Findings Log

Findings produced by automated audits and setup verification. Findings are
recorded here first, then become GitHub issues. **Nothing here is fixed in the
run that discovered it** — see `dev/AUDIT.md`.

---

## 2026-09-18 — AI scaffolding setup (playbook v2) verification

Discovered while running the Step 9 verification pass. None were fixed, per the
playbook scope guard.

### 1. 33 exported topics missing from the pkgdown reference index

`pkgdown::check_pkgdown()` fails. Every missing topic belongs to the benchmarking
framework introduced in commit `1dd3973`, which did not update `_pkgdown.yml`.
None of these functions currently appear on the package website.

```
adjustment_threshold, benchmark, benchmark_get_comparison, benchmark_get_entry,
benchmark_get_prediction, benchmark_plot_comparison,
benchmark_plot_composite_exposures, benchmark_plot_duplicate_exposures,
benchmark_plot_exposures, benchmark_plot_signatures, create_benchmark,
description, example_predicted_exp, example_predicted_sigs, final_comparison,
final_pred, full_benchmark-class, full_benchmark_example, ground_truth,
indv_benchmarks, initial_comparison, initial_pred, intermediate_comparison,
intermediate_pred, method_id, method_view_summary, predict_and_benchmark,
sig_view_summary, single_benchmark-class, single_summary,
synthetic_breast_counts, synthetic_breast_true_exposures, threshold
```

**Judgment required:** these are not all public API. Accessor generics such as
`description`, `threshold`, `adjustment_threshold`, `method_id` and
`ground_truth` may belong under `@keywords internal` instead of the index. A
maintainer should split the list before anyone edits `_pkgdown.yml`.

**Action:** open an issue. The `pkgdown-check` CI job is `continue-on-error`
until it is closed, then flip it to required.

### 2. Tarball is at the Bioconductor size limit

`R CMD build --no-build-vignettes` produces a **5.0 MB** tarball. Bioconductor
guidance is under 5 MB, and that figure excludes built vignettes.

Avoidable payload, all shipping because `vignettes/articles/` is not in
`.Rbuildignore`:

| Path | Size | Referenced by |
|---|---|---|
| `vignettes/ui_screenshots/` | 4.7 MB | `vignettes/articles/tutorial_tcga_ui.Rmd` only (pkgdown-only) |
| `vignettes/figures/` | 496 KB | nothing — stale knitr chunk output |
| `vignettes/articles/MANIFEST.txt` | 92 KB | nothing — rsconnect deploy artifact |

The only real vignette, `vignettes/musicatk.Rmd`, references no external images.

**Action:** open an issue. Adding `^vignettes/articles$` to `.Rbuildignore` and
deleting `vignettes/figures/` is a small, separate PR.

### 3. Duplicate pkgdown output committed under `docs/`

`docs/articles/articles/` (~15 MB) is stale output from a pre-2.0 pkgdown, which
mirrored source paths. pkgdown 2.x flattens `vignettes/articles/` to
`docs/articles/`. Because `docs/` is committed and pkgdown never deletes files it
did not generate, rebuilding does not clear it — this is why the attempt recorded
on the `temp2` branch ("Attempt to fix double articles by redoing build") was
reverted. Requires `clean_site()` then rebuild, or the gh-pages migration.

Also present: ~48 committed macOS `" 2"` duplicate files across `docs/` and
`vignettes/`, and one orphaned `docs/articles/musicatk_help.html` with no source.

**Action:** subsumed by the deferred pkgdown migration (playbook Step 5).

### 4. `bioc` git remote points at the wrong package

`git@git.bioconductor.org:packages/celda.git` rather than `packages/musicatk.git`.
A release push would go to the wrong repository. See `dev/RELEASE.md`.

**Action:** `git remote set-url bioc git@git.bioconductor.org:packages/musicatk.git`
— per-clone, so each maintainer fixes their own.

### 5. 2,922 lints across a never-linted codebase

`lintr` has never been run on this package. `make lint` reports **2,922 lints**
— 1,583 in `R/`, 1,336 in `inst/`, 3 in vignettes.

| Linter | Count |
|---|---|
| `trailing_whitespace_linter` | 713 |
| `object_usage_linter` | 541 |
| `line_length_linter` | 421 |
| `indentation_linter` | 384 |
| `return_linter` | 227 |
| `brace_linter` | 210 |
| `paren_body_linter` | 82 |
| `object_name_linter` | 80 |
| `commas_linter` | 72 |
| `infix_spaces_linter` | 62 |
| others | ~130 |

**1,536 of these (53%) are mechanically fixable by `make style`** — whitespace,
indentation, braces, commas, spacing. Those should land as one formatting-only
PR that touches nothing else, so it can be reviewed by confirming the diff is
whitespace.

The remainder need judgment. `object_usage_linter` (541) is the interesting
group — it flags undefined globals and unused variables, which can indicate real
bugs rather than style problems. `line_length_linter` (421) overlaps with what
BiocCheck enforces independently, so that one has to be fixed regardless.

**Action:** open two issues — one for the mechanical `make style` pass, one for
triaging `object_usage_linter`. The `lint` CI job is `continue-on-error` until
the count is down; flip it to required then.

**Note on config:** `.lintr` is DCF format and does not accept `#` comments — a
comment line makes lintr abort with "Malformed config file". Also,
`cyclocomp_linter` is not in lintr >= 3.4 defaults, so it must not be listed in
`linters_with_defaults()`.
