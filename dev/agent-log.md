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

**Resolved (2026-09-25).** The maintainer decided all 33 are public API and
none go under `@keywords internal`. They were added to `_pkgdown.yml` in new
"Benchmarking" and "Accessing benchmark data" sections, with the classes and
datasets under "Data objects and classes". `check_pkgdown()` passes, and the
`pkgdown-check` CI job is now blocking.

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
triaging `object_usage_linter`.

**Update (2026-09-25).** The `lint` CI job now lints only the R files a PR
changes, and it is blocking. New and edited files must be clean, and the
backlog shrinks one file at a time. `.lintr` also disables `return_linter`
(explicit `return()` is common Bioconductor style) and `object_usage_linter`
(false positives from data.table and ggplot column names). That brings the
full backlog (`make lint`) to 2,154. With `object_usage_linter` off, nothing
checks for undefined globals or unused variables any more, so the triage issue
is the only remaining guard. Re-enable the linter once the false positives are
silenced with `utils::globalVariables()` or `# nolint`.

**Note on config:** `.lintr` is DCF format and does not accept `#` comments — a
comment line makes lintr abort with "Malformed config file". Also,
`cyclocomp_linter` is not in lintr >= 3.4 defaults, so it must not be listed in
`linters_with_defaults()`.

---

## 2026-11-04 — Whole-package review (R/, tests/, DESCRIPTION)

Read-only review of `R/` (12,024 lines), `tests/`, DESCRIPTION and NAMESPACE,
split into six areas, each reviewed by a separate agent against the same
checklist (bugs, code quality, efficiency, coverage gaps). `inst/shiny/` is out
of scope and gets its own review. **Nothing was fixed in this run** (see
`dev/AUDIT.md`); findings become GitHub issues after maintainer triage.

**Evidence levels.** *Verified* = reproduced independently by the lead after
the area review. *Confirmed* = reproduced by the area reviewer with a script
(scripts kept outside the repo). *Plausible* = traced in code, not run. Line
numbers are against `devel` at 2b8cb17.

**Baseline (devel 2b8cb17, R 4.5, ggplot2 4.0.3).** `make test`: 41 pass, 0
fail, 2 warnings (both #95). `make coverage`: **9.62%**; 13 of 28 files at 0%
(annotate_variants, benchmarking, class_count_table, class_full_benchmark,
class_single_benchmark, compare_results, k_val_assit, misc_functions,
musicatk, plot_sample_counts, plot_sample_reconstruction_error,
plot_signatures, standard_tables); table_utils 2.1%, discovery_prediction
2.3%, load_data 7.7%. `make lint`: 2,154 lints.

**Results change.** Many fixes below change counts, exposures, p-values or
plots. Every PR built from them needs the "Scientific correctness" section
filled in by a person, and the shipped `data/*.rda` objects need rebuilding
(see CLS-1, TAB-10).

### Cross-cutting themes

1. **Samples are aligned by position, not by name.** Annotation rows,
   count-table columns and exposure columns can be in different orders, and
   nothing enforces agreement. All seven shipped musica datasets are already
   out of sync (`res_annot`: annotations start TCGA-94, exposures start
   TCGA-56; verified). Positional consumers then give wrong statistics and
   mislabelled plots with no error. Sites: CLS-1, PLT-1, PLT-5, DIS-4, DIS-5,
   TAB-6, TAB-8, TAB-10, BEN-7, CLS-4. Fix once: a validity method plus one
   helper that returns annotations aligned to a vector of sample names.
2. **The string `"NA"` instead of `NA`.** Strand annotation fills unmatched
   rows with `"NA"`, which then passes `!is.na()` and string comparisons
   (IMP-2/TAB-9, IMP-15).
3. **In-place mutation of the caller's object.** About 25 sites use
   `eval.parent(substitute())` or `deparse(substitute()) + assign(parent.frame())`.
   The second form breaks or creates junk variables when the argument is an
   expression (`L$m`, `f()`); see IMP-19, CLS-2, BEN-17. Changing this is an
   API change and needs an ADR.
4. **Copy-pasted blocks.** The result/modality/model validation block appears
   about 17 times across 11 files, and "store a model in `result_list`" five
   times, all with raw `@` access. One `.get_model_checked()` and one
   `.store_model()` helper would remove most `@` uses outside `R/class_*.R`.
5. **Examples run hg38 data through hg19.** `data(musica)`,
   `musica_sbs96_tiny` and friends are hg38 (100% ref match vs 21-49% for
   hg19), but roxygen examples use `select_genome("19")` and
   `annotate_transcript_strand(musica, "19")`, so documented outputs count
   only a fraction of variants (IMP-17, TAB-11).

### Proposed fix PRs

One category per PR, as the standards require. Each adds regression tests for
its findings, which is also how coverage rises (item 3 of the plan).

| PR | Scope | Findings | Changes results | Needs ADR |
|----|-------|----------|-----------------|-----------|
| A | Sample alignment invariant, validity methods, rebuild `data/*.rda` | CLS-1, CLS-4, CLS-7, CLS-8, CLS-9, PLT-1, PLT-5, DIS-5, TAB-6, TAB-8, TAB-10, BEN-7 | Yes | No |
| B | Variant import correctness | IMP-3–IMP-13, IMP-18, IMP-20 | Yes | No |
| C | Variant annotation correctness, hg38 examples | IMP-1, IMP-2/TAB-9, IMP-14–IMP-17, TAB-11 | Yes | No |
| D | IND83 counting (#9) | TAB-1–TAB-5 (with IMP-5 coordinates from B) | Yes | No |
| E | SBS/DBS table robustness | TAB-7 (#52 follow-up), TAB-12, TAB-13, TAB-14 (`rc` speed) | No (crash fixes) | No |
| F | Discovery/prediction correctness | DIS-1–DIS-4, DIS-6–DIS-10, DIS-13, CLS-11 | Yes | No |
| G | Differential analysis and clustering statistics (#34) | PLT-2–PLT-4, PLT-6–PLT-8, PLT-19 | Yes | No |
| H | Plotting fixes (#95) | PLT-9–PLT-18, PLT-20, PLT-21 | Plots only | No |
| I | Benchmarking correctness (#94) | BEN-1–BEN-16 | Yes (benchmark summaries; rebuild `full_benchmark_example`) | No |
| J | Accessor and setter bugs (#48) | CLS-2, CLS-3, CLS-5, CLS-6, CLS-12, CLS-13 | No | No |
| K | Shared helpers, remove `@` outside class files, dead code | Theme 4, CLS-10, CLS-14, CLS-15, DIS-12, PLT-22, TAB-15, IMP-21, BEN-17 | No | Only if it changes the mutate-in-place API (theme 3) |
| L | Dependencies | Dependency audit below | No | Yes (moving Imports to Suggests) |
| M | Efficiency | DIS-11, TAB-14, BEN-9, BEN-18, IMP-22 | No | No |
| — | Stan prerequisites in `result_model` | Stan notes below | No | Yes, with the Stan ADR |

Suggested order: A first (other fixes depend on name-based alignment), then
J, B, C, D, F, G, I, E, H in any order, then K, M and L.

### Findings: variant import and annotation (`R/load_data.R`, `R/annotate_variants.R`)

- **IMP-1 | bug | high | annotate_variants.R:66-75 | Verified (code).**
  `add_flank_to_variants` reverse-complements the *same-side* window for
  purine-reference SBS; the pyrimidine-strand flank is the mirrored window.
  Reverse-change SBS matched the SBS96 context for 128/439 variants (chance);
  forward-change 445/445. Fix: fetch `center - flank_end .. center - flank_start`
  before reverse-complementing.
- **IMP-2 | bug | high | annotate_variants.R:233,245,248 | Verified (code).**
  `annotate_transcript_strand` fills non-genic SBS with the string `"NA"`, so
  `ifelse("NA" == "+", "U", "T")` labels every intergenic SBS "T", and
  `drop_na_variants` never drops them. On `data(musica)` all 297 intergenic
  SBS are "T"; SBS192_Trans totals 879 against 587 genic SBS (TAB-9). Fix:
  `NA_character_`.
- **IMP-3 | bug | high | load_data.R:288-289,348 | Confirmed.**
  `multiallele = "expand"` doesn't recode GT, so every ALT is attributed to any
  sample carrying any ALT (`C>A,T`, S1 `0/1`, S2 `0/2` → both samples get
  both). Fix: select allele indices per sample before `expand()`.
- **IMP-4 | bug | high | load_data.R:856-867 | Verified (code).**
  `convert_dbs` merges `diff(start) == 1` across the whole table without
  checking sample, sorting, or runs of 3+. Sample A chr1:1000000 + sample B
  chr1:1000001 became one DBS for A; a 3-base MNV lost its third base;
  unsorted adjacent SBS weren't merged. Fix: order by sample/chr/start and
  require same sample.
- **IMP-5 | bug | high | load_data.R:824-849 | Verified (code).**
  `standardize_indels` strips the VCF anchor base but leaves `start`/`end`
  unchanged, so the same deletion is `DEL_T_1_5+` from VCF and `DEL_T_1_3`
  from MAF. Root cause of part of TAB-4.
- **IMP-6 | bug | medium | load_data.R:824-849 | Confirmed.** Anchor base is
  assumed, never checked: ref `A`, alt `GGT` is rewritten to `-`>`GT`, losing
  A>G. Treat as complex (relates to #7).
- **IMP-7 | bug | medium | load_data.R:287-294,368 | Confirmed.**
  `multiallele = "exclude"` returns a list column of DNAStringSet for
  CollapsedVCF, and for ExpandedVCF drops every multi-base ALT (all indels),
  121 → 117 rows on the LUAD example.
- **IMP-8 | bug | medium | load_data.R:348 | Confirmed.** VCFs without GT
  (Strelka, tumour-only) return NULL with "All variants matched the reference
  allele".
- **IMP-9 | bug | medium | load_data.R:470-526,542 | Confirmed.** The
  `fix_vcf_errors = TRUE` (default) path never assigns `dt`, so it returns
  `stats::dt`, the t-density function. `skip = 7` is hard-coded.
- **IMP-10 | bug | medium | load_data.R:201,214 | Confirmed.**
  `extract_variants("x.vcf.gz")` fails: `file_ext()` returns `"gz"`.
- **IMP-11 | bug | medium | load_data.R:742; annotate_variants.R:143 | Confirmed.**
  Matrix input is rejected though documented; factor ref/alt crash `nchar()`.
- **IMP-12 | bug | medium | load_data.R:1124,1230 | Confirmed.** Variants on
  contigs absent from the genome pass the check and then crash `getSeq`.
- **IMP-13 | bug | medium | load_data.R:938,943,1002,1007 | Confirmed.**
  `create_musica_from_counts` overwrites rownames with its own motif order, so
  SigProfiler-ordered input is silently mislabelled (reversed rows:
  `C>A_ACA` = 0, truth 4). Compounds DIS-3.
- **IMP-14 | bug | medium | annotate_variants.R:204,212-213 | Confirmed.**
  The custom-TxDb path errors before it is reached, tests a class *name*, and
  would store the TxDb rather than `genes(txdb)`.
- **IMP-15 | bug | medium | annotate_variants.R:295-300 | Confirmed.**
  Replication-strand vector is SBS-length but indexed by full-table row;
  crashes when the last rows aren't SBS. Also `"NA"` strings and `@listData`.
- **IMP-16 | bug | medium | annotate_variants.R:108 | Confirmed.**
  `annotate_variant_length` selects `Variant_Type == "indel"`, but types are
  INS/DEL, so all indel lengths are NA; `-` also counts as one base.
- **IMP-17 | bug | medium | annotate_variants.R:25,200; table_utils.R:177,285; vignettes/musicatk.Rmd:156 | Confirmed.**
  hg38 data run through hg19 resources in examples and the vignette (theme 5).
- **IMP-18 | bug | low | load_data.R:394 | Confirmed.** One `extra_fields`
  column comes back named `V2`; `Number=A`/`.` INFO fields error.
- **IMP-19 | bug | low | load_data.R:1059-1090 | Confirmed.** `add_result`
  with `musica = L$m` creates a variable named `` `L$m` ``; no dimension
  checks; overwrites existing `model_id`; doc says exposures are samples ×
  signatures (stored transposed).
- **IMP-20 | bug | low | load_data.R:810-822 | Confirmed.** Lowercase bases
  are treated as compound and dropped.
- **IMP-21 | quality | low | load_data.R (various) | Confirmed.** Duplicate
  insertion block (833-840, unreachable), dead `used_fields`,
  `length(comp_ins > 0)`, commented-out block 510-526, placeholder roxygen
  ("Aaron - Need to describe...", "Add Description" ×4), progress bar shown
  with `verbose = FALSE` and never closed, wrong defaults in
  `extract_variants` docs.
- **IMP-22 | efficiency | low | load_data.R:1225-1228,1105 | Plausible.**
  `.check_variant_ref_in_genome` copies the whole table four times to read
  four columns.

### Findings: count tables (`R/standard_tables.R`, `R/table_utils.R`, `R/class_count_table.R`)

- **TAB-1 | bug | high | standard_tables.R:509-510,523-555,589,811-816 | Confirmed.**
  `create_ind83_table` crashes unless every indel subclass is present
  (only-DEL input stops in `subset_variant_by_type`; zero-row branches rbind
  undefined `del1_ta`/`ins1_ta`/`ins2_ta`; one ≥2bp deletion class gives a
  `table()` length error). Fix: `table(factor(x, levels = motifs), factor(sample, levels = samples))`.
- **TAB-2 | bug | high | standard_tables.R:833-844 | Verified.** `.micro_left`
  compares against the far end of the flank and returns the wrong length:
  `.micro_left("ACGTT", "TTTT…TTGTT")` = 4, expected 3. Moves deletions
  between MH and repeat classes.
- **TAB-3 | bug | high | standard_tables.R:754-755,796-797 | Verified.** Left
  repeats are counted with a reversed flank but an unreversed unit, so any
  non-palindromic multi-base unit has 0 left repeats:
  `.count_repeat("AC", stri_reverse("TTACACAC"))` = 0, expected 3.
- **TAB-4 | bug | high | standard_tables.R:707-712,746-751,772-779 | Confirmed.**
  Indel flanks assume `start..end` is exactly the deleted bases, which holds
  only for MAF deletions. MAF insertions and all VCF-derived indels skip an
  adjacent base. On `indel_musica`, 12 of 56 1bp insertions are
  miscategorized (e.g. chr1:153056355 +T is `INS_T_1_4`, correct `5+`).
- **TAB-5 | bug | medium | standard_tables.R:746-751 | Confirmed.** Insertion
  flanks are fixed at 10bp, so `INS_repeats_4_5+` and `INS_repeats_5+_5+` are
  nearly unreachable. Use `len * 5` as deletions do.
- **TAB-6 | bug | medium | standard_tables.R:513,585-595 | Confirmed.** IND83
  colnames come from `levels(sample_names())`: character `Samples` gives an
  83×0 table; a different level order silently mislabels columns.
- **TAB-7 | bug | medium | standard_tables.R:261-263,390-401 | Verified.**
  `rc(character(0))` errors, so DBS sets whose refs are all canonical (only
  CC>TT) crash `build_standard_table(..., "DBS78")`. This is the remaining
  crash behind #52. DBS78 collapsing itself is correct (exhaustive test).
- **TAB-8 | bug | medium | table_utils.R:199,313 | Confirmed.**
  `combine_count_tables` rbinds by position without aligning sample columns;
  `build_custom_table` uses a different column order.
- **TAB-9** — same defect as IMP-2, seen from SBS192_Trans.
- **TAB-10 | bug | medium | standard_tables.R:74,344-346,585 | Confirmed.**
  Columns come from stale factor levels: `dbs_musica` DBS78 is 78×668 for 239
  samples, `indel_musica` IND83 has 530 columns for 5 samples.
- **TAB-11 | bug/quality | medium | standard_tables.R:72-74,199-201 | Confirmed.**
  Variants whose motif isn't in the schema are dropped by `xtabs` with no
  message; with the hg19 examples, 47 of 192 SBS are counted for one sample.
- **TAB-12 | quality | low | standard_tables.R:125,237-242,465-470 | Confirmed.**
  SBS192 without `strand_type` errors "argument is of length zero"; the
  overwrite check uses `SBS192` instead of `SBS192_Trans`/`_Rep` (#17).
- **TAB-13 | bug | low | table_utils.R:334 | Confirmed.** `setdiff()` with one
  argument in the missing-table error path.
- **TAB-14 | efficiency | low | standard_tables.R:390-401,135-140,318-338 | Confirmed.**
  Recursive `rc()` takes 5.6 s vs 0.009 s vectorized on 20k strings
  (identical output); SBS192 makes two `getSeq` calls; DBS78 loops over
  samples.
- **TAB-15 | quality | low | various | Confirmed (code).** Dead code
  (`table_96`, `drop_count_table`, `return_table = FALSE` branch, empty `if`,
  `ifelse(ins, plus <- 0, plus <- 0)`), self-referential default in
  `create_dbs78_table`, all-NA IND83 `features`, vignette uses deprecated
  `table_name =`.

### Findings: discovery, prediction, COSMIC (`R/discovery_prediction.R`, `R/compare_results.R`, `R/cosmic_data.R`, `R/k_val_assit.R`)

- **DIS-1 | bug | high | discovery_prediction.R:534-537,581,591-592 | Confirmed.**
  `generate_result_grid` ignores `annotation` and always reads `Tumor_Type`;
  with another column nothing is stored, and with `Tumor_Type` every group
  overwrites the same `lda<k>` model.
- **DIS-2 | bug | high | discovery_prediction.R:574-583 | Confirmed.** A grid
  re-run fetches the old `lda2` after discovery renamed the new fit to
  `lda2.1`; the new fit is silently discarded. Same pattern in BEN-10.
- **DIS-3 | bug | high | discovery_prediction.R:269-277,404,436-452 | Confirmed.**
  `predict_exposure` never checks motif names or modality: a row-permuted
  table predicts silently (top signature SBS1 → SBS9), and DBS78 signatures on
  an SBS96 table give only a recycling warning.
- **DIS-4 | bug | high | discovery_prediction.R:74-75,271-272,585-588; k_val_assit.R:269-272,564-599 | Confirmed.**
  Zero-count samples are dropped without a message, so exposures have fewer
  columns than the table; downstream crashes ("subscript out of bounds", "NA
  in probability vector") or mislabels. `prob[prob == NaN]` is a no-op; no
  `drop = FALSE`.
- **DIS-5 | bug | medium | discovery_prediction.R:887-901,907-911 | Confirmed.**
  `combine_predict_grid` cbinds by group, so columns aren't in sample order;
  the combined model has empty `num_signatures` and `model_id`.
- **DIS-6 | bug | medium | discovery_prediction.R:781-789 | Confirmed.**
  `auto_predict_grid(..., result_name = "pred")` always fails: the temporary
  prediction is stored under `"result"` and read back from `result_name`.
- **DIS-7 | bug | medium | discovery_prediction.R:715,730-739 | Confirmed.**
  Factor annotation columns fail (`model_id` gets a factor).
- **DIS-8 | bug | low | discovery_prediction.R:794-807 | Confirmed.** No
  signature passing thresholds gives "invalid first argument".
- **DIS-9 | bug | low | discovery_prediction.R:82-103,297-318,919-940; benchmarking.R:144-151 | Confirmed.**
  ID de-dup strips two characters; the 12th run is `lda2..11` (also BEN-15).
- **DIS-10 | bug | low | k_val_assit.R:602-615 | Confirmed.** Invalid
  `error_type` only messages and returns NULL after all fits run.
- **DIS-11 | efficiency | medium | k_val_assit.R:404-420 | Confirmed (timing).**
  Vendored `ckmeans` computes a must-link closure with an O(n³) R loop although
  musicatk never passes must-links: n=505 takes 28 s; defaults extrapolate to
  ~11 min of clustering plus ~5,000 NMF fits (#26).
- **DIS-12 | quality | medium | discovery_prediction.R (5 sites); compare_results.R:94-200,233-331,359-434 | Confirmed (code).**
  Model-storage block copied five times; `add_result()` exists but isn't used;
  the three compare functions repeat ~60 lines; `compare_cosmic_v2` title
  default references a non-existent argument; NMF internals read via `@fit@W`
  instead of `basis()`/`coef()`; `.multi_modal_discovery` is dead.
- **DIS-13 | quality | low | discovery_prediction.R:360-421,585-588 | Confirmed.**
  LDA prediction has no `seed`; single-sample prediction emits 100
  non-convergence warnings; the grid's reconstruction error is an
  undocumented squared MSE, not comparable with `compare_k_vals`.

### Findings: classes and accessors (`R/class_musica.R`, `R/class_result_model.R`, `R/class_result_collection.R`, `R/misc_functions.R`)

- **CLS-1 | bug | high | class_musica.R:12,356-381; data/*.rda | Verified.**
  No invariant keeps sample order consistent (theme 1). `res_annot`
  annotations TCGA-94, -56, -77…; exposures TCGA-56, -77, -94….
  `exposure_differential_analysis` gives p = 0.571/0.095/0.857 vs
  1.0/0.190/1.0 when matched by name.
- **CLS-2 | bug | high | class_musica.R:186-189 | Verified.**
  `name_signatures()` returns the character vector, so
  `res <- name_signatures(res, ...)` destroys `res`.
- **CLS-3 | bug | high | class_musica.R:146 | Verified. #48.**
  `subset_variants_by_samples(musica, s[1:2])` returns 161 rows; correct
  (`%in%`) is 319. Silent when lengths divide evenly.
- **CLS-4 | bug | medium | class_musica.R:41-67,91-113 | Confirmed.**
  `subset_musica_by_*` don't subset `result_list`; exposures keep all samples
  and plots gain an NA group.
- **CLS-5 | bug | medium | class_result_model.R:333-354,442-463,553-574,873-893,979-999 | Verified.**
  15 replacement methods are registered for `value = "matrix"` and can never
  dispatch (`num_signatures<-`, `other_parameters<-`, `credible_intervals<-`,
  `model_id<-`, `modality<-`).
- **CLS-6 | bug | medium | class_result_collection.R:189 | Confirmed.**
  `hyperparameter<-` reads the misspelled slot `hperparameter` and never
  assigns.
- **CLS-7 | bug | medium | class_musica.R:414,417 | Confirmed.**
  `sample_names()` reads `$Sample` (partial match of `Samples`); a user column
  named `Sample` hijacks it. Returns a factor, documented as character.
- **CLS-8 | bug | medium | class_musica.R:361-365,415 | Confirmed.** Objects
  from `create_musica_from_counts()` have no sample names, so they can't be
  annotated (tests work around it with `@`).
- **CLS-9 | quality | medium | class files | Confirmed.** No validity methods
  on `musica`, `result_model`, `result_collection`.
- **CLS-10 | quality | medium | class_result_collection.R:50,88,159,233 | Confirmed.**
  Bad names give NULL or "no applicable method for `@`"; the compensating
  validation block is copied about 17 times (theme 4).
- **CLS-11 | bug | low | misc_functions.R:48 | Verified.** `.cosine` drops to a
  vector for one signature; `compare_cosmic_v2` on a 1-signature model errors.
- **CLS-12 | bug | low | class_musica.R:101-104 | Confirmed.**
  `paste(collapse = TRUE)` errors and hides the intended message.
- **CLS-13 | bug | low | class_musica.R:239,358 | Confirmed.** `samp_annot<-`
  rejects Date/POSIXct; `variants<-` rejects data.frame with a dispatch error.
- **CLS-14 | quality | low | discovery_prediction.R (4 sites); load_data.R:1073-1087 | Confirmed (code).**
  No public setter for adding a model; insertion duplicated with raw `@`.
- **CLS-15 | quality | low | misc_functions.R:33; class_result_model.R:1017-1031; umap.R:71 | Confirmed.**
  `.cosineDist` dead; `table_selected()` duplicates `modality()`; a
  `Matrix::Matrix()` dummy call justifies an Import.
- **CLS-16 | quality | low | class_musica.R:415; load_data.R:871 | Confirmed.**
  `gtools::mixedsort` reads `-94` as negative, so TCGA barcodes sort in
  reverse.
- **CLS-17 | quality | low | class files | Confirmed.** No `show()` methods;
  printing a musica dumps the whole variant table.

`@` access outside `R/class_*.R` (code lines): discovery_prediction 24,
kmeans 7, compare_results 6, plotting 4, umap 4, load_data 4, table_utils 3,
and 1-2 in nine other files; tests/testthat/test-kmeans.R 2.

### Findings: downstream analysis and plotting

- **PLT-1** — instance of CLS-1 in differential_analysis.R:75,136-138,158,
  plot_heatmap.R:114-116, kmeans.R:259-263. Confirmed: Kruskal p =
  0.121/1.000/0.699 positional vs 0.439/0.439/0.245 by name; heatmap labels
  the wrong column.
- **PLT-2 | bug | high | differential_analysis.R:182-189 | Verified (code).**
  glm.nb coefficients are flattened column-major (4k + 1 values), and the
  slice takes the last k, i.e. p₂…p_k plus the ANOVA p. Every `(fdr)` column
  is shifted by one group. `A(coef)` is the intercept. The test asserts a
  "NaNs produced" warning from degenerate fits.
- **PLT-3 | bug | medium | differential_analysis.R:145-146 | Confirmed.**
  wilcox with 3+ comparisons crashes (`-1:1` subscript).
- **PLT-4 | bug | medium | differential_analysis.R:150-151 | Confirmed.**
  wilcox output is a character matrix (p-values are strings); column named
  `Pr(>|z|)`.
- **PLT-5 | bug | medium | kmeans.R:259-262 | Confirmed.**
  `plot_cluster(group = "annotation")` ignores `annotation` and plots the
  first annotation column.
- **PLT-6 | bug | medium | kmeans.R:114-115 | Confirmed.** `hcut(diss,
  isdiss = FALSE)` recomputes Euclidean distance on the distance matrix, so
  `dis.method` is ignored.
- **PLT-7 | bug | medium | kmeans.R:40,97-102 | Confirmed. #34.**
  `method = "kmeans"` runs `cluster::fanny` with `maxit = 10`, `tol = 1e-15`,
  which is the cause of #34; `k_select` uses `stats::kmeans`.
- **PLT-8 | bug | medium | differential_analysis.R:216-231 | Confirmed.**
  `plot_differential_analysis` melts character columns for wilcox, plots
  intercepts as effects for glm.nb, computes `signif` on coefficients, and
  returns NULL for other methods. Redesign or deprecate.
- **PLT-9 | bug | medium | plot_heatmap.R:126-139 | Confirmed.**
  `subset_tumor` uses regex selection but exact-match annotation; errors on
  overlapping labels; missing `all_of()`.
- **PLT-10 | bug | medium | plot_heatmap.R:104-105 | Confirmed.**
  `scale = TRUE` z-scores each sample across signatures, not each signature
  across samples.
- **PLT-11 | bug | medium | plot_signatures.R:270-274,289,304-305 | Confirmed.**
  Motif header colours are positional; a custom `color_mapping` mismatches
  header and bars.
- **PLT-12 | bug | medium | plotting.R:45-48 | Confirmed.** A vector
  `color_variable` recycles across signatures (64 of 288 rows right with 3
  signatures); a column name is ignored.
- **PLT-13 | bug | medium | plot_exposures.R:220-250 | Confirmed.**
  `plot_type = "box"`/`"violin"` with `group_by = "none"` silently returns a
  bar plot; the existing test only checks `is.ggplot`.
- **PLT-14 | bug | low | umap.R:317-318,334-335 | Confirmed.** `signif()`
  colour limits drop the top samples to grey.
- **PLT-15 | bug | low | umap.R:263 | Confirmed.** `y = UMAP_1` typo stretches
  the y-axis.
- **PLT-16 | bug | low | umap.R:206-221 | Confirmed.**
  `color_by = "cluster"` demands `annotation`.
- **PLT-17 | bug | low | umap.R:74,85,203,299 | Confirmed.** Unclear errors
  for `plot_umap` before `create_umap` and for zero-exposure samples (0/0
  also in kmeans.R:89, plot_exposures.R:141, plot_heatmap.R:102).
- **PLT-18 | bug | low | plot_exposures.R:156-161 | Confirmed.** `threshold`
  drops rows instead of setting them to 0 as documented.
- **PLT-19 | bug | low | kmeans.R:329,373 | Confirmed.** `k_select` default
  `n = 10` errors on ≤10 samples, including the package's own example.
- **PLT-20 | quality | low | kmeans.R:265,276 | Confirmed. #95.** `aes_string()`.
- **PLT-21 | bug | low | plot_signatures.R:247 | Confirmed.**
  `exists("annotation_text")` picks up a global variable.
- **PLT-22 | quality | low | all nine files | Confirmed.** Validation block
  copied nine times; dead `.addSmallLegend`, `head(sample_annot)`, unused
  `header`; misleading messages.

### Findings: benchmarking (`R/benchmarking.R`, `R/class_single_benchmark.R`, `R/class_full_benchmark.R`)

- **BEN-1 | bug | high | benchmarking.R:1140,1145-1231 | Confirmed.** A
  predicted signature in two duplicate groups is merged into both; total
  exposure grows from 1,514,031 to 1,576,239.
- **BEN-2 | bug | high | benchmarking.R:1347-1359,1383-1388 | Verified (code).**
  Composite correction handles only two components (`factor1`, `factor2`); a
  3-component composite aborts `benchmark()`.
- **BEN-3 | bug | high | benchmarking.R:1321,1345,1354 | Verified (code).**
  Composite correction hard-codes 96 rows and `count_tables$SBS96`; any DBS78
  composite crashes.
- **BEN-4 | bug (methodology) | high | benchmarking.R:1319,1329,1341,1611-1628 | Confirmed.**
  Composites are replaced by scaled ground-truth signatures, so headline
  metrics reward composites (Mean CS 0.972 for a 7-signature prediction with a
  composite vs 0.954 for the clean 8). Must be fixed before comparing methods,
  including the Stan models.
- **BEN-5 | bug | medium | benchmarking.R:1260,1270,1460,1470 | Confirmed.**
  Missing `drop = FALSE` loses signature names when one column survives.
- **BEN-6 | bug | medium | benchmarking.R:625-643 | Confirmed.**
  `benchmark_get_comparison` accepts stage aliases but dispatches only on
  exact names: `"final"` errors.
- **BEN-7 | bug | medium | benchmarking.R:21-26,802-807,901-903 | Confirmed.**
  True exposures are stored samples × signatures (transposed vs everywhere
  else) and never matched by sample name; shuffled rows drop correlation to
  ~0 silently.
- **BEN-8 | bug | medium | benchmarking.R:23-26 | Confirmed.** Modality
  aliases (`"SBS"`) stored as typed, then `benchmark(..., "SBS96")` fails.
- **BEN-9 | bug/efficiency | medium | benchmarking.R:177,202,229 | Confirmed.**
  `plot = FALSE` still draws via `compare_results()` → `grid.arrange`.
- **BEN-10 | bug | medium | benchmarking.R:384-398 | Plausible.**
  `predict_and_benchmark` benchmarks the old model after a rename (DIS-2).
- **BEN-11 | bug | medium | benchmarking.R:1551,1627-1628 | Confirmed.**
  Found and Missing come from different stages and count rows, not unique
  signatures ("Num Found 9" of 8).
- **BEN-12 | bug/quality | medium | benchmarking.R:1521; k_val_assit.R:602 | Confirmed.**
  `.get_reconstruction_error` defined twice; `benchmark()` works only because
  of Collate order.
- **BEN-13 | bug | medium | compare_results.R:41 via benchmarking.R:177 | Confirmed.**
  A prediction with no match above threshold aborts `benchmark()` instead of
  recording a total miss.
- **BEN-14 | quality | low | benchmarking.R:111,117,432,548,615,678,768,856,998 | Confirmed. #94.**
  Nine `class(x)[1] !=` checks, not one.
- **BEN-15** — same as DIS-9.
- **BEN-16 | quality | low | benchmarking.R:247,1586,1636-1657,1739-1770,1795 | Confirmed.**
  Summaries are character matrices; dead `is.infinite`; `min(NULL)` warnings;
  `exists()` without `inherits = FALSE`.
- **BEN-17 | quality | medium | benchmarking.R (file) | Confirmed (code).**
  About a third duplicated or dead: validation block in seven exports,
  discarded ggplots, `benchmark_plot_comparison` copies compare_results.R,
  `assign(parent.frame())`, doc errors (nonexistent `make_copy`, IND83
  claimed but unsupported).
- **BEN-18 | efficiency | low | benchmarking.R:122,191,248,259,1163-1167,1349-1362,1795 | Confirmed (code).**
  Repeated `get_model()` calls; O(runs²) summary rebuild across a k sweep.

### Dependency audit (proposed PR L; needs an ADR)

- Remove: `tidyverse` (no uses; `@import tidyverse` in four plot files),
  `base`.
- `Matrix (>= 1.6.1)`: only a dummy call (umap.R:71), likely a uwot pin.
- Imports used once, candidates for Suggests + `requireNamespace()`:
  `decompTumor2Sig`, `maftools`, `MCMCprecision`, `philentropy`, `shiny`
  (launcher only; all other app packages are already Suggests).
- Unused Suggests: `survival`, `XVector`, `cowplot`.
- Missing: `DT` (22 uses) and `BiocManager` (1) in `inst/shiny`.

### Open issues: status on devel 2b8cb17

| Issue | Status | Notes |
|-------|--------|-------|
| #7 compound INDELs | Still current | Enhancement; COSMIC ID83 has no complex class. IMP-6 shows some complex variants are mangled, not removed. |
| #9 indel miscalculations (labelled invalid) | **Valid** | TAB-2–TAB-5 are concrete miscounts. Recommend removing the label. |
| #17 non-standard table names | Partial | SBS192 naming and the overwrite check (TAB-12). |
| #25 indel plotting | Still relevant | Fix PLT-11 first. |
| #26 parallel gridding | Still relevant | DIS-11 makes `compare_k_vals` the biggest beneficiary. |
| #34 cluster non-convergence | Reproducible | Root cause PLT-7. |
| #38 simple file format | Not implemented | IMP-10 blocks compressed files too. |
| #48 `==` vs `%in%` | Reproducible | CLS-3. |
| #52 NN alleles | Naming is by design | The input still fails, now in `rc(character(0))` (TAB-7). Close and file TAB-7. |
| #72 outdated COSMIC | Valid | Bundled sets are v3.0/3.1 era (67 SBS); COSMIC is v3.5 (Nov 2025). |
| #80 error in plot | Not a musicatk bug | User called `maftools::plotSignatures`. Close with pointer to `plot_signatures()`. |
| #81 user-built COSMIC results | Valid | No exported `result_model` constructor; compare_cosmic functions hard-code bundled sets. |
| #94 class() comparison | Reproducible | Nine sites (BEN-14). |
| #95 deprecated ggplot2 | Reproducible | PLT-20 and `is.ggplot` in two test files. |
| #85, #47, #93, #96, #97 | Not reviewed | Outside `R/` logic (navbar, CI, tarball, Rd anchors) or pure enhancement. |

Finding 1 of the 2026-09-18 entry (33 pkgdown topics) is confirmed resolved:
all 33 are in `_pkgdown.yml`.

### Coverage gaps (for item 3 of the plan)

No tests touch: count-table construction (all four modalities, including the
IND83 internals), annotation, any musica/result_model/result_collection
accessor or setter, subsetting, `add_result`, `create_musica_from_counts`,
discovery or prediction (beyond two input-type errors), COSMIC comparison,
`compare_k_vals`, benchmarking (no test file), `plot_signatures`,
`plot_sample_counts`, `plot_sample_reconstruction_error`, `plot_umap`,
`plot_differential_analysis`. Existing plot and differential tests assert
only `is.ggplot`/`expect_type`, which none of the PLT bugs would trip.
Highest-value invariants: per-sample `colSums` equal variant counts per type;
SBS96/DBS78 reverse-complement symmetry; hand-built IND83 cases on a tiny
synthetic genome; exposures sum to counts; sample order identical across
annotations, tables and exposures; benchmark Found + Missing = number of true
signatures and total exposure preserved.

### Prerequisites for the Stan inference models (item 5 of the plan)

- `result_model` has no slot for algorithm, version, seed or settings; the
  `model_id` string is the only record, and grid code depends on its format.
- `credible_intervals`, `other_parameters` and `metrics` are untyped and never
  filled, and their setters can't dispatch (CLS-5). Define a structure (e.g.
  lower/upper arrays with the same dims as exposures, `probs`) before use.
- No place for posterior draws, a fit object, or diagnostics (Rhat, ESS,
  divergences); `parameter`/`hyperparameter` live per `result_name`, not per
  model.
- Methods are hard-coded `if/else` branches with no registry and no `...` or
  `control` pass-through; fitting and storage are fused. A
  `.fit_<method>(counts, k, ...) -> list(signatures, exposures, extras)`
  contract plus one `.store_model()` would let a Stan backend plug in.
- Prediction must validate motif alignment centrally (DIS-3), and models
  need validity and sample alignment (CLS-9, CLS-1), or every new method
  inherits those bugs.
- Benchmarking needs BEN-4 and BEN-3 fixed before comparing models, accepts a
  stored model only by name, and has no exposure-accuracy metric or room for
  runtime and uncertainty.

---

## 2026-11-06: rep_range work (#220)

Found by `/code-review` on the #220 branch. Not fixed there.

### 1. `inst/shiny/server_tables.R` is dead code

Filed as #226.

`app.R` sources only `ui.R` and `server.R`, and `server.R` defines its own
`add_tables()`. The copy in `server_tables.R` is never used. It has also
drifted: it passes `table_name =` to `build_standard_table()`, which no
longer has that argument, so it would fail if it were ever sourced. Remove
it, or source it in place of the copy in `server.R`.

### 2. Shiny: replication-strand genome check comes after the overwrite prompt

Filed as #227.

In `server.R`, the check that the selected genome has a replication strand
object (hg19 or hg38 only) runs inside `add_tables()`. That is after the
"overwrite existing table?" modal. A user with an SBS192_Rep table who
picks mm10 confirms the overwrite, and only then sees the error. The table
isn't changed.

### 3. Shiny: genome-to-replication-object lookup lives in the app

`server.R` `add_tables()` maps the selected genome to `rep_range` (hg19) or
`rep_range_hg38` (hg38) and shows an error for any other genome. That is
package logic in the app, against the standards' no-logic-in-the-app rule.
`inst/` is invisible to R CMD check, so it is untested. It should move into
the package: for example, `annotate_replication_strand()` could accept a
build name, as `annotate_transcript_strand()` does. The app would then pass
`input$table_genome_list` straight through.

**Resolved in the same branch:** the maintainer approved having
`annotate_replication_strand()` take `"hg19"`/`"hg38"`, and the app now
passes the selected genome straight through.

## 2026-11-10: differential analysis and clustering (PR G of #208)

Found by the spec review of `fix/diff-analysis-clustering`; pre-existing,
left out of that PR.

### 1. `exposure_differential_analysis()` rejects factor `group1`/`group2`

`group1` and `group2` must be character, although the annotation itself
may be a factor. Passing `levels(f)[1:2]` works; `f[1:2]` errors with
"must be character vectors". The Shiny app passes character vectors.

### 2. Misleading message for kruskal and glm.nb with group1/group2

With three or more groups and `method = "kruskal"` or `"glm.nb"`, giving
`group1`/`group2` prints "'annotations' is of length 2. 'group1' and
'group2' were ignored.", though the annotation has more than two values.

### 3. Duplicate wilcox pairs are accepted

The error says pairs "must be unique", but the check only rejects
`group1[i] == group2[i]`, so `group1 = c("A", "A")`,
`group2 = c("B", "B")` runs the same comparison twice.

### 4. `k_select()` doesn't match `cluster_exposure()`'s distances

`k_select()` has no `"fanny"` option, and its `"hclust"` and `"pam"` paths
call `fviz_nbclust()` on the exposure matrix with Euclidean distance,
ignoring any `dis.method`. Now that `cluster_exposure()` uses `dis.method`
for hclust (#173), k picked with `k_select()` may not suit a cosine or
Manhattan clustering. Found by `/code-review` of PR G.

