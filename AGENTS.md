# musicatk: notes for coding agents

The shared development standards (r-bioc-dev-standards) load automatically
at session start. This file adds only what is specific to this package.

## About

musicatk (MUtational SIgnature Comprehensive Analysis ToolKit) is a
Bioconductor package for discovering and analyzing mutational signatures in
cancer genomes. It extracts variants from VCF/MAF/matrix sources, builds
SBS/DBS/INS/DEL motif count tables with genomic context (including
replication and transcription strand), discovers signatures de novo or
predicts exposures against known signature sets, and supports comparison
against COSMIC v2/v3 plus clustering, UMAP, differential analysis and
visualization. A Shiny GUI in `inst/shiny` exposes the same workflow.

## Layout

- `R/class_*.R`: S4 class definitions (`musica`, `count_table`,
  `result_model`, `result_collection`, `single_benchmark`, `full_benchmark`)
- `R/load_data.R`, `R/annotate_variants.R`: variant import and annotation
- `R/standard_tables.R`, `R/table_utils.R`: motif count table construction
- `R/discovery_prediction.R`: discovery (NMF, LDA) and exposure prediction
- `R/compare_results.R`, `R/cosmic_data.R`: COSMIC v2/v3 comparison
- `R/differential_analysis.R`, `R/kmeans.R`, `R/umap.R`, `R/k_val_assit.R`:
  downstream analysis
- `R/plot_*.R`, `R/plotting.R`: visualization
- `R/benchmarking.R`: benchmarking framework
- `inst/shiny/`: Shiny app, `app.R` plus paired `ui_*.R` / `server_*.R` per
  tab. In scope for agents, under the standards' no-logic-in-the-app rule.
- `inst/extdata/`: example VCF/MAF/count-table fixtures
- `vignettes/musicatk.Rmd`: the package vignette
- `vignettes/articles/`: pkgdown-only articles (not built by R CMD check)
- `dev/agent-log.md`: findings awaiting triage into GitHub issues
- Pure R: no `src/`, no Rcpp, no compilation step.

**Object model.** `musica` is the primary container. Slots: `variants`
(data.table), `count_tables` (list), `sample_annotations` (data.frame),
`result_list` (SimpleList). Since v2.0.0 it also holds discovery/prediction
results, replacing the former `musica_result` class. `count_table` holds
per-sample motif counts; `result_model` / `result_collection` hold signature
results; `single_benchmark` / `full_benchmark` back `R/benchmarking.R`.

Accessors are S4 generics with replacement forms: `variants()`, `tables()`,
`samp_annot()`, `result_list()`, `sample_names()`, `signatures()`,
`exposures()`, `modality()`, and others. Use these, never `@`, outside
`R/class_*.R`.

**Naming.** Exported functions and arguments are snake_case (`.lintr`
enforces snake_case). Internal helpers are prefixed with `.`. Code is
indented with 2 spaces (see Overrides).

**pkgdown.** `docs/` is committed pkgdown output; the migration to a
`gh-pages` branch hasn't been done. Treat it as generated: never hand-edit
it, and never run `pkgdown::clean_site()` or `pkgdown::deploy_to_branch()`.
To preview one page, use `pkgdown::build_article("<name>")` or
`pkgdown::build_reference_index()`.

**Known repo-hygiene issues**, tracked separately; don't fix them in
unrelated PRs:
- Duplicate pkgdown output under `docs/articles/articles/`.
- 48 committed macOS `" 2"` duplicate files, all under `docs/articles/`.
- `vignettes/ui_screenshots` (4.7 MB) and `vignettes/figures` (496 KB) ship
  in the tarball though only `vignettes/articles/` references them.

## Tests

- `tests/testthat/` is thin: differentialanalysis, heatmap, kmeans,
  load-vcf, plotting, signatures. Coverage starts low; raise it, never let
  it drop.
- `inst/` is invisible to R CMD check, so tests and lintr are the only
  guards on the Shiny app. `make lint` covers `inst/shiny`.
- The Shiny golden path is meant to be covered by a small shinytest2 smoke
  suite in `tests/app/`, run with `make test-app`. It doesn't exist yet, so
  `make test-app` is a placeholder. Creating it needs shinytest2 in
  DESCRIPTION Suggests, which needs an ADR.
- Lint backlog: `make lint` reports about 2,150 lints. CI
  (`.github/workflows/lint.yaml`) lints only the R files a PR changes, so
  new and edited code must be clean.

## Extra make targets

- `site-check` is a standard target; see the standards.
- `test-app`: runs the shinytest2 suite (placeholder for now). Safe.
- `build`: builds the source tarball in the repo root. Ask first.
- `app`: launches the Shiny app for visual verification. Ask first.
- `style`: people only (denied in `.claude/settings.json`). It restyles the
  whole package, which conflicts with "styler on new files only".
- `site`, `site-deploy`: people only (denied). Agents never build or deploy
  the pkgdown site.
- `clean`: people only (denied). It deletes files.

## Setup in a new worktree

- R >= 4.4.0. `NMF` is a hard Depends, not just an Import.
- Dependencies: `BiocManager::install("musicatk", dependencies = TRUE)`.
  The BSgenome hg19/hg38/mm9/mm10 and TxDb hg19/hg38 annotation packages are
  large and the first install is slow; don't reinstall them casually.
- renv is local-only: `renv/`, `renv.lock` and `.Rprofile` are gitignored,
  so there is no shared lockfile to restore. A fresh worktree uses whatever
  R library the developer's setup provides.

## Related packages

None.

## Overrides

- Indentation is 2 spaces, not Bioconductor's 4: `.lintr` sets
  `indentation_linter(indent = 2L)` to match the existing code. Converting
  would be a separate restyling PR.
- In addition to the standards' review step, run `/code-review` on the
  branch before hand-off (it's a PR template item).
