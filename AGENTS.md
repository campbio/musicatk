# AGENTS.md

## Campbell Lab Playbook (common across lab packages — v2.0, do not edit per-repo)

### Common commands
make test / make check / make bioccheck / make docs / make lint / make site
(See Makefile for definitions. These are the ONLY sanctioned entry points.
Run `make test` after every change; `make check` before opening a PR.)

### Git and PR workflow
- Branch from devel; all work lands via PR. Never push to devel or master.
- Use plan mode for any non-trivial change.
- Run /code-review before requesting human review.
- Every user-facing change gets a NEWS.md entry.

### Coding conventions
- Style enforced by lintr/styler (config in repo); <= 80-char lines (BiocCheck).
- roxygen2 owns man/ and NAMESPACE — NEVER hand-edit them.
- Use accessor functions, not @ slot access, outside class definition files.

### Documentation (pkgdown)
- The website is GENERATED. Improve docs by editing roxygen comments, vignettes,
  and _pkgdown.yml — never files under docs/ or the gh-pages branch.
- New exported functions MUST be added to the _pkgdown.yml reference index;
  verify with pkgdown::check_pkgdown().
- To preview one changed page: pkgdown::build_article("<name>") or
  build_reference_index(). NEVER run a full build_site() as verification —
  full site builds/deploys are a local maintainer action (make site-deploy).
- Files under vignettes/articles/ are pkgdown-only and NOT checked by
  R CMD check — knit locally when you edit them.

### Shiny app rules (packages with inst/shiny only)
- The app contains NO analysis logic. Server code only wires inputs to
  exported package functions and renders results. New app features are
  implemented as tested, exported functions first.
- Reactive logic is tested with shiny::testServer(); the golden path is
  covered by a small shinytest2 smoke suite (make test-app).
- UI changes are verified with a screenshot of the RUNNING app
  (make app + browser), not just passing tests.
- inst/ code is invisible to R CMD check — tests and lintr are the only
  guards; inst/shiny is included in the lint paths.

### Versioning and releases
- Bioconductor even/odd x.y.z scheme; releases ~April and ~October.
- Follow dev/RELEASE.md for the release checklist.

### Safety rules
- No structural refactors (file splits, DESCRIPTION dependency changes,
  class redesign) without an approved ADR — propose via a GitHub issue.
- Never commit secrets, tokens, or absolute local paths.
- Architectural decisions are recorded in dev/adr/ (see template and index
  there). Never store anything in docs/ — that is pkgdown build output.
- Maintainer docs (release, roadmap, audits) live in dev/, not the root.

## This package: musicatk

### Project overview
musicatk (MUtational SIgnature Comprehensive Analysis ToolKit) is a Bioconductor
package for discovering and analyzing mutational signatures in cancer genomes. It
extracts variants from VCF/MAF/matrix sources, builds SBS/DBS/INS/DEL motif count
tables with genomic context (including replication and transcription strand),
discovers signatures de novo or predicts exposures against known signature sets,
and supports comparison against COSMIC v2/v3 plus clustering, UMAP, differential
analysis and visualization. A Shiny GUI in inst/shiny exposes the same workflow.

### Repository map
- `R/class_*.R` — S4 class definitions (musica, count_table, result_model,
  result_collection, single_benchmark, full_benchmark)
- `R/load_data.R`, `R/annotate_variants.R` — variant import and annotation
- `R/standard_tables.R`, `R/table_utils.R` — motif count table construction
- `R/discovery_prediction.R` — discovery (NMF, LDA) and exposure prediction
- `R/compare_results.R`, `R/cosmic_data.R` — COSMIC v2/v3 comparison
- `R/differential_analysis.R`, `R/kmeans.R`, `R/umap.R`, `R/k_val_assit.R` — downstream
- `R/plot_*.R`, `R/plotting.R` — visualization
- `R/benchmarking.R` — benchmarking framework
- `inst/shiny/` — Shiny app: `app.R` plus paired `ui_*.R` / `server_*.R` per tab
- `inst/extdata/` — example VCF/MAF/count-table fixtures
- `vignettes/musicatk.Rmd` — the package vignette
- `vignettes/articles/` — pkgdown-only articles (not built by R CMD check)
- `tests/testthat/` — testthat suite

### Object model
`musica` is the primary container. Slots: `variants` (data.table), `count_tables`
(list), `sample_annotations` (data.frame), `result_list` (SimpleList). Since
v2.0.0 it also holds discovery/prediction results, replacing the former
`musica_result` class.

Accessors are S4 generics with replacement forms — `variants()`, `tables()`,
`samp_annot()`, `result_list()`, `sample_names()`, `signatures()`, `exposures()`,
`modality()`, and others. Use these, never `@` slot access, outside `R/class_*.R`.

`count_table` holds per-sample motif counts. `result_model` / `result_collection`
hold signature results. `single_benchmark` / `full_benchmark` back `R/benchmarking.R`.

### Environment setup
- R >= 4.4.0. `NMF` is a hard Depends, not just an Import.
- `BiocManager::install("musicatk", dependencies = TRUE)`
- The repo uses renv (`renv/` and `.Rprofile` are gitignored locally).
- Heavy annotation deps: BSgenome hg19/hg38/mm9/mm10 and TxDb hg19/hg38. Large
  downloads; the first install is slow. Do not reinstall casually.
- No `src/` — pure R, no compilation step.

### Package-specific notes
- The Shiny app in `inst/shiny` IS in scope for agents, subject to the
  no-logic-in-server rule above.
- No shinytest2 suite exists yet. `make test-app` is a placeholder until
  `tests/app/` is created, which requires adding shinytest2 to DESCRIPTION
  Suggests — out of scope for the scaffolding PR that introduced this file.
- No Rcpp/C++.
- The testthat suite is thin (differentialanalysis, heatmap, kmeans, load-vcf,
  plotting, signatures). Start the coverage ratchet low and raise it.
- Known repo-hygiene issues, tracked separately and deliberately NOT fixed here:
  duplicate pkgdown output under `docs/articles/articles/`; ~48 committed macOS
  `" 2"` duplicate files; `vignettes/ui_screenshots` (4.7 MB) and
  `vignettes/figures` (496 KB) shipping in the tarball though only
  `vignettes/articles/` references them.
- `docs/` is committed pkgdown output for this package and the gh-pages
  migration has not been done yet. Treat it as generated: never hand-edit it,
  and never run `pkgdown::clean_site()` or `deploy_to_branch()` — both are
  maintainer actions (`make site-deploy`).
- The `bioc` git remote currently points at `packages/celda.git`, not
  `packages/musicatk.git`. Fix before any release push — see dev/RELEASE.md.
