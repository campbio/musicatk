# musicatk Roadmap

Direction for the package, so contributors and agents know which way proposals
should lean. Maintained by the Campbell Lab; update it when priorities shift.

## Current focus

- **Benchmarking framework.** Added in the v2 line (`R/benchmarking.R`,
  `class_single_benchmark.R`, `class_full_benchmark.R`,
  `class_result_collection.R`, `class_result_model.R`) to compare discovery and
  prediction methods systematically. Still settling; expect API movement.
- **Post-v2 consolidation.** v2.0.0 folded the old `musica_result` class into
  `musica`. Remaining work is making accessor conventions uniform and documenting
  the object model.

## Known technical debt

Tracked but deliberately out of scope for feature work:

- Duplicate pkgdown output committed under `docs/articles/articles/` (~15 MB),
  left from a pre-pkgdown-2.0 build. Needs a clean rebuild.
- Migration of the pkgdown site off the main branch to `gh-pages`.
- ~48 committed macOS `" 2"` duplicate files across `docs/` and `vignettes/`.
- `vignettes/ui_screenshots` (4.7 MB) and `vignettes/figures` (496 KB) ship in
  the package tarball although only pkgdown articles reference them. Relevant to
  the Bioconductor 5 MB tarball guidance.
- Thin test coverage — six testthat files for 28 source files.
- Shiny app has no automated smoke test and no module structure.

## Directions under consideration

<!-- Maintainer: replace this list as priorities are set. Do not let agents
     invent roadmap items — this section should reflect real lab decisions. -->

- Precomputed vignettes (`.Rmd.orig` pattern) so full site builds are fast.
- Shiny module refactor (requires an ADR).

## Out of scope

- Support for mutation types beyond SBS/DBS/INS/DEL without a specific use case.
- Breaking changes to the `musica` object outside a Bioconductor major release.
