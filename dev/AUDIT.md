# Dependency and Deprecation Audit

A periodic sweep for dependency drift and deprecated API usage. Run it at least
once per Bioconductor cycle, ideally a month before the freeze so findings can be
triaged without schedule pressure.

## How to run

Give an agent this prompt from a clean checkout of `devel`:

> Run `BiocCheck::BiocCheck()` and `BiocManager::valid()` on this package.
> Use r-lib lifecycle practices to find deprecated functions or S4 methods
> and report replacements compatible with the current Bioconductor release.
> Write findings to `dev/agent-log.md`. Do not change code — findings become
> GitHub issues (and ADRs where structural).

## Rules

- **The audit never changes code.** It produces findings only. Mixing discovery
  with remediation makes both harder to review.
- Findings go to `dev/agent-log.md`, then become GitHub issues.
- Anything structural — dropping a dependency, replacing a deprecated class,
  restructuring to satisfy a new check — needs an ADR before implementation.
  See `dev/adr/README.md`.

## What to look at

- Deprecated or defunct functions in Bioconductor dependencies, especially
  `SummarizedExperiment`, `GenomicRanges`, `VariantAnnotation`, `BSgenome`,
  and `GenomicFeatures`.
- Packages in DESCRIPTION that are no longer used, or used only in one place.
  musicatk has a large Imports list; verify each still earns its place.
- Suggests-only packages used unconditionally in tests or vignettes.
- New BiocCheck rules since the last cycle — these appear as NOTEs first and
  become WARNINGs later.
- R version floor (`Depends: R (>= 4.4.0)`) versus what Bioconductor currently
  requires.

## Log

| Date | Run by | Findings | Issues opened |
|---|---|---|---|
| _(none yet)_ | | | |
