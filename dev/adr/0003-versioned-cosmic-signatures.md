# 0003. Ship each COSMIC signature version as its own data objects

- **Status:** Accepted
- **Date:** 2026-11-12
- **Deciders:** Joshua Campbell (maintainer), 2026-11-12. Proposed by
  Natasha Gurevich on #253, following a suggestion from the maintainer.

## Context

musicatk ships the COSMIC v3 signatures as data objects
(`cosmic_v3_sbs_sigs`, `cosmic_v3_dbs_sigs`, `cosmic_v3_indel_sigs`,
`cosmic_v3_sbs_sigs_exome`) and compares results to them with
`compare_cosmic_v3()`. Users also pass the objects directly to other
functions, for example
`predict_exposure(musica, "SBS96", cosmic_v3_sbs_sigs, "lda")`.

Until now these objects held COSMIC v3.0 (GRCh37). #253 first updated them
in place to v3.6 (May 2026), which adds signatures (SBS 67 -> 101,
DBS 11 -> 22, Indel 17 -> 25) and revises the probabilities of existing
ones. Updating in place silently changes the results of every existing
script after an upgrade, with no way to get the earlier results back from
the package. Published analyses cite a specific COSMIC version, and
reproducing them needs that version. COSMIC releases a new version every
year or two, so this will recur.

Three more facts shape the design:

- **Builds.** COSMIC publishes SBS and DBS catalogues normalized to
  several genomes: GRCh37, GRCh38, mm9 and mm10 for both v3.0 and v3.6,
  plus mm39 for v3.6 SBS. The objects musicatk has shipped are GRCh37
  only, COSMIC's default.
  Indel signatures are not normalized to genome sequence content, so
  COSMIC publishes one Indel catalogue per version (labelled GRCh37).
- **Exome.** musicatk's genome catalogues come from COSMIC's downloads
  page. Exome-normalized SBS and DBS catalogues for each v3 version and
  build are published in SigProfilerAssignment, maintained by the COSMIC
  signature authors. Its v3.6 genome files are identical to COSMIC's
  downloads.
- **Size.** One full version (SBS and DBS for two builds, plus Indel) is
  about 180 KB of `.rda` files. The source tarball is 8.2 MB against
  Bioconductor's 10 MB limit, so several versions fit.

## Decision

We will keep every COSMIC version musicatk has shipped as its own data
objects, add each new version alongside the old ones, and let functions
choose a version with an argument that defaults to the latest.

1. **Existing objects are frozen.** `cosmic_v3_sbs_sigs`,
   `cosmic_v3_dbs_sigs`, `cosmic_v3_indel_sigs` and
   `cosmic_v3_sbs_sigs_exome` stay COSMIC v3.0 permanently and keep their
   names. Their documentation states the version and build.
2. **New versions get versioned names**, with the build always explicit for
   SBS and DBS and no build for Indel:

   ```
   cosmic_v3_6_sbs_sigs_grch37    cosmic_v3_6_sbs_sigs_grch38
   cosmic_v3_6_dbs_sigs_grch37    cosmic_v3_6_dbs_sigs_grch38
   cosmic_v3_6_indel_sigs
   ```

   Exome catalogues add `_exome` (`cosmic_v3_6_sbs_sigs_grch37_exome`).

   A later version follows the same pattern (`cosmic_v3_7_...`), as would
   other genomes (`..._mm10`). v3.0 stays exactly the four frozen objects:
   it is kept to reproduce earlier musicatk results, which only ever used
   its GRCh37 catalogues.
3. **`compare_cosmic_v3()` gains `cosmic_version`.** The default is the
   latest version that has the requested catalogue (sample type, modality
   and build), and the function reports which version and build it used.
   Today that is v3.6 for genome and exome data. An explicit
   `cosmic_version` for which musicatk ships no catalogue matching the
   request, such as `cosmic_version = "3.0", genome_build = "hg38"`, is an
   error that lists
   what is available; there is no silent fallback to another build or
   version. Which of COSMIC's catalogues musicatk ships
   for each version is decided per release, not by this ADR.
4. **The Shiny app** defaults to the latest version the same way and offers
   the earlier versions.

## Consequences

- Upgrading musicatk no longer changes results for code that uses the data
  objects directly. Code that calls `compare_cosmic_v3()` without
  `cosmic_version` does see new matches when a new version ships, but the
  function says which version it used, and `cosmic_version` restores the
  earlier results.
- Adding a COSMIC release becomes additive: new data files, their
  documentation and a new entry in the version table, with no change to
  existing objects and no new function.
- `compare_cosmic_v3()` needs an internal table of available catalogues
  (version, modality, sample type, build) in place of its current
  `if`/`else` chain. Its error message and its "latest available" default
  both come from that table.
- More exported data objects to document and list in `_pkgdown.yml`, about
  five per release. The tarball grows by about 180 KB per release; when it
  approaches the size limit, a later ADR decides whether to drop old
  versions.
- The unversioned names `cosmic_v3_*` read as "current" but mean v3.0.
  Their documentation must say so; whether to mark them superseded in
  favour of the versioned names is left for later.
- #253 changes accordingly: restore the v3.0 objects, add the v3.6 objects
  under the new names (renaming the unreleased `cosmic_v3_sbs_sigs_grch38`
  and `cosmic_v3_dbs_sigs_grch38`), add `cosmic_version`, and update the
  tests, the Shiny app, the vignette and NEWS. Mouse catalogues are a
  separate feature, since `compare_cosmic_v3()` would also need to accept
  mouse builds.

## Alternatives considered

- **Update the existing objects in place to the latest version** (the first
  version of #253). Simplest, but silently changes results for every
  existing script on upgrade, and earlier versions are no longer available
  from the package.
- **One function per version** (`compare_cosmic_v3_0()`,
  `compare_cosmic_v3_6()`). Pins the version in the call, but adds a
  function with duplicated logic for every COSMIC release and gives the
  Shiny app nothing to offer as a choice.
- **Default `cosmic_version` to the oldest version (v3.0).** Nothing changes
  for anyone unless they opt in, but new users get 2019 signatures unless
  they know to ask, and the gap widens with every release.
- **Leave the build implicit for the default build**
  (`cosmic_v3_6_sbs_sigs` for GRCh37). Shorter, but readers must remember
  that a missing suffix means GRCh37, and it breaks the pattern once other
  genomes are added. COSMIC's own file names always state the build.
- **Add the v3.0 GRCh38 catalogues** (`cosmic_v3_0_sbs_sigs_grch38`,
  `cosmic_v3_0_dbs_sigs_grch38`), which COSMIC publishes. This would let
  `cosmic_version` and `genome_build` combine freely, but v3.0 is kept to
  reproduce earlier musicatk results, which never had a GRCh38 option, so
  these objects would reproduce nothing. They would also sit awkwardly next
  to the frozen, unversioned GRCh37 names. New analyses on hg38 get v3.6
  GRCh38 by default.
- **Dots in the version (`cosmic_v3.6_indel_sigs`).** Reads exactly like
  COSMIC's version number, but the package's lint configuration requires
  snake_case names and flags any assignment to a dotted name, every
  existing data object uses underscores, and a dot in an R name
  conventionally marks an S3 method. Users mostly give the version as
  `cosmic_version = "3.6"`, which keeps the dot.
- **Fall back to GRCh37 (with a warning) when the requested build is
  missing for a version.** More forgiving, but returns a comparison against
  a build the user did not ask for.
- **Store the versions internally (`R/sysdata.rda`) and expose them only
  through an accessor function.** Keeps the export list short, but users
  could no longer pass a specific version to `predict_exposure()` by name,
  and it would be a different pattern from every other dataset in the
  package.

## Corrections

- **2026-11-12, before this ADR was merged (#253).** The accepted text said
  COSMIC stopped publishing exome-normalized catalogues after v3.0, and gave
  `cosmic_version = "3.6", sample_type = "exome"` as an example error. Both
  were wrong: SigProfilerAssignment publishes v3.6 exome catalogues for
  GRCh37 and GRCh38. The maintainer decided to ship the v3.6 SBS exome
  catalogues in #253 (DBS exome is #256). The Context bullet, the naming
  example and decision point 3 were corrected in place. The decision itself
  is unchanged.
