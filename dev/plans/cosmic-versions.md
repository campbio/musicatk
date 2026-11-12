# Spec: Finish #253, versioned COSMIC signatures (ADR 0003)

Changes results: yes. `compare_cosmic_v3()` with default arguments on
genome data now compares against COSMIC v3.6 instead of v3.0. Code that
uses the data objects directly (`cosmic_v3_sbs_sigs` and the others, for
example in `predict_exposure()`) gives the same results as before. Devel
only.

## Where #253 stands

- Its code still updates `cosmic_v3_sbs_sigs`, `cosmic_v3_dbs_sigs` and
  `cosmic_v3_indel_sigs` in place to v3.6.
- It adds `cosmic_v3_sbs_sigs_grch38` and `cosmic_v3_dbs_sigs_grch38`.
- It adds `genome_build` to `compare_cosmic_v3()`.
- ADR 0003, which the maintainer accepted on 2026-11-12, replaces this
  with versioned objects. Its "Consequences" section lists the changes
  #253 needs.
- Conflicts with `devel` (9967143) are only `DESCRIPTION` and `NEWS.md`.

## Approach

**Branch.** Work on `feature/cosmic-v3.6` and push to Natasha's fork, so
#253 stays the PR. Bring `devel` in with a merge, not a rebase: the
standards don't let agents force-push.

**ADR.**
- Set ADR 0003's status to Accepted.
- Deciders: Joshua Campbell (maintainer), on 2026-11-12.
- Update the index row in `dev/adr/README.md`.
- Its decision text is unchanged.

**Data.**
- Restore `data/cosmic_v3_sbs_sigs.rda`, `cosmic_v3_dbs_sigs.rda` and
  `cosmic_v3_indel_sigs.rda` from `devel`, byte for byte. They are v3.0
  again.
- `data-raw/make_cosmic_v3_6.R` writes five new objects. It keeps its
  checks: every motif present, and each signature sums to 1.

  | Object | Model id |
  |---|---|
  | `cosmic_v3_6_sbs_sigs_grch37` | `Cosmic_v3.6_SBS_Genome_GRCh37` |
  | `cosmic_v3_6_sbs_sigs_grch38` | `Cosmic_v3.6_SBS_Genome_GRCh38` |
  | `cosmic_v3_6_dbs_sigs_grch37` | `Cosmic_v3.6_DBS_Genome_GRCh37` |
  | `cosmic_v3_6_dbs_sigs_grch38` | `Cosmic_v3.6_DBS_Genome_GRCh38` |
  | `cosmic_v3_6_indel_sigs` | `Cosmic_v3.6_Indel_Genome` |

- It reads the motif order from the v3.0 objects, which it no longer
  overwrites.
- This removes the unreleased `data/cosmic_v3_sbs_sigs_grch38.rda` and
  `cosmic_v3_dbs_sigs_grch38.rda`, which are renamed as the ADR says.
  `make docs` removes their Rd files.
- The four frozen v3.0 objects have no `data-raw/` script, and none is
  added; they are frozen as shipped.

**`compare_cosmic_v3()`.**
- New argument `cosmic_version = "latest"`, placed after `genome_build`.
  It accepts `"latest"`, `"3.0"` or `"3.6"` (`3.6` as a number and
  `"v3.6"` also work).
- `genome_build = "hg19"` and the build spellings are as in #253, case
  insensitive:
  - `hg19`, `19` or `GRCh37`;
  - `hg38`, `38` or `GRCh38`.
- An internal table `.cosmic_catalogues()` in `R/compare_results.R`
  replaces the `if`/`else` chain. Each row gives a version, modality,
  sample type, build and object name. Indel's build is NA, meaning any
  build.

  | Version | Modality | Sample type | Build | Object |
  |---|---|---|---|---|
  | 3.0 | SBS96 | exome | GRCh37 | `cosmic_v3_sbs_sigs_exome` |
  | 3.0 | SBS96 | genome | GRCh37 | `cosmic_v3_sbs_sigs` |
  | 3.0 | DBS78 | genome | GRCh37 | `cosmic_v3_dbs_sigs` |
  | 3.0 | IND83 | genome | any | `cosmic_v3_indel_sigs` |
  | 3.6 | SBS96 | genome | GRCh37 | `cosmic_v3_6_sbs_sigs_grch37` |
  | 3.6 | SBS96 | genome | GRCh38 | `cosmic_v3_6_sbs_sigs_grch38` |
  | 3.6 | DBS78 | genome | GRCh37 | `cosmic_v3_6_dbs_sigs_grch37` |
  | 3.6 | DBS78 | genome | GRCh38 | `cosmic_v3_6_dbs_sigs_grch38` |
  | 3.6 | IND83 | genome | any | `cosmic_v3_6_indel_sigs` |

- Objects are fetched with `getExportedValue("musicatk", name)`, which
  works for lazy data.
- **`"latest"`** picks the highest version with a row matching the sample
  type, modality and build. For genome data that is v3.6. For exome it is
  v3.0.
- **An explicit version without a matching row is an error.** Examples:
  `"3.0"` with `hg38` genome, or `"3.6"` with exome. The error lists the
  shipped catalogues for that sample type and modality, as version and
  build. There is no fallback.
- **Exome with `hg38` is also an error**, even with `"latest"`. musicatk
  ships only the v3.0 GRCh37 exome catalogue. #253 ignored the build for
  exome; the ADR rules that out ("no silent fallback").
- **Every call reports what it used with `message()`.** Example: "Comparing
  to COSMIC v3.6 SBS96 genome signatures (GRCh37)." For Indel it says
  "(one catalogue for all builds)". For exome it adds that COSMIC stopped
  publishing exome catalogues after v3.0.
- This replaces #253's exome `warning()` and its Indel build message.
  With the version always reported, they would repeat it, and a warning
  in every exome call (as in the vignette) is noise.
- The checks for modality and sample type keep their messages.

**Docs.**
- `R/cosmic_data.R`: the four v3.0 objects get titles and text that state
  "COSMIC v3.0" and GRCh37. For example, the exome object says it is the
  last exome catalogue COSMIC published. Each points to the v3.6 objects
  and `compare_cosmic_v3(cosmic_version =)`.
- The five v3.6 objects get one Rd each.
- `_pkgdown.yml` lists the five new objects and drops the two `_grch38`
  names.
- `compare_cosmic_v3()` documents `genome_build` and `cosmic_version`,
  including the default, the message and the errors. Its example keeps
  `genome_build = "hg38"`, because `res` is hg38.

**Shiny (`inst/shiny/server.R`).**
- The comparison tab keeps #253's genome build input, set from the
  import genome. It adds a "COSMIC version" input with the choices
  "Latest (3.6)" (`"latest"`) and `"3.0"`, passed to
  `compare_cosmic_v3()` as is.
- The prediction tab lists the v3.6 objects (SBS and DBS for each build,
  plus Indel) and the v3.0 objects, with the version and build in each
  label. The default is v3.6 SBS for the import genome's build, or
  GRCh37 if the build is unknown.
- #253's `fill_aet()` fallback for aetiology labels, and the SBS22a/b and
  SBS40a/b/c labels, apply to all the objects.
- No other changes to the app.

**Vignette.** In `musicatk.Rmd`, rewrite #253's sentences:
- the exome comparison uses the v3.0 exome catalogue, the last one COSMIC
  published, and the message says so;
- for genome data, the default is the latest version (v3.6), with
  `genome_build` set to the variants' build;
- `cosmic_version = "3.0"` reproduces earlier results.

No code chunk changes, so the vignette's results don't change. The
article `tutorial_tcga_cnsl.Rmd` only mentions `compare_cosmic_v3()` in
text, and is unchanged.

**NEWS and version.**
- Replace #253's 2.7.14 section with a "# musicatk 2.7.16" section that
  describes the versioned design. Keep #253's Shiny "NULL" label fix.
- Version 2.7.15 -> 2.7.16.

## Tests

Rework `tests/testthat/test-compare-cosmic.R`, keeping #253's helpers
(`.quiet_compare()` and `.expected()`):

- The v3.0 objects are unchanged:
  - SBS 67, DBS 11, Indel 17 and exome 65 signatures;
  - model ids `Cosmic_v3_SBS_Genome`, `Cosmic_v3_DBS_Genome`,
    `Cosmic_v3_Indel_Genome` and `Cosmic_v3_SNV_Exome`.
- The v3.6 objects:
  - SBS 101, DBS 22 and Indel 25 signatures, with the model ids above;
  - the two builds have the same signature names but different
    probabilities;
  - empty exposures.
- Default (`"latest"`) genome SBS:
  - `hg19` gives `cosmic_v3_6_sbs_sigs_grch37`;
  - `hg38` and `GRCh38` give `_grch38`;
  - the message names v3.6 and the build.
- `cosmic_version = "3.0"` with `hg19` gives `cosmic_v3_sbs_sigs`, the old
  result; `3.0` as a number works too.
- DBS: v3.6 for each build, and v3.0 for `hg19`.
- Indel: v3.6 for both `hg19` and `hg38`, and v3.0 on request; the
  message says one catalogue.
- Exome: `"latest"` gives the v3.0 exome object, with a message naming
  v3.0. Exome with `"3.6"` errors, and so does exome with `hg38`. Each
  error lists what is shipped.
- `"3.0"` with `hg38` errors and lists v3.0 GRCh37 and v3.6
  GRCh37/GRCh38.
- An unknown version (`"2.0"`), a bad build (`"mm10"`) and a `NULL`
  build each error.

These are tests for a new feature: the objects and the argument don't
exist yet. Before the code changes, they fail with "object not found" or
"unused argument" rather than on an assertion. The ones that check the
v3.0 objects keep their v3.0 contents, and that "3.0" reproduces the old
result, fail on their assertions on the current branch.

## Checks

- `make test`.
- `make coverage`, which must not drop below `devel`.
- `make lint` on changed lines.
- `make docs`.
- `make site-check`.
- `make article FILTER=musicatk` to render the vignette.
- `make check-full` and `make bioccheck` before handing off.
- Shiny: you check the two tabs in `make app`. I won't launch the app
  without asking.

Out of scope:
- mouse catalogues;
- `data-raw/` scripts for the frozen v3.0 objects;
- marking the unversioned names superseded (the ADR leaves this for
  later).

## Added in review

Spec review:
- Shiny: a failed comparison no longer ends the session.
- The comparison tab has one "COSMIC V3" choice.
- The plot title names the COSMIC version.
- SBS22c has a label.
- Stronger tests and fuller argument docs.

`/code-review`:
- `genome_build` and `cosmic_version` are the last arguments, not after
  `modality`, so existing positional calls keep working.
- The result's `"cosmic_catalogue"` attribute names the object used. The
  app labels the comparison table with it, instead of with "latest".
- `genome_build` also takes the numbers 19 and 38, as `select_genome()`
  does.
- The app's version choice reads "Latest" rather than "Latest (3.6)", so
  the label doesn't go stale when a version is added.

Declined in `/code-review`, for the maintainer:
- **Exome with `hg38` is an error.** The maintainer approved this as
  specified (ADR: no silent fallback). The reviewer noted that users with
  hg38 exome data must pass `hg19` or omit the argument.
- **`genome_build` defaults to `"hg19"`, so hg38 data compared without it
  uses GRCh37 v3.6.** This is #253's and the ADR's default. The message
  names the build.
- **The app's build default falls back to hg19 when `vals$genome` is
  unset** (for example after an `.rds` upload). The input is visible and
  editable. `vals$genome` is set before the import succeeds, which is
  existing app behaviour.
- **The app's catalogue list and default are hard-coded rather than read
  from `.cosmic_catalogues()`.** That would mean the app calling an
  internal function. A new release needs app edits, as today.
- **A shared build-alias helper for `select_genome()`,
  `.select_rep_range()` and `.cosmic_build()`.** This is cleanup across
  files (group K).
- **Lazy loading of the 8 COSMIC objects in the app**, about 200 KB per
  session.
