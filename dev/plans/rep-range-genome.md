# Spec: build-tagged replication strand data, hg38 object (#220)

The maintainer approved this on 2026-11-06. Tag `rep_range` as hg19, fix
its 0-based starts, add an hg38 liftover as a second object, and use the
object that matches the build in the examples and the Shiny app. The
purine-reference strand flip (#222) gets its own PR. Storing the genome
build in `musica` is #225.

## Problem

`rep_range` has hg19 coordinates, but it doesn't record that: there is no
genome and no chromosome lengths in its seqinfo. Every bundled dataset is
hg38. The examples and the Shiny app pair the two, so replication strand
comes from the wrong coordinates and nothing warns about it. There is no
script that builds the object, so nobody could rebuild or check it.

## What the data is (verified 2026-11-06)

`rep_range` is GEO GSM4305465 (`WT_siLIG1_rep1_rfd_10kb.bw`, GSE134225):
HCT116 cells with ligase 1 depleted, GLOE-seq, hg19, 10 kb bins.

- Every one of its 280,760 ranges matches a range in that file.
- Each label is the sign of the bin's RFD: `lagging` where RFD > 0,
  `leading` where RFD < 0. Bins where RFD = 0 were dropped.
- GEO defines RFD = (REV - FWD) / (REV + FWD).
- In these cells the nicks lie on the newly made lagging strand. A
  `lagging` bin is therefore one where the plus (reference) strand is
  the lagging-strand template. The full argument is in the #222 comment.

The starts are 0-based. For example, bins read 80000-90000 and
90000-100000, so 278,674 positions fall in two bins, and 62,530 of those
pairs have different labels.

## Changes

- **`data-raw/rep_range.R`** (new, build-ignored). It takes the place of
  the missing original script.
  - Downloads GSM4305465 and labels each bin by the sign of RFD. The
    file has 280,770 non-zero ranges: the old object had merged 10 pairs
    of neighbouring bins with the same label, so the check compares the
    positions each label covers.
  - Keeps chr1-22, X, Y and M, as the current object does, and drops
    bins where RFD = 0.
  - Makes the starts 1-based.
  - Checks the result equals the current `rep_range` with start + 1.
  - Sorts the bins and sets seqinfo from `BSgenome.Hsapiens.UCSC.hg19`
    (genome `hg19`, chromosome lengths).
  - Makes `rep_range_hg38` with rtracklayer `liftOver()` and the UCSC
    `hg19ToHg38.over.chain`, downloaded into a temporary folder:
    - Bins are lifted as `+` strand. Pieces that land on `-` sit in an
      inverted region, so their label is swapped: fork direction
      reverses relative to the coordinates.
    - Abutting pieces of the same bin are merged.
    - Where pieces of different bins overlap, the overlap is kept only
      if the labels agree. Bins that don't map are dropped.
    - Neighbouring ranges with the same label are merged, which keeps every
      label and cuts the file to about a third.
    - Seqinfo is set from `BSgenome.Hsapiens.UCSC.hg38`.
  - Both objects keep the metadata column name `X`, so existing code
    that reads `rep_range$X` still works.
  - rtracklayer is used only by this script, so DESCRIPTION is
    unchanged.
- **`data/rep_range.rda`** is regenerated (105 KB), and
  **`data/rep_range_hg38.rda`** is new (272 KB, after neighbouring ranges
  with the same label are merged).
- **`R/cosmic_data.R`**: the docs for `rep_range` give the source sample,
  the build, the meaning of the labels (RFD, not timing), and the
  1-based bins. New docs for `rep_range_hg38` describe the liftover and
  what was dropped.
- **`R/annotate_variants.R`, `R/standard_tables.R`**: the `rep_range` docs
  name both objects. The examples use `rep_range_hg38` with the hg38
  `musica`, so the hg19 caveat comment goes.
- **Shiny** (`inst/shiny/server.R`; `server_tables.R` is never sourced, so it is left alone): use `rep_range` for
  hg19 and `rep_range_hg38` for hg38. For other genomes, show an error
  instead of annotating.
- **NEWS, version 2.7.6.**

## Effect on results

- For hg19 variants, labels change only at bin edges, about 1 in 10,000
  positions. Before, an edge position could match two bins and take
  whichever label came last. Now it has one bin.
- With hg38 variants and `rep_range_hg38`, SBS192_Rep tables differ from
  the old (wrong) hg19-on-hg38 tables. The examples' output changes.
- Nothing else changes.

## Tests

In `tests/testthat/test-rep-range.R`:
- `genome()` is `hg19` for `rep_range` and `hg38` for
  `rep_range_hg38`, with chromosome lengths set.
- Neither object has overlapping ranges.
- `rep_range` labels chr8:80001-90000 `lagging`. A variant at
  chr8:90000 gets that bin's label, and annotation stays the same for an
  interior position.
- In `rep_range_hg38`, a bin that lifts over 1:1 keeps its label.
- A bin in an inverted region gets the swapped label at its hg38
  position.

## Checks

- `make docs`, `make test`, `make coverage`, `make check-full` and
  `make bioccheck`.
- Lint the changed lines.
- Confirm the tarball size stays within Bioconductor limits.

## Out of scope

- The strand flip for purine-reference SBS (#222).
- Storing the genome build in `musica` (#225).
- Smoothing or combining replicates. A single replicate is shipped, and
  rep2 agrees with it on about 80% of bins.
