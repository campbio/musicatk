# Spec: IND83 counting (PR D of #208)

## Problem

`create_ind83_table()` in `R/standard_tables.R` crashes on common inputs
and puts some indels in the wrong ID83 bin:

- **#126.** It crashes unless every indel class is present. Input with
  only deletions or only insertions stops in `subset_variant_by_type()`
  ("No variants of type: INS"). A missing 1 bp deletion, 1 bp insertion
  or ≥2 bp insertion class leads to `rbind()` of an undefined `del1_ta`,
  `ins1_ta` or `ins2_ta`. When ≥2 bp deletions are all in repeat bins or
  all in MH bins, `.count2_del()` calls `paste()` on zero-length vectors,
  gets one string, and `table()` fails ("all arguments must have the same
  length"). Reproduced on `indel_musica` subsets.
- **#127.** `.micro_left()` compares the deletion's suffix with the far
  (5') end of the left flank rather than the bases next to the deletion:
  `.micro_left("ACGTT", "TTTTTTTTTTTTTTTTTTTTTTGTT")` gives 4, not 3.
- **#128.** Left-side repeats reverse the flank but not the unit, so a
  non-palindromic multi-base unit never matches on the left:
  `.count_repeat("AC", stri_reverse("TTACACAC"))` gives 0, not 3.
- **#129.** Flanks are `start - w .. start - 1` and `end + 1 .. end + w`.
  That is right for deletions (since #108 both MAF and VCF deletions span
  the deleted bases), but not for insertions, where `start` is the base
  before the insertion: the left flank skips that base, and for MAF
  insertions (`end = start + 1`) the right flank skips the base after the
  insertion too. Example in `indel_musica`: chr1:153056355 +T (C|TTTTTGCT)
  is counted `INS_T_1_4`; it is `INS_T_1_5+`.
- **#130.** ≥2 bp insertion flanks are a fixed 10 bp, so a 3 bp unit can
  match at most 3 copies on each side and `INS_repeats_4_5+` and
  `INS_repeats_5+_5+` are nearly unreachable.
- **#9** (2020) reported miscounts against SigProfiler; these are the
  causes found. End-to-end validation against SigProfilerMatrixGenerator
  is #218.

## Effect on earlier results

IND83 tables change. On `indel_musica` (300 indels, all from MAF), a
prototype of the fix moves 16 indels to another bin: 12 of the 56 one-bp
insertions (7 of them `INS_T_1_4` → `INS_T_1_5+`), 3 ≥2 bp insertions
(e.g. chr16:71647917 +TC, `INS_repeats_2_0` → `_2_4`) and 1 MH deletion
(chr9:128695510 del AATA, `DEL_MH_4_3` → `_4_2`). Each was checked
against the hg38 sequence. Totals are unchanged. VCF input with insertions or
deletions inside multi-base repeats changes more. The stored
`indel_musica` object has no IND83 table, so nothing in `data/` is stale
from this PR (rebuilding the object itself is #217). SBS and DBS tables
are unchanged.

## Approach

Replace the four counters with one per-variant classifier, then build
the table with one `table()` call.

- **Coordinates (#129).** One helper fetches flanks of width `w` for
  every indel. Deletions of `n` bases: left `start - w .. start - 1`,
  right `start + n .. start + n + w - 1`. Insertions: left
  `start - w + 1 .. start`, right `start + 1 .. start + w`. `end` is not
  used, so MAF (`end = start + 1`) and VCF (`end = start`) insertions get
  the same flanks. One `getSeq()` call per side, as now.
- **Flank width (#130).** `w = 7` for 1 bp indels (as now) and `w = 5n`
  for longer ones, insertions and deletions alike.
- **Repeats (#128).** `.count_repeat()` keeps counting copies at the
  start of the right flank. A new `.count_repeat_left(unit, flank)`
  counts copies at the end of the left flank by reversing both the unit
  and the flank (`stringi::stri_reverse`, as now, so `stringi` stays in
  use; dropping it would be a DESCRIPTION change). Total repeats = left +
  right, as now.
- **Microhomology (#127).** `.micro_left()` compares the last `k` bases
  of the deletion with the last `k` bases of the left flank, for the
  largest `k < n`. `.micro_right()` is already right.
- **Classification.** Unchanged rules: 1 bp indels use the pyrimidine
  base (A/G reverse-complemented with swapped flanks) and cap repeats at
  5+; ≥2 bp insertions → `INS_repeats_<len>_<repeats>`; ≥2 bp deletions
  with 0 repeats and MH > 0 → `DEL_MH_<len>_<mh>`, otherwise
  `DEL_repeats_<len>_<repeats>`; lengths, repeats and MH capped at 5+.
- **Table (#126).** The classifier returns one label per indel; the table
  is `table(factor(label, levels = <83 motifs>), <sample factor>)`, so any
  mix of classes works and every motif and sample has a row and column.
  Insertions and deletions are taken with `Variant_Type %in% c("INS",
  "DEL")`, not `subset_variant_by_type()`. With no indels at all it stops
  with "No insertions or deletions found", as SBS96 does with no SBS.
- Row order, motif names, `annotation`, `features` and `type` slots are
  unchanged.

Out of scope, with follow-up issues: indels within `w` bases of a
chromosome end make `getSeq()` fail (unchanged); IND83 from input built
with `standardize_indels = FALSE` keeps the VCF anchor base in
`ref`/`alt` and is miscounted (unchanged; the docs already say indels are
standardized); the all-`NA` IND83 `features` and the dead
`ifelse(ins, plus <- 0, plus <- 0)` stay for #139 (PR K).

## Tests

One regression test per issue in `tests/testthat/test-tables.R`, each
failing on its assertion on current code:

- #126: IND83 from deletions only, insertions only, no 1 bp deletions,
  and a single ≥2 bp deletion (all from `indel_musica`, hg38): 83 rows,
  one column per sample, total equals the number of indels.
- #127: `.micro_left("ACGTT", "TTTTTTTTTTTTTTTTTTTTTTGTT")` is 3.
- #128: `.count_repeat_left("AC", "TTACACAC")` is 3, and a 2 bp deletion
  in the middle of an AC repeat on hg38 gets the same bin as the
  left-aligned one.
- #129: chr1:153056355 +T as a MAF insertion (`end = start + 1`), as a
  VCF insertion (`end = start`) and shifted one base right inside the
  T run all give `INS_T_1_5+`.
- #130: a 3 bp insertion in a hg38 tract of ≥5 copies of that unit gives
  `INS_repeats_3_5+`.

The existing "IND83 columns follow sample_names()" test must still pass.

## Docs, NEWS, version

A `@details` paragraph on `build_standard_table()` saying how indels are
binned (pyrimidine base for 1 bp, copies of the unit on both sides,
microhomology on either side for deletions without repeats), then
`make docs`. NEWS under 2.7.9 with the `indel_musica` numbers measured
from the final code; DESCRIPTION 2.7.9. PR says `Fixes #126`, `#127`,
`#128`, `#129`, `#130` and `#9`. Before hand-off: `make test`,
`make coverage`, `make lint` on the changed files.
