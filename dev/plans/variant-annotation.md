# Spec: variant annotation correctness, hg38 examples (PR C of #208)

Decisions from the maintainer on 2026-11-06: indel lengths are signed,
`rep_range` stays hg19 and is documented as such (follow-up issue), and the
fixes go to devel only.

## Problem

The variant annotation functions in `R/annotate_variants.R` give wrong
results or crash: purine-reference SBS get the flank from the wrong side,
intergenic SBS are labelled as transcribed-strand, a custom TxDb can't be
used, replication strand crashes when the table doesn't end in SBS, and
every indel length is `NA`. Count tables drop variants whose motif isn't
in the schema without saying so, and the examples and vignette run hg38
data through hg19 resources, so they count only a fraction of variants.

Issues: #104, #105, #117, #118, #119, #120, #135.

While evaluating #120: `rep_range` is in hg19 coordinates (chr17 ends at
the hg19 length; 16 chromosomes run past their hg38 length), but every
bundled variant dataset is hg38. The Shiny app also passes the hg19
`rep_range` whatever genome the user picks.

## Effect on earlier results

- SBS192_Trans counts only SBS inside genes. Before, every intergenic SBS
  was counted as "T"; on `data(musica)` the table totals 807 rather
  than 884. SBS192_Rep and SBS96 are unchanged.
- `add_flank_to_variants()` changes for A>* and G>* SBS (half of all SBS
  before collapsing). Nothing in the package or the app calls it.
- `Variant_Length` is filled for INS/DEL; before, it was `NA`.
- Replication strand values are unchanged where it didn't crash. SBS
  that don't overlap `rep_range` are `NA`, as they already were in
  effect: the string `"NA"` wasn't among the table levels.

## Approach

- **Flank side (#104).** Reverse-change SBS (A>*, G>*) fetch the mirrored
  window and then reverse-complement it: positions
  `start - flank_end .. start - flank_start` for a right flank and
  `end - flank_end .. end - flank_start` for a left flank. One `getSeq`
  call with per-row start and end.
- **Missing transcript strand (#105).** Fill with `NA_character_`, so
  intergenic SBS stay `NA` and `drop_na_variants()` removes them before
  SBS192_Trans is built.
- **Custom TxDb (#117).** Test `methods::is(genome_build, "TxDb")` first
  and use `genes(genome_build)`; otherwise the "19"/"hg19"/"38"/"hg38"
  branches as before. Fix the "provde" typo in the error.
- **Replication strand (#118).** A full-length `NA_character_` vector
  filled at `sbs_index[queryHits]`; read the metadata with
  `S4Vectors::mcols(rep_range)[[1]]`, not `@listData`.
- **Indel length (#119).** For `INS`/`DEL`, length is
  `nchar(alt) - nchar(ref)` with `-` counted as 0, so insertions are
  positive and deletions negative (as the code meant, and as VCF SVLEN).
  The docs say so. MNV keeps `NA`, as now.
- **hg38 examples (#120).** Examples in `R/annotate_variants.R`,
  `R/discovery_prediction.R`, `R/standard_tables.R`, `R/table_utils.R`
  and the vignette use `"hg38"`. The `rep_range` docs and the `rep_range`
  argument of `annotate_replication_strand()` say it is hg19 and must
  match the variants' build; the replication examples get a comment
  saying so. Follow-up: #220 (an hg38 `rep_range`, and the Shiny app
  using `rep_range` for every genome).
- **Dropped motifs (#135).** One helper, `.warn_unmatched_motifs()`, used
  by `create_sbs96_table()` and `create_sbs192_table()`: it warns with
  the number of variants not in the schema and suggests checking that
  `g` is the build the variants were called on.

Out of scope: an SBS in two overlapping genes on opposite strands keeps
the last gene's strand (unchanged, #221); the in-place mutation pattern (#205);
the misplaced "it to the musica" in the `add_flank_to_variants` roxygen
(#124).

## Tests

One regression test per issue, each failing on its assertion on current
code, in a new `tests/testthat/test-annotate-variants.R` (and
`test-tables.R` for #135):

- #104: a G>T and a C>A at known hg38 positions; left and right flanks
  match the pyrimidine-strand bases.
- #105: on `data(musica)`, intergenic SBS have `NA` strand, and the
  SBS192_Trans total equals the number of SBS with a strand.
- #117: a `TxDb` object gives the same strands as `"hg38"`.
- #118: a table whose last rows are not SBS annotates without error;
  non-SBS rows are `NA`.
- #119: INS `-`>`AT` is 2, DEL `ATG`>`-` is -3; SBS 1 and DBS 2.
- #135: an SBS whose ref doesn't match the genome warns for SBS96 and
  SBS192.

## Docs, NEWS, version

roxygen for the changed examples and `rep_range`, then `make docs`.
NEWS under 2.7.5; DESCRIPTION 2.7.5. Render the vignette with
`make article FILTER=musicatk`. Before hand-off: `make test`,
`make coverage`, `make lint` on the changed files.
