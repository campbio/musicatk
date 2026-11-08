# Spec: variant import and IND83 follow-ups

Follow-ups from PR B (#216) and PR D (#237) of #208, and from the
SigProfilerMatrixGenerator validation (#218). Fixes #214, #241, #215,
#213, #233, #234 and #236. The example data (#238, #217) is a separate PR.

## Problem

- **#214.** `extract_variants_from_vcf()` with `filter = TRUE` keeps only
  `FILTER == "PASS"`. FILTER `.` (no filters applied; `readVcf()` reads it
  as `NA`) is dropped, so VCFs that mark passing calls with `.`, such as
  the PCAWG consensus VCFs, return nothing.
- **#241.** `extract_variants_from_maf_file()` calls `maftools::read.maf()`,
  which stops with "No non-synonymous mutations found" on a MAF with no
  records or only silent and non-coding records. In `extract_variants()`
  one such file stops the whole call.
- **#215.** ALT `*` (spanning deletion) and symbolic alleles (`<DEL>`,
  breakends) are kept as alleles; `add_variant_type()` calls `C>*` an SBS.
- **#213.** `VariantAnnotation::expand()` splits multiallelic sites without
  trimming, so REF `ATG`, ALT `A,AT` gives `ATG>AT`. Standardizing strips
  one anchor base only, so the row is removed as a compound deletion.
  The same cause makes REF `AT`, ALT `A,GT` give the DBS `AT>GT`, which is
  the SBS A>G.
- **#233.** `.classify_indels()` requests flanks beyond the chromosome
  ends, and `getSeq()` stops the whole IND83 table.
- **#234.** `.classify_indels()` assumes the `-` form. With
  `standardize_indels = FALSE` the anchor base is counted as part of the
  indel (chr1:153056355 `C>CT` is `INS_repeats_2_0`, not `INS_T_1_5+`).
- **#236.** Labels outside the 83 motifs (for example `INS_N_1_0`) become
  `NA` in `factor()` and `table()` drops them without a warning.

## Effect on earlier results

Counts change only for inputs that have the cases above: VCFs with
FILTER `.`, silent-only or empty MAFs, `*` or symbolic ALT, multiallelic
or otherwise non-minimal VCF alleles, indels near chromosome ends. The
bundled example files have none of them (no multiallelic, `.`, `*` or
non-minimal records; the LUSC MAF has non-synonymous records), and
`data/` objects are not rebuilt, so example counts don't change. IND83
from unstandardized indels now stops instead of returning wrong counts.

## Approach

All in `R/load_data.R` and `R/standard_tables.R`.

- **#214.** With `filter = TRUE`, keep rows whose FILTER is `PASS`, `.`
  or `NA`. When none pass, the warning also lists the FILTER values seen.
  `filter` docs in `extract_variants()`, `extract_variants_from_vcf()` and
  `extract_variants_from_vcf_file()` say `.` counts as passing.
- **#241.** In `extract_variants_from_maf_file()`:
  - Read the file with `data.table::fread()` using the arguments
    `read.maf()` itself uses (`skip = "Hugo_Symbol"`, `quote = ""`,
    `fill = TRUE`), and pass the table to `read.maf()`, so the file is
    parsed once and exactly as before.
  - A MAF with no records warns, naming the file, and returns `NULL`, as
    VCFs do.
  - If `read.maf()` stops with "No non-synonymous", call it again on the
    same table with `vc_nonSyn` set to the classes in the file, so
    silent-only files (and files with class names `read.maf()` doesn't
    know) are read. Normal files give the same rows in the same order;
    any other error is raised unchanged.
  - `@return` says `NULL` is returned for a MAF without records.
- **#215.** `extract_variants_from_vcf()` leaves out alleles whose ALT has
  characters other than A, C, G, T and N (`*`, `<DEL>`, breakends, `.`),
  with a message giving the count. If nothing is left, the existing "No
  variants were found" warning and `NULL`.
- **#213.** After expansion and genotype selection,
  `extract_variants_from_vcf()` trims each row's alleles to their minimal
  form with an internal `.trim_alleles()`: first the shared suffix
  (keeping at least one base of each allele), then the shared prefix,
  keeping one anchor base when the lengths differ and none when they are
  equal. `start` moves by the trimmed prefix and `end` is
  `start + nchar(ref) - 1`. Alleles are uppercased first. Only rows with
  both alleles longer than one base can change. So `ATG>AT` becomes
  `TG>T` one base on, a simple deletion; `AT>GT` becomes the SBS `A>G`; `ATG>AG` becomes `AT>A`.
  Complex alleles with different first bases (`A>GTT`) are unchanged and
  are still removed as complex later. VCF input only.
- **#233.** `.classify_indels()` clamps flanks to `1 .. seqlengths(g)`.
  A flank that falls wholly outside the chromosome is empty. Shorter
  flanks just count fewer repeats or less microhomology. Lengths are
  looked up by chromosome name (`chr` can be a factor), and indels on
  chromosomes missing from `g` stop with an error naming them.
- **#234.** `create_ind83_table()` stops when an insertion's `ref` or a
  deletion's `alt` isn't `-`, saying how many and to build the musica
  object with `standardize_indels = TRUE`. The `standardize_indels` docs
  say IND83 needs it.
- **#236.** `create_ind83_table()` calls `.warn_unmatched_motifs()` on the
  labels before counting, as SBS96 does, with a hint that fits IND83
  (an N in the inserted or deleted bases).

## Tests

One regression test per issue, each failing on its assertion on current
code. Tiny fixtures from `write_test_vcf()` and temporary MAFs.

- #214 (`test-load-vcf.R`): `PASS`, `.` and `LowQual` rows keep 2 rows;
  only failing rows warn naming `LowQual` and return `NULL`.
- #241 (`test-load-variants.R`): a header-only MAF warns and returns
  `NULL`; a silent-only MAF (3 Silent rows of the LUSC MAF) returns 3
  rows; `extract_variants()` on (empty MAF, LUSC MAF) returns the LUSC
  rows.
- #215 (`test-load-vcf.R`): a site with ALT `A,*` and a `<DEL>` row keep
  only the `A` allele, with a message.
- #213 (`test-load-vcf.R`): REF `ATG`, ALT `A,AT` with GT `1/2` gives
  deletions `TG>T` at start+1 and `ATG>A` at start; REF `AT`, ALT `A,GT`
  with GT `1/2` gives `AT>A` and the SBS `A>G`; REF `ATG`, ALT `AG` gives
  `AT>A`. Through `create_musica_from_variants()` (no ref check) none is
  removed as compound and the types are DEL, DEL, DEL, SBS, DEL.
- #233 (`test-tables.R`): insertions at hg38 chrM:2 and chrM:3 get
  `INS_T_1_1` and `INS_T_1_0`, and a deletion of the last base of chrM is
  counted. A factor `chr` gets the same bins, and a chromosome missing
  from `g` stops with an error naming it.
- #234 (`test-tables.R`): chr1:153056355 `C>CT` with
  `standardize_indels = FALSE` stops with an error naming
  `standardize_indels`.
- #236 (`test-tables.R`): an `N` insertion warns "not in the IND83
  schema".

The existing tests must still pass.

## Docs, NEWS, version

Roxygen changes above, then `make docs`. NEWS under 2.7.10, one entry
per issue, saying which inputs change. DESCRIPTION 2.7.10. Devel only, as
with PRs B and D. PR says `Fixes` for each issue. Before hand-off:
`make test`, `make coverage`, `make lint` on changed files,
`make check-full`, `make bioccheck`.
