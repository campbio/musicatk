# Spec: variant import correctness (PR B of #208)

Approved by the maintainer on 2026-11-05.

## Problem

Variant import (`R/load_data.R`) mis-attributes multiallelic alleles to
samples, merges adjacent SBS across samples into one DBS, keeps the VCF
anchor base in indel coordinates, rewrites complex variants as indels,
and fails or returns the wrong object on several valid inputs
(compressed files, matrices, factors, lowercase bases, VCFs without GT,
contigs missing from the genome, malformed VCFs). Count tables built from
these variants are wrong or never built.

Issues: #106, #107, #108, #109, #110, #111, #112, #113, #114, #115, #116,
#121, #123. Also removes the duplicate insertion block and
`length(comp_ins > 0)` from #124 (the rest of #124 stays in PR K).

## Effect on earlier results

Only ID83 tables depend on indel coordinates. VCF deletions read their
left flank one base too far left, so a deletion is mis-binned when the
base beyond the anchor extends the repeat or microhomology; with
left-aligned calls this is occasional. MAF insertions and VCF insertions
are mis-read by the ID83 counter itself (#129, PR D). Multi-sample MAFs
(#107) and multi-sample multiallelic VCFs (#106) give wrong SBS/DBS
counts routinely.

## Approach

Prefer functions from packages already imported; no new internal helpers.

- **Multiallelic sites (#106, #110).** For a `CollapsedVCF`, the allele
  numbers each row would get from `VariantAnnotation::expand()` are
  `sequence(lengths(alt(vcf)))`. `expand()` splits the rows, and a row is
  kept only when its allele number is among the selected sample's GT
  allele numbers (`0/2` keeps allele 2 only). This handles haploid and
  partly missing GTs, which `genotypeCodesToNucleotides()` doesn't.
  `"exclude"` drops rows with more than one ALT before expanding; for an
  `ExpandedVCF` it drops rows that share chromosome, position and REF. An
  `ExpandedVCF` passed in directly keeps the current rule (any non-zero GT),
  since its allele numbers can't be recovered; the docs say so. ALT is
  always a character column.
- **No GT (#111).** If `geno(vcf)$GT` is absent, every variant passing
  FILTER is kept, with a message. A VCF with no samples needs `rename`;
  otherwise it stops with a clear error.
- **extra_fields (#121).** INFO columns are taken with `drop = FALSE`;
  `Number=A` fields are already split by `expand()`; other list columns are
  joined with commas (`S4Vectors::unstrsplit()`).
- **fix_vcf_errors (#112).** Finish the existing branch: read the file as
  lines, keep the `##` header and `#CHROM` line, drop body lines with REF or
  ALT `0`, stray `#` lines, or the wrong number of fields, write a temporary
  VCF and read it with `readVcf()`. Warn how many lines were dropped; stop
  if it still can't be read. No hard-coded `skip`.
- **Compressed files (#113).** Strip a trailing `.gz`/`.bgz` before
  `tools::file_ext()`, so `x.vcf.gz` is read as a VCF and `x.maf.gz` as a
  MAF.
- **create_musica_from_variants input (#114, #123).** Accept matrices
  (`data.table::as.data.table()`); convert factor columns to character and
  start/end to integer; uppercase ref and alt. `add_variant_type()` uses
  `as.character()` so factors don't crash `nchar()`.
- **Missing contigs (#115).** After style mapping, variants on chromosomes
  not in `seqnames(genome)` are dropped with a warning naming them.
- **Indels (#108, #109).** The anchor base is stripped only when it matches
  (`substr(alt, 1, 1) == ref` for insertions, `substr(ref, 1, 1) == alt`
  for deletions, the rule VariantAnnotation's `isInsertion()` and
  `isDeletion()` use). Others are complex and removed with a message, as
  compound indels are now. Deletions get `start + 1`, so `start..end` is
  exactly the deleted bases, as in MAF. Insertion coordinates are left
  as they are: in both VCF and MAF, `start` is the base before the
  insertion, and PR D makes the ID83 counter read insertion flanks from
  `start`.
- **DBS and MNV (#107).** Per sample and chromosome,
  `GenomicRanges::reduce()` on the SBS positions finds runs of adjacent
  SBS. A run of two becomes one DBS; a run of three or more becomes one
  row with `Variant_Type` `"MNV"`, counted as a single mutation and not
  included in any standard count table. Other rows keep their order; the
  merged row takes the place of the run's first SBS.
- **create_musica_from_counts (#116).** If `x` has row names, they must
  be the standard motif names; rows are reordered to the standard order by
  name. Any other row names stop with an error. Without row names the
  standard order is assumed, with a message.

## Deviations found in review

- `create_musica_from_counts()` uses one small internal helper,
  `.match_count_rows()`, because the SBS96 and DBS78 branches need the same
  check; there is no other new helper.
- The missing-contig check (#115) runs in `create_musica_from_variants()`
  whether or not `check_ref_chromosomes` is set, since count tables can't
  be built for those variants either.
- Deletions also get `end` set from the deleted bases, so hand-made tables
  with `end == start` are handled.
- Only the six required columns are converted from factors; extra fields
  are left as given.
- `create_musica_from_counts()` converts COSMIC/SigProfiler row names
  (`A[C>A]A`, `AC>CA`) to musicatk motif names before matching. Without
  this, COSMIC-ordered tables that loaded correctly before would be
  rejected. Approved by the maintainer on 2026-11-05.

## Tests

Test first, one regression test per issue, in `tests/testthat/test-load-vcf.R`
and a new `tests/testthat/test-load-variants.R`. Tiny VCF text fixtures are
written to `tempdir()`; the hg38 BSgenome is used only where a genome is
required. Each test must fail on its assertion before the fix.

## Docs, NEWS, version

- Roxygen: `multiallele` behaviour, GT handling, the indel coordinate
  convention, DBS/MNV merging, row names in `create_musica_from_counts()`.
  Fix the "Default" wording only where touched. Then `make docs`.
- Vignette: note MNV handling next to the `convert_dbs` sentence.
- NEWS: one entry per fix under 2.7.4. DESCRIPTION 2.7.3 → 2.7.4.

## Not in scope

- Shipped `dbs_musica`/`indel_musica`: no `data-raw/` script builds them
  from variants, so they aren't rebuilt. Whether they contain cross-sample
  DBS merges is unverified.
- ID83 counter fixes: PR D.

## Follow-up check

Compare SBS96/DBS78/ID83 counts from the example VCF/MAF files against
SigProfilerMatrixGenerator after PR D, to confirm the conventions match.

## Release

Devel only.
