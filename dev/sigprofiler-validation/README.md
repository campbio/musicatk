# SBS96, DBS78 and IND83 counts against SigProfilerMatrixGenerator

One-off validation for issue #218, run 2026-11-08 after the import fixes
(#216) and the IND83 counting fixes (#237). It checks that musicatk's
standard count tables match SigProfilerMatrixGenerator (SPMG) on the same
variants. Nothing here is part of the package (`dev/` is in
`.Rbuildignore`), and SPMG is not a dependency.

**Result:** every difference has a known cause. Where both tools read a
variant the same way, all SBS96, DBS78 and IND83 counts agree. The
differences come from:

1. **A difference in convention.** SPMG also counts the two bases of each
   DBS in SBS96; musicatk does not.
2. **An SPMG error with GDC MAF input.** SPMG puts each MAF insertion one
   base to the left of where it is. Given the same insertions as VCF records,
   SPMG bins them as musicatk does.
3. **Bad example data.** Six insertions in the bundled LUAD and SKCM VCFs
   have a REF base that doesn't match hg38. SPMG skips them; musicatk warns
   and counts them.

## Versions

| | |
|---|---|
| musicatk | 2.7.9 (`devel` at fb46ef9), loaded from the working tree |
| R | 4.6.1, BSgenome.Hsapiens.UCSC.hg38 1.4.5 |
| SigProfilerMatrixGenerator | 1.4.0, GRCh38 reference from `install("GRCh38")` |
| Python | 3.12 |

## Inputs

The bundled files in `inst/extdata`: the LUSC MAF (3 samples) and the LUAD
and three SKCM VCFs (4 samples), all hg38.

- **musicatk:** `extract_variants()`, `create_musica_from_variants()` and
  `build_standard_table()` for SBS96, DBS78 and IND83, all with default
  arguments.
- **SPMG:** `SigProfilerMatrixGeneratorFunc()` with `"GRCh38"`,
  `exome = False`, `plot = False` and `seqInfo = True`. `seqInfo` writes the
  bin of each variant, which the per-variant indel comparison uses. SPMG runs
  once on the VCFs and once on the MAF.

**Filtering.** SPMG ignores the VCF FILTER column, while musicatk keeps only
PASS records from a VCF by default. To compare counting rather than
filtering, `run_spmg.py` gives SPMG copies of the VCFs with only their PASS
records (8 records removed). Neither tool filters MAF records. All 602 LUSC
records are counted by both, including 28 marked `panel_of_normals`.

**Sample names.** These are matched by TCGA barcode: SPMG names a VCF sample
after its file, and musicatk names it after the tumour column.

## Results

| Table | musicatk | SPMG | Samples identical | Cells different |
|---|---|---|---|---|
| SBS96 | 870 | 884 | 4 of 7 | 12 |
| DBS78 | 7 | 7 | 7 of 7 | 0 |
| IND83 | 27 | 21 | 4 of 7 | 9 |

Both tools report all 7 samples.

### SBS96: SPMG also counts DBS bases as SBS

Both tools have the same DBS78 counts (7 DBS). musicatk removes the 14 bases
of these DBS from the SBS (`convert_dbs = TRUE`); SPMG keeps them in SBS96 as
well. Every one of the 14 extra SPMG counts is a base of a DBS: SPMG's
per-variant output has an SNV row at both positions of each DBS. For example,
`C[CC>TT]A` in TCGA-EE-A3J5 adds `C>T_CCC` and `C>T_CCA`. The other 870 SBS
are counted the same in every motif.

This is a difference in convention, not an error. musicatk's SBS96 and DBS78
don't share any variants; with SPMG, a DBS is counted in both.

### DBS78: identical

Every sample has the same count in every DBS78 motif.

### IND83: 19 of 27 indels agree; the other 8 are explained

All 16 deletions agree, as do 3 of the 11 insertions: 2 from the MAF and 1
from a VCF. The other 8 are all insertions.

**SPMG puts GDC MAF insertions one base to the left (2 insertions).** SPMG
reads every MAF indel as starting at `Start_Position - 1`
(`convert_input_to_simple_files.py`, `convertMAF()`). That is right for
deletions, where `Start_Position` is the first deleted base. In a GDC MAF,
`Start_Position` of an insertion is already the base before the inserted
sequence, so SPMG places the insertion one base too far to the left. musicatk
follows the GDC convention. For two insertions in TCGA-77-7335 this changes
the repeat count:

| Insertion | Flanking hg38 sequence | musicatk | SPMG |
|---|---|---|---|
| chr3:52053983 +CTC | `...TCTG G [CTC] CTCCAG...` | `INS_repeats_3_1` | `INS_repeats_3_0` |
| chr3:52404546 +CT | `...CTGGA A [CT] CTCGG...` | `INS_repeats_2_1` | `INS_repeats_2_0` |

As a check, `run_spmg.py` also gives SPMG the 4 MAF insertions as VCF
records anchored on the base at `Start_Position`, with REF set to that hg38
base. SPMG then bins all four as musicatk does (`INS_repeats_3_1`,
`INS_repeats_2_1`, `INS_T_1_1` and `INS_repeats_5+_0`). So the two tools
classify insertions the same way, and the difference comes only from SPMG's
reading of MAF insertion coordinates.

**Six VCF insertions have a REF that doesn't match hg38 (6 insertions).** In
the bundled VCFs, the REF base of 6 of the 7 insertions isn't the hg38 base
at POS: four in LUAD TCGA-97-7938 and two in SKCM TCGA-ER-A197. For example,
`chr7 48036323 T TA` has `C` at that position in hg38. There is no consistent
offset that would explain the REF bases, and the other variants in these
files match hg38, so where these insertions really are is unknown.
- **SPMG** skips these records ("The reference base does not match the
  reference chromosome position").
- **musicatk** warns ("Reference bases for 6 out of 911 variants did not
  match") and then counts them using POS, so its bins for them can't be
  trusted.

## Not covered by the example data

The bundled files have no run of three or more adjacent SBS, no complex
indel and no multiallelic record. `run_spmg.py` and `compare.R` check the
first two with a small synthetic VCF (hg38 bases):

- **Run of three adjacent SBS** (chr8:8328070-8328072 `GCC>TAA`): SPMG
  reports it as an MNS and still counts all three bases in SBS96. musicatk
  makes it one `"MNV"` variant, which no standard table counts, so it is
  missing from SBS96 and DBS78. That is documented in
  `create_musica_from_variants()`.
- **Complex indel** (chr3:52404546 `AC>GTT`): SPMG counts it as an indel but
  puts it in no ID83 bin. musicatk removes it with "Removing 1 compound
  insertions". Neither tool counts it in ID83; including such indels is #7.

## Rerun

You need SPMG and its GRCh38 reference in a Python environment (Python 3.12
works):
- `pip install SigProfilerMatrixGenerator`
- `from SigProfilerMatrixGenerator import install; install.install("GRCh38")`,
  a download of about 800 MB that unpacks to 2.9 GB.

Then, from the package root:

```sh
python dev/sigprofiler-validation/run_spmg.py inst/extdata <out_dir>
Rscript dev/sigprofiler-validation/compare.R inst/extdata <out_dir>
```

`compare.R` loads musicatk from the working tree with `pkgload`. It prints the
summary, the differing cells and indels, the anchored MAF insertion check and
the synthetic variants, and writes these files to `<out_dir>`:

- `summary.csv`
- `cells_different.csv`: every motif count that differs
- `indels_per_variant.csv`: the IND83 bin of each indel from both tools
- `variants_not_counted.csv`: variants of other types, such as MNV
