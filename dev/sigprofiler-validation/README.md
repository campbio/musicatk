# SBS96, DBS78 and IND83 counts against SigProfilerMatrixGenerator

One-off validation for issue #218, run 2026-11-08 after the import fixes
(#216) and the IND83 counting fixes (#237). It checks that musicatk's
standard count tables match SigProfilerMatrixGenerator (SPMG) on the same
variants. Nothing here is part of the package (`dev/` is in
`.Rbuildignore`), and SPMG is not a dependency.

**Result:** every difference has a known cause. Where both tools read a
variant the same way, all SBS96, DBS78 and IND83 counts agree. The
differences come from three conventions and one defect in the example data:

1. SPMG also counts the two bases of each DBS in SBS96; musicatk does not.
2. SPMG puts each GDC MAF insertion one base to the left of where it is.
3. SPMG skips variants whose REF doesn't match the genome; musicatk warns and
   keeps them.
4. Six insertions in the bundled SKCM and LUAD VCFs have a REF base that
   doesn't match hg38.

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
- **SPMG:** `SigProfilerMatrixGeneratorFunc(genome = "GRCh38", exome = FALSE)`,
  run once on the VCFs and once on the MAF.

SPMG ignores the VCF FILTER column, and musicatk keeps only PASS records by
default. To compare counting rather than filtering, `run_spmg.py` gives SPMG
copies of the VCFs with only their PASS records (8 records removed). Sample
names are matched by TCGA barcode: SPMG names a VCF sample after its file,
and musicatk after the tumour column.

## Results

| Table | musicatk | SPMG | Samples identical | Cells different |
|---|---|---|---|---|
| SBS96 | 870 | 884 | 4 of 7 | 12 |
| DBS78 | 7 | 7 | 7 of 7 | 0 |
| IND83 | 27 | 21 | 4 of 7 | 9 |

### SBS96: SPMG also counts DBS bases as SBS

Both tools find the same 7 DBS. musicatk removes their 14 bases from the SBS
(`convert_dbs = TRUE`); SPMG keeps them in SBS96 as well. Every one of the 14
extra SPMG counts is a base of a DBS: the per-variant SPMG output has an SNV
row at both positions of each DBS. For example, `C[CC>TT]A` in
TCGA-EE-A3J5 adds `C>T_CCC` and `C>T_CCA`. The other 870 SBS are counted the
same in every motif.

This is a difference in convention, not an error. musicatk's SBS96 and
DBS78 don't share any variants; with SPMG, a DBS is counted in both.

### DBS78: identical

All 7 DBS are in the same bin in both tools.

### IND83: 19 of 27 indels agree; the other 8 are explained

All 16 deletions agree, as do 3 of the 11 insertions: 2 from the MAF and 1
from a VCF. The other 8 insertions:

**SPMG puts GDC MAF insertions one base to the left (2 insertions).** SPMG
reads every MAF indel as starting at `Start_Position - 1`. That is right for
deletions, where `Start_Position` is the first deleted base. In a GDC MAF,
`Start_Position` of an insertion is already the base before the inserted
sequence, so SPMG places the insertion one base too far to the left. musicatk
follows the GDC convention. For two insertions in TCGA-77-7335 this changes
the repeat count:

| Insertion | Flanking hg38 sequence | musicatk | SPMG |
|---|---|---|---|
| chr3:52053983 +CTC | `...TCTG G [CTC] CTCCAG...` | `INS_repeats_3_1` | `INS_repeats_3_0` |
| chr3:52404546 +CT | `...CTGGA A [CT] CTCGG...` | `INS_repeats_2_1` | `INS_repeats_2_0` |

As a check, the 4 MAF insertions were given to SPMG as VCF records anchored
on the base at `Start_Position` (REF taken from hg38). SPMG then gave all
four the same bin as musicatk (`3:Ins:R:1`, `2:Ins:R:1`, `1:Ins:T:1`,
`5:Ins:R:0`). So the two tools classify insertions the same way and differ
only in how they read MAF insertion coordinates.

**Six VCF insertions have a REF that doesn't match hg38 (6 insertions).** In
the bundled VCFs, the REF base of 6 of the 7 insertions isn't the hg38 base
at POS: four in LUAD TCGA-97-7938 and two in SKCM TCGA-ER-A197. For example,
`chr7 48036323 T TA` has `C` at that position in hg38. The SNVs and the one
deletion in these files match hg38, so it is the insertions that are wrong,
probably from converting MAFs to VCF (the SKCM `##reference` lines point to
`output_maf2vcf/hg38.fasta`). SPMG skips these
records ("The reference base does not match the reference chromosome
position"). musicatk warns ("Reference bases for 6 out of 911 variants did
not match"), then counts them using POS. Their musicatk bins can't be trusted,
because the true position of each insertion isn't known.

## Not covered by this data

The bundled files have no run of three or more adjacent SBS, no complex
indel and no multiallelic record. A small synthetic VCF (hg38 bases) shows
how SPMG treats the first two:

- **Run of three adjacent SBS** (chr8:8328070-8328072 `GCC>TAA`): SPMG
  reports it as an MNS and still counts all three bases in SBS96. musicatk
  makes it one `"MNV"` variant, which no standard table counts, so it is
  missing from SBS96 (#107).
- **Complex indel** (chr3:52404546 `AC>GTT`): SPMG counts it as an indel but
  puts it in no ID83 bin. musicatk removes it with a message (#109). Neither
  tool counts it; including such indels is #7, and #213 covers multiallelic
  indels that are removed as complex.

## Rerun

From the package root, with SPMG and its GRCh38 reference installed in a
Python environment (`pip install SigProfilerMatrixGenerator`, then
`from SigProfilerMatrixGenerator import install; install.install("GRCh38")`,
a download of about 800 MB that unpacks to 2.9 GB):

```sh
python dev/sigprofiler-validation/run_spmg.py inst/extdata <out_dir>
Rscript dev/sigprofiler-validation/compare.R inst/extdata <out_dir>
```

`compare.R` loads musicatk from the working tree with `pkgload`, prints the
summary, and writes these files to `<out_dir>`:

- `summary.csv`
- `cells_different.csv`: every motif count that differs
- `indels_per_variant.csv`: the IND83 bin of each indel from both tools
- `variants_not_counted.csv`: variants of other types, such as MNV
