# SBS96, DBS78 and IND83 counts against SigProfilerMatrixGenerator

One-off validation for issue #218, run 2026-11-08 after the import fixes
(#216) and the IND83 counting fixes (#237). It checks that musicatk's
standard count tables match SigProfilerMatrixGenerator (SPMG) on the same
variants. Nothing here is part of the package (`dev/` is in
`.Rbuildignore`), and SPMG is not a dependency.

The "bundled example files" below are the 7-sample files that
`inst/extdata` held before #217 replaced them with 12 new samples, so
`compare.R` and `run_spmg.py` no longer find them; the results are kept
as a record.

It was run on three sets of variants:

- the bundled example files in `inst/extdata`;
- the GDC open-access MAFs for all 33 TCGA projects (10,549 samples, 2.57
  million records);
- 801 PCAWG consensus indel VCFs from TCGA donors (1.28 million passing
  indels). PCAWG data is controlled access, so only aggregate counts are
  reported.

**Result:** every difference has a known cause. Where both tools read a
variant the same way, all SBS96, DBS78 and IND83 counts agree. In the TCGA
and PCAWG data, every one of the 1.41 million indels musicatk counts gets
the same IND83 bin from SPMG. The differences come from:

1. **Conventions.**
   - SPMG also counts the bases of each DBS in SBS96; musicatk does not.
   - SPMG calls a DBS only when no other SNV is within 5 bp; musicatk pairs
     any two adjacent SBS.
   - SPMG bins complex records such as `GC>T` and `C>TT` as plain indels;
     musicatk removes them (#7).
2. **An SPMG error with GDC MAF input.** SPMG puts each MAF insertion one
   base to the left of where it is. 92% of TCGA insertions get a different
   bin as a result. Given the same insertions as VCF records, SPMG bins them
   as musicatk does. Reported upstream as
   [SigProfilerMatrixGenerator#264](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/issues/264).
3. **Bad example data.** Six insertions in the bundled LUAD and SKCM VCFs
   have a REF base that doesn't match hg38. SPMG skips them; musicatk warns
   and counts them. Fixing the data is #238.

The large runs also found two musicatk usability problems, filed as
separate issues:

- **#241:** `extract_variants()` stops on a MAF with no non-synonymous
  variants. GDC has empty per-aliquot MAFs, and one of them stops the whole
  call.
- **The default `filter = TRUE` keeps only records whose FILTER is `PASS`.**
  PCAWG marks passing records with `.`, so musicatk reads nothing from those
  files unless `filter = FALSE`.

## Versions

| | |
|---|---|
| musicatk | 2.7.9 (`devel` at fb46ef9), loaded from the working tree |
| R | 4.6.1, BSgenome.Hsapiens.UCSC.hg38 1.4.5 |
| SigProfilerMatrixGenerator | 1.4.0, GRCh38 and GRCh37 references from `install()` |
| TCGAbiolinks | 2.40.0; GDC data downloaded 2026-11-08 |
| BSgenome.Hsapiens.UCSC.hg19 | for PCAWG |
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

The bundled VCFs are to be corrected in #238; musicatk's warn-and-keep
behaviour stays as it is.

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

## All TCGA projects (GDC MAFs)

`download_tcga.R` downloads the GDC open-access masked somatic mutation MAFs
(WXS, one per aliquot) for every TCGA project with TCGAbiolinks.
`compare_tcga.R` then runs each project.

- **musicatk:** `extract_variants_from_matrix()` on the MAF records, as in the
  TCGA tutorial article, then the standard tables with default arguments.
  SBS96 is also built with `convert_dbs = FALSE`, which counts every base as
  an SBS, as SPMG does. (`extract_variants()` on the files would stop at the
  first empty MAF; see #241.)
- **SPMG, MAF:** the same MAFs.
- **SPMG, VCF:** the same records written as VCF. Deletions are anchored on
  the base before `Start_Position` and insertions on the base at
  `Start_Position` (the GDC convention), with the anchor base from hg38. This
  shows what SPMG does once the MAF coordinate problem is out of the way.

| Comparison | musicatk | SPMG | Samples identical |
|---|---|---|---|
| SBS96, musicatk with `convert_dbs = FALSE` | 2,432,355 | 2,432,355 | 10,549 of 10,549 |
| SBS96, defaults | 2,392,890 | 2,432,355 | 6,986 |
| DBS78 | 19,602 | 19,592 | 10,541 |
| IND83 vs SPMG on VCF | 136,657 | 137,405 | 9,883 |
| IND83 vs SPMG on MAF | 136,657 | 136,987 | 5,062 |

- **SBS96:** with every base counted, the two tools agree in every sample.
  With defaults, musicatk leaves out the bases of DBS and MNV, as above.
- **DBS78:** 8 samples differ, by 10 DBS in all, and musicatk always has
  more. SPMG treats SNVs within 5 bp of each other as one MNS cluster
  (`mnv = 5` in `MutationMatrixGenerator.py`), so it calls a DBS only when no
  other SNV is within 5 bp. For example, at chr3:31981081, 31981082 and
  31981084 in TCGA-GS-A9TT, musicatk calls a `GC>AT` DBS and keeps the third
  SBS, while SPMG reports one MNS.
- **IND83 against SPMG on VCF:** SPMG gives every one of the 105,711
  deletions and 30,946 insertions musicatk counts the same bin. SPMG's other
  748 indels (642 deletions, 106 insertions) are all complex MAF records,
  where both alleles are present and differ in length (for example
  `GC>T`). SPMG looks only at the allele lengths, so it bins `GC>T` as a
  deletion of `C` and drops the G>T change. musicatk removes such records
  with a message. Records with both alleles longer than one base are left
  out of ID83 by both tools.
- **SPMG on MAF against SPMG on VCF:** the 105,711 deletions are identical.
  Of the 30,948 insertions, 28,558 (92%) are in a different bin. The most
  common changes put homopolymer insertions with 5 or more repeats in the
  0 bin: `1:Ins:T:0` instead of `1:Ins:T:5` (6,638) and `1:Ins:C:0` instead
  of `1:Ins:C:5` (5,473).
- **`Tumor_Seq_Allele1`:** musicatk takes REF from this column of a MAF. In
  all 2.57 million GDC records it equals `Reference_Allele`, so that choice
  makes no difference here.

## PCAWG consensus indels (aggregate counts only)

`compare_pcawg.R` runs on a folder of PCAWG consensus indel VCFs (GRCh37,
without sample columns; here the 801 files for TCGA donors). PCAWG data is
controlled access: the script writes only totals and bin-by-bin counts, with
no sample names or positions, and deletes its working copies.

- **Filtering.** PCAWG marks passing records with FILTER `.`, not `PASS`.
  Both tools read the same copies, which keep records with `PASS` or `.`
  (1,278,326 of 2,939,937), and musicatk reads them with `filter = FALSE`.
- **Inputs.** musicatk adds `chr` to the chromosome names and uses hg19.
  SPMG uses GRCh37. Records on MT would be left out, but there are none.

| | musicatk | SPMG |
|---|---|---|
| Samples | 801 | 801 |
| IND83 / ID83 total | 1,274,460 | 1,275,079 |

462 of the 801 samples have identical counts.

- **The same bin from both tools:** every one of the 388,042 insertions and
  886,418 deletions musicatk counts.
- **SPMG's 619 extra deletions** are exactly the 619 complex records (first
  base differs) that musicatk removes.
- **Left out by both tools:** 3,247 compound deletions.

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

For the large runs, add TCGAbiolinks (for the download) and BSgenome hg19
(for PCAWG), and install SPMG's GRCh37 reference as well:

```sh
Rscript dev/sigprofiler-validation/download_tcga.R <tcga_dir>
Rscript dev/sigprofiler-validation/compare_tcga.R <tcga_dir> <tcga_out> <python>
Rscript dev/sigprofiler-validation/compare_pcawg.R <pcawg_vcf_dir> <pcawg_out> <python>
```

`<python>` is the Python with SPMG. The TCGA download is about 1 GB, and the
comparison takes about 17 minutes. `compare_tcga.R` writes per-project
tables under `<tcga_out>/projects/` and the combined `tcga_*.csv`, and skips
projects it has already done.

`compare.R` loads musicatk from the working tree with `pkgload`. It prints the
summary, the differing cells and indels, the anchored MAF insertion check and
the synthetic variants, and writes these files to `<out_dir>`:

- `summary.csv`
- `cells_different.csv`: every motif count that differs
- `indels_per_variant.csv`: the IND83 bin of each indel from both tools
- `variants_not_counted.csv`: variants of other types, such as MNV
