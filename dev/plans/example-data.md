# Spec: rebuild the example data from GDC hg38 samples

Fixes #238 and #217. Branched from `fix/import-ind83-followups` (PR #243),
because the rebuild needs its import fixes; this PR goes in after #243.

## Problem

- **#238.** Six PASS insertions in two bundled VCFs have a REF base that
  isn't the hg38 base at POS. The VCFs were converted from older GDC
  MuTect2 MAFs with `maf2vcf`. Two of the six (LUAD chr13:52024294 +G,
  chr14:93612828 +A) are in today's GDC ensemble MAF for that aliquot at
  the same `Start_Position` with the same inserted base, so only the REF
  base is wrong: `maf2vcf` took the anchor base from a FASTA that didn't
  match hg38. The other four aren't in today's MAFs.
- **#217, and the objects built from the example files.** No shipped
  example object has a `data-raw/` script, and all of them predate the
  import and table fixes of #208:
  - `musica` and `musica_annot` (911 variants, 7 samples) were built from
    the bundled extdata files. The 8 VCF indels were never standardized
    (they are typed `"unknown"`, including the 6 above) and adjacent SBS
    were not merged into DBS, so SBS96 has 884 SBS where the current code
    gives 870 SBS and 7 DBS. `musica_sbs96` and `musica_sbs96_tiny` are
    SBS-only versions. `res` and `res_annot` are a 3-signature NMF model on
    that SBS96 table: two smoking-like signatures in the lung samples and
    one UV signature in the melanoma samples.
  - `dbs_musica` mixes genome builds: 400 SBS from the 3 LUSC samples of
    the bundled MAF (hg38) and 400 DBS from 236 other LUAD samples (TCGA
    and Broad), which match hg19, from an unknown source.
  - `indel_musica` has 300 indels from 5 TCGA UCEC aliquots that match
    hg38, but only 4 of them are in those aliquots' current GDC MAFs, so
    its source can't be recovered.
- `inst/extdata/sample_annotations.txt` labels the three SKCM (melanoma)
  samples "Breast".

## Effect on earlier results

Every example object changes, and so do the examples, vignette output and
tests that use them:

- `musica`, `musica_annot`: indels standardized, the 6 insertions
  corrected, 7 DBS; SBS96 drops from 884 to 870.
- `res`, `res_annot`: a new NMF fit on the new SBS96. Signature and
  exposure values change; the pattern (two smoking-like signatures in the
  lung samples, one UV signature in the melanoma samples) must hold, and
  is checked against COSMIC v3 (SBS4/SBS92-like and SBS7-like) before
  saving.
- `dbs_musica`, `indel_musica`: new samples or variants, all hg38.
- The vignette's "Reference bases for 6 out of 911 variants did not
  match" warning goes away.

Synthetic and COSMIC objects (`synthetic_breast_*`, `example_predicted_*`,
`full_benchmark_example`, `cosmic_*`, `rep_range*`) are unchanged.

## Approach

All variants are hg38, from GDC open-access data, built with the current
`extract_variants()` and `create_musica_from_variants()` defaults.

- **#238 data fix.** In the two VCFs, set REF of the six insertions to the
  hg38 base at POS and ALT to that base plus the inserted sequence. Keep
  POS and the inserted sequence. The two adjacent SKCM insertions stay two
  records, as in the source MAF.
- **Annotations.** Label the SKCM samples "Skin" in
  `inst/extdata/sample_annotations.txt`.
- **`data-raw/example_musica.R`** builds, from the bundled extdata files:
  - `musica`: all seven files, default settings, SBS96 table.
  - `musica_annot`: `musica` with `sample_annotations.txt`.
  - `musica_sbs96`: the SBS of `musica` with its SBS96 table.
  - `musica_sbs96_tiny`: the same two samples as now, first 40 SBS of
    each by position, SBS96 table.
  - `res`, `res_annot`: `discover_signatures(..., "SBS96",
    num_signatures = 3, algorithm = "nmf", model_id = "res")` with a fixed
    seed and `nstart = 10`, on `musica` and `musica_annot`.
  - It stops if the signatures don't show the expected pattern (best
    COSMIC v3 match by cosine for each signature; one UV-like signature
    carrying most of the melanoma exposure, smoking-like signatures
    carrying most of the lung exposure).
- **`data-raw/gdc_example_musica.R`** downloads GDC open-access aliquot
  ensemble masked MAFs (hg38) by fixed file ID and builds:
  - `indel_musica`: the insertions and deletions of the same five UCEC
    aliquots as now (about 211 indels), no count tables.
  - `dbs_musica`: all variants of five TCGA-LUAD aliquots, with SBS96 and
    DBS78 tables. The aliquots are the five with the most DBS among the
    open TCGA-LUAD MAFs, chosen once during development; the script lists
    the file IDs and says how they were chosen. If the object is much
    larger than the current one (about 10 KB), use fewer samples, since
    the tarball is near the Bioconductor size limit.
- **Docs.** `R/test_data.R` gets a description and `@source` for each
  object, saying which script builds it. `make docs`.
- Both scripts set a seed (only discovery uses random numbers), run from
  the package root with `Rscript`, and save with `usethis::use_data()` or
  `save(..., compress = "xz")`, matching the current files.

## Tests, examples, vignette

- Update tests that depend on exact values of these objects; never weaken
  an assertion that checks behaviour, only expected numbers that come from
  the data.
- Run all examples (`make check-full`) and the vignette; the vignette and
  examples must still show the same patterns.
- `test-load-vcf.R` "Testing VCF Input" reads the LUAD VCF and expects
  121 rows; the count doesn't change (records are corrected, not
  removed).

## NEWS, version

NEWS under 2.7.11: the corrected example VCFs, the "Skin" labels, and
the rebuilt objects with their new sizes and sources. DESCRIPTION 2.7.11.
Before hand-off: `make test`, `make coverage`, `make lint` on changed
files, `make check-full`, `make bioccheck`, and the tarball size.
