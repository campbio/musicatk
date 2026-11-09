# Spec: rebuild the example data from GDC hg38 samples

Fixes #238 and #217. Branched from `fix/import-ind83-followups` (PR #243,
now merged), because the rebuild needs its import fixes.

## Problem

- **#238.** Six PASS insertions in two bundled VCFs have a REF base that
  isn't the hg38 base at POS. The VCFs were converted from older GDC
  MuTect2 MAFs with `maf2vcf`, which took the anchor base from a FASTA
  that didn't match hg38.
- **#217, and the objects built from the example files.** No shipped
  example object has a `data-raw/` script, and all of them predate the
  import and table fixes of #208:
  - `musica` and `musica_annot` (911 variants, 7 samples) were built from
    the bundled extdata files without standardized indels or DBS merging.
    `musica_sbs96` and `musica_sbs96_tiny` are SBS-only versions. `res`
    and `res_annot` are a 3-signature NMF model on that SBS96 table.
  - `dbs_musica` mixes genome builds: 400 SBS from the 3 LUSC samples of
    the bundled MAF (hg38) and 400 DBS from 236 other LUAD samples, which
    match hg19, from an unknown source.
  - `indel_musica` has 300 indels from 5 TCGA UCEC aliquots that match
    hg38, but only 4 of them are in those aliquots' current GDC MAFs, so
    its source can't be recovered.
- **Weak, mislabeled signatures in the vignette.** The bundled files hold
  1 LUAD, 3 LUSC and 3 SKCM samples, about 870 SBS in all, of which about
  500 are in the 4 lung samples. The vignette's 3-signature NMF finds UV
  cleanly (closest COSMIC v3 SBS7b, cosine 0.90), but its two lung
  signatures are weak: closest SBS13 at 0.67 and SBS40 at 0.74 (SBS4
  0.73). The vignette names them `c("Smoking", "APOBEC", "UV")` by
  position, while the order is UV, APOBEC-like, smoking-like, and its text
  ("Signatures 1 and 2 were most highly correlated to COSMIC Signature 7
  and 4") is hard-coded. This was already so on devel: the 6 corrected
  records are insertions, which SBS96 doesn't count.
- `inst/extdata/sample_annotations.txt` labels the SKCM (melanoma)
  samples "Breast".

## Effect on earlier results

Every example object changes, and so do the examples, vignette output and
tests that use them:

- The bundled example files are replaced by 12 new TCGA samples (see
  Approach), so `musica`, `musica_annot`, `musica_sbs96`,
  `musica_sbs96_tiny`, `res` and `res_annot` have new samples, variants,
  counts and signatures. The pattern changes from "two weak lung
  signatures and UV" to three clean signatures (APOBEC, smoking, UV),
  checked before saving.
- `dbs_musica`, `indel_musica`: new samples or variants, all hg38.
- File names in `inst/extdata` change; code that used the old names must
  change too (examples, tests, vignette, articles).
- The vignette's "Reference bases for 6 out of 911 variants did not
  match" warning goes away.

Synthetic and COSMIC objects (`synthetic_breast_*`, `example_predicted_*`,
`full_benchmark_example`, `cosmic_*`, `rep_range*`) and
`luad_tcga_count_table.csv` are unchanged.

## Approach

All variants are hg38, from GDC open-access "Aliquot Ensemble Somatic
Variant Merging and Masking" MAFs (WXS), built with the current
`extract_variants()` and `create_musica_from_variants()` defaults. Every
script runs from the package root with `Rscript`, sets a seed (only
exposure prediction and discovery use random numbers), and saves `.rda`
files with `save(..., compress = "xz")`.

### Choosing the example samples: `data-raw/choose_example_samples.R`

Run as `Rscript data-raw/choose_example_samples.R [download_dir]`. Not
needed to rebuild anything (the build script lists the chosen file IDs);
it shows how the samples were chosen and checks the choice against the
current GDC release.

1. List the open ensemble MAFs of TCGA-LUAD, TCGA-LUSC and TCGA-SKCM from
   the GDC files endpoint (file ID and project), and download the ones not
   already in `download_dir` (about 1,500 files, roughly 350 MB).
2. For each aliquot, count its SBS: SNP records that aren't adjacent to
   another SNP of the same sample (adjacent pairs become DBS in
   `create_musica_from_variants()`).
3. Keep aliquots that:
   - have 150 to 500 SBS: enough for a clean SBS96 profile, few enough to
     keep the files small;
   - aren't whole-genome-amplified DNA (analyte W or X);
   - are primary tumors (sample type 01) for LUAD and LUSC; for SKCM,
     primary (01) or metastatic (06), since most TCGA melanomas are
     metastases;
   - are the first by barcode when a patient has several.
4. For those aliquots, build SBS96 (hg38) and estimate exposures with
   `predict_exposure(..., signature_res = cosmic_v3_sbs_sigs_exome,
   signatures_to_use = c("SBS1", "SBS2", "SBS4", "SBS5", "SBS7a",
   "SBS7b", "SBS13", "SBS40"), algorithm = "lda")`. Convert to shares of
   each sample's total:
   - smoking = SBS4,
   - APOBEC = SBS2 + SBS13,
   - UV = SBS7a + SBS7b.
5. Pick four per group, ranked by the group's share, ties by barcode:
   - **Smoking**: LUAD or LUSC, APOBEC share below 0.1.
   - **APOBEC**: LUAD or LUSC, smoking share below 0.1.
   - **UV**: SKCM.

   Stop if any pick has a group share below 0.5. Print the 12 with file
   ID, project, SBS count and the three shares.

### Writing the example files: `data-raw/example_files.R`

Downloads the 12 MAFs by file ID (hard-coded, with group and project),
checks each holds exactly its aliquot, and writes `inst/extdata`:

- `public_<project>_<aliquot>.vcf`: the smoking lung sample with the
  highest smoking share, for the single-VCF examples.
- `public_TCGA_lung.maf`: the other 7 lung samples. The first line of the
  GDC header (`#version gdc-1.0.0`), then the 34 standard TCGA MAF columns
  plus `t_depth`, `t_ref_count`, `t_alt_count` and `n_depth`; the current
  MAF keeps all GDC columns (about 1.4 KB per record), which would make a
  7-sample file several MB.
- `public_SKCM_<aliquot>.vcf`: one per UV sample (4).
- `sample_annotations.txt`: `Sample_Names` and `Tumor_Subtypes` (`Lung`
  or `Skin`) for the 12 samples, in file order.

The VCFs are VCFv4.2 with `##reference=GRCh38`, contig lines for chr1-22,
X, Y and M with hg38 lengths from the BSgenome, FILTER `PASS`, INFO `.`,
FORMAT `GT:AD:DP`, then a tumor column (`0/1:<t_ref_count>,<t_alt_count>:
<t_depth>`) and a normal column (`0/0:.:<n_depth>`), named by
`Tumor_Sample_Barcode` and `Matched_Norm_Sample_Barcode`. Substitutions
use the MAF alleles at `Start_Position`. Indels use the VCF convention:
POS is the base before the event (the MAF `Start_Position` for an
insertion, `Start_Position - 1` for a deletion), REF and ALT start with
that hg38 base. Records are sorted by chromosome order and position.

The old files (`public_LUAD_TCGA-97-7938.vcf`, the three
`public_SKCM_*.vcf`, `public_TCGA.LUSC.maf`) are deleted. The script then
checks that `extract_variants()` on the new files gives as many variants
as the MAF records, with no reference-mismatch warning on hg38, and that
`inst/extdata` is no larger than now (about 0.9 MB).

### Rebuilding the objects

- **`data-raw/example_musica.R`** builds, from the new files:
  - `musica`: all 12 samples, default settings, SBS96 table.
  - `musica_annot`: `musica` with `sample_annotations.txt`.
  - `musica_sbs96`: the SBS of `musica` with its SBS96 table.
  - `musica_sbs96_tiny`: the first 40 SBS by position of one smoking and
    one APOBEC sample, SBS96 table.
  - `res`, `res_annot`: `discover_signatures(..., "SBS96",
    num_signatures = 3, algorithm = "nmf", seed = 1, nstart = 10)` on
    `musica` and `musica_annot`, model IDs `"res"` and `"res_annot"`.

  Signatures are ordered APOBEC, smoking, UV (as the `name_signatures()`
  example expects). The script stops unless:
  - each signature's closest match in `cosmic_v3_sbs_sigs_exome` is in
    its group (SBS2 or SBS13; SBS4; SBS7a or SBS7b) with cosine at least
    0.8;
  - each sample's largest exposure is its group's signature (groups are
    listed in the script, from `choose_example_samples.R`).

  If the chosen samples can't pass this, stop and bring the numbers back
  for a decision; don't loosen the check.
- **`data-raw/gdc_example_musica.R`** is unchanged: `dbs_musica` (five
  TCGA-LUAD aliquots with the most DBS, SBS96 and DBS78) and
  `indel_musica` (the indels of the earlier five UCEC aliquots).
- `data-raw/choose_dbs_aliquots.R` is unchanged.

### Vignette

- Import section: the new file names; `luad_vcf` becomes the single lung
  VCF, `lusc_maf` becomes the lung MAF, `melanoma_vcfs` the four SKCM
  VCFs (variable names follow the files).
- Signature names are assigned from the data, not by position:

  ```r
  matches <- compare_cosmic_v3(musica, "ex_result", "exome",
                               threshold = 0.5)
  best <- matches[!duplicated(matches$x_sig_name), ]
  etiology <- c(SBS4 = "Smoking", SBS2 = "APOBEC", SBS13 = "APOBEC",
                SBS7a = "UV", SBS7b = "UV")
  ```

  Each signature gets the etiology of its best match, or that match's
  COSMIC name when it isn't in `etiology`, and the result goes to
  `name_signatures()`.
- The paragraph after `compare_cosmic_v2()` reports the matches with
  inline R from that call's result, instead of hard-coded signature
  numbers.
- `plot_sample_counts()` shows one sample of each group, by name.
- Other text that quotes numbers or sample types from the data is
  checked and updated.

### Other references to the old files

- `R/load_data.R` examples, tests, and
  `vignettes/articles/tutorial_tcga_cnsl.Rmd`: new file names.
- Commented-out code in `inst/shiny/server.R` that names the old files is
  left alone.

### Docs

- `R/test_data.R`: descriptions of `musica`, `musica_annot`,
  `musica_sbs96`, `musica_sbs96_tiny`, `res` and `res_annot` updated (12
  samples, three groups), `@source` naming the scripts. `make docs`.
- `data-raw/repair_sample_alignment.R` keeps its note that the objects it
  repaired are now rebuilt.

## Tests, examples, vignette

- Update tests that depend on exact values of these objects or files;
  never weaken an assertion that checks behaviour, only expected numbers
  that come from the data.
- `test-load-vcf.R`: "Testing VCF Input" reads the single lung VCF and
  expects its record count; "Reference bases of the example files match
  hg38" reads all new files and expects no warning.
- Run all examples (`make check-full`) and the vignette; the vignette must
  show three signatures labeled Smoking, APOBEC and UV, with the labels
  matching their COSMIC matches.

## NEWS, version

NEWS under 2.7.11: the new example files and samples (12 TCGA samples
chosen for strong smoking, APOBEC and UV signatures), the "Skin" labels,
the signature labels in the vignette, and the rebuilt objects with their
sizes and sources. DESCRIPTION stays 2.7.11. Before hand-off: `make test`,
`make coverage`, `make lint` on changed files, `make check-full`,
`make bioccheck`, and the size of `inst/extdata` and the tarball.
