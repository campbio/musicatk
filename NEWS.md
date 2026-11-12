# musicatk 2.7.16 (2026-11-12)
* Added the COSMIC v3.6 (May 2026) signatures as new data objects, following ADR 0003: `cosmic_v3_6_sbs_sigs_grch37`, `cosmic_v3_6_sbs_sigs_grch38`, `cosmic_v3_6_dbs_sigs_grch37`, `cosmic_v3_6_dbs_sigs_grch38` and `cosmic_v3_6_indel_sigs`. v3.6 has 101 SBS, 22 DBS and 25 Indel signatures, and revises the probabilities of signatures that were in v3.0. COSMIC publishes one Indel catalogue per version, so the Indel object has no build. `data-raw/make_cosmic_v3_6.R` rebuilds them from the COSMIC files.
* `cosmic_v3_sbs_sigs`, `cosmic_v3_dbs_sigs`, `cosmic_v3_indel_sigs` and `cosmic_v3_sbs_sigs_exome` are unchanged and stay COSMIC v3.0 (GRCh37), so code that uses them, for example in `predict_exposure()`, gives the same results. Their documentation now states the version and build.
* `compare_cosmic_v3()` gains `genome_build` (`"hg19"`/`"GRCh37"` or `"hg38"`/`"GRCh38"`, default `"hg19"`) and `cosmic_version` (`"latest"`, `"3.0"` or `"3.6"`, default `"latest"`). `"latest"` uses the newest version with a catalogue for the sample type, modality and build: v3.6 for genome data and v3.0 for exome data, since COSMIC stopped publishing exome catalogues after v3.0. **With the defaults, genome comparisons now use v3.6 instead of v3.0, so their matches and similarities change**; `cosmic_version = "3.0"` gives the earlier results. Each call prints a message naming the version, modality, sample type and build used. A combination musicatk doesn't ship, such as `cosmic_version = "3.0"` with `genome_build = "hg38"`, or exome with `"hg38"`, is an error that lists the catalogues available.
* In the Shiny app, the prediction tab offers the COSMIC v3.6 sets for each build and the v3.0 sets, labelled with version and build, and defaults to v3.6 SBS for the build the variants were imported with. The comparison tab has a genome build choice, set from the import genome, and a COSMIC version choice. Signatures without a hard-coded aetiology label show their name instead of "NULL", and SBS22a/b and SBS40a/b/c have labels.

# musicatk 2.7.15 (2026-11-11)
* `plot_heatmap(scale = TRUE)`, the default, z-scores each signature across samples. Before, it z-scored each sample across signatures, so default heatmaps look different; a signature with the same value in every sample is shown as 0 (#176).
* `plot_heatmap(subset_tumor =)` works. It keeps samples whose value in any sample annotation column exactly matches one of the given values (several are allowed), stops naming values found in no column or matching no sample in the model, and shows the `annotation` columns, or the columns the values were found in. Before, it always stopped with an error about the number of observations in the top annotation (#175).
* `plot_signatures()` colours the motif header boxes by name from `color_mapping`, matching the bars, and uses white text on dark boxes. Before, header colours were taken in the order of `color_mapping`, so a reordered mapping coloured headers and bars differently. An unnamed `color_mapping` now works, used in the order the values first appear along the motifs (before, it stopped with "arguments imply differing number of rows"); one missing a value stops with an error naming it (#177).
* `plot_signatures(color_variable =)` uses the named column of the count table annotation, and a vector with one value per motif is matched to the motifs. Before, a column name was ignored and a vector was recycled across signatures. Without `color_mapping`, colours are generated when the table's default mapping doesn't cover the values (#178).
* `plot_signatures()` header boxes follow runs of equal values along the motifs, so a `color_variable` whose values recur in several places gets one box per run instead of overlapping boxes (#177, #178).
* `plot_signatures()` no longer adds a global variable named `annotation_text` to the plot (#186).
* `plot_exposures(plot_type = "box")` and `"violin"` with the default `group_by = "none"` give one box or violin per signature. Before, they silently gave a bar plot (#179).
* `plot_exposures(threshold =)` sets exposures below the threshold to 0, as documented. Before, it dropped them, which removed points from box and violin plots (#184).
* `plot_umap()` colour scales include the largest exposure; before, rounding could leave the top samples grey (#180). Annotation labels no longer stretch the y-axis (#181). `color_by = "cluster"` needs `clust`, not `annotation` (#182). Before `create_umap()`, it says to run `create_umap()` first (#183).
* `create_umap()`, and with `proportional = TRUE` `cluster_exposure()`, `k_select()`, `plot_cluster(group = "signature")`, `plot_umap()`, `plot_heatmap()` and `plot_exposures()`, stop with an error naming samples whose exposures are all 0. Before, they failed with unclear errors or plotted missing values (#183).
* `plot_cluster()` no longer uses `aes_string()`, deprecated in ggplot2 (#95).

# musicatk 2.7.14 (2026-11-11)
* `build_standard_table(..., "DBS78")` works when no DBS has a reference that needs reverse-complementing (for example, only CC>TT); before, it stopped with "Must be character or character vector" (#132, #52). `rc()` returns an empty vector for an empty input, and is vectorised, so it is much faster on long vectors.
* `build_standard_table(..., "SBS192")` without `strand_type` stops with an error saying it is needed; before, it failed with "argument is of length zero" (#136). The check for an existing table uses the table's name, `SBS192_Trans` or `SBS192_Rep`, as for the other tables: building it again stops with "already exists" instead of "Table names must be unique", and with `overwrite = TRUE` it warns that the table is overwritten, which it didn't before.
* `build_standard_table()` messages name the table they create, for example "Building count table 'SBS96' ..." for `modality = "SNV96"`, so the name to use afterwards is clear (#17).
* `combine_count_tables()` names the tables in `to_comb` that don't exist; before, that error failed with 'argument "y" is missing' (#137).

# musicatk 2.7.13 (2026-11-10)
* `exposure_differential_analysis(method = "glm.nb")` adjusted the wrong p-values: each `<group>(fdr)` column held the next group's adjusted p-values, and the last held the ANOVA's. Each non-reference group's p-values are now BH-adjusted across signatures. The reference group (the first level of the annotation) has no group effect: its columns, which hold the intercept, are named `Intercept(coef)`, `Intercept(Std. Error)`, `Intercept(z)` and `Intercept(Pr(>|z|))`, and it has no `(fdr)` column. The other groups' `(coef)` is their log fold change against the reference. Coefficients are read by name, so an aliased group gives NA instead of shifting the columns (#169).
* `exposure_differential_analysis(method = "wilcox")` works with any number of comparisons in `group1` and `group2`; three or more failed with "only 0's may be mixed with negative subscripts" (#170). The p-values and FDRs are numbers instead of character strings, in columns renamed `p-value` and `fdr` (from `Pr(>|z|)` and `(fdr)`), as for `method = "kruskal"`. The FDR is adjusted within each comparison, across signatures. For two groups (one comparison) this gives the same values as before; with two comparisons, the FDR was adjusted across both together (#171).
* `cluster_exposure(method = "kmeans")` runs k-means with `stats::kmeans()`, which `k_select()` already used, and needs `dis.method = "euclidean"`. It uses 10 random starts (new `nstart` argument) seeded by a new `seed` argument, default `1`, applied with `withr::with_seed()`, so the clusters are the same each time and the caller's random number state is left unchanged; `seed = NULL` uses it. Before, it ran fuzzy clustering (`cluster::fanny()`) with at most 10 iterations, which often warned that it hadn't converged; that is now `method = "fanny"`, with up to 500 iterations by default (`iter.max` now defaults to `NULL`: 10 for `"kmeans"` and `"hkmeans"`, 500 for `"fanny"`). The Shiny app offers only Euclidean distance for k-means (#34).
* `k_select()` takes a `seed` argument, default `1`, so its plot is the same each time, and runs k-means with 10 random starts.
* `cluster_exposure(method = "hclust")` clusters the dissimilarity matrix that `dis.method` chooses. Before, it clustered the Euclidean distances between rows of that matrix, so `dis.method` acted only indirectly and the clusters differ from those of the chosen distance (#173).
* `cluster_exposure(dis.method = "cosine")` uses 1 - cosine similarity as the distance. Before, `"pam"` and `"kmeans"` used the cosine similarity itself, treating the most similar samples as the farthest apart.
* `plot_differential_analysis()` plots -log10(FDR) for each signature, for all three methods, with a dashed line at its new `threshold` argument (default 0.05). The method is read from the table's columns. Before, it plotted only the first columns of the table as given, so for `"wilcox"` it plotted text and for `"glm.nb"` the intercept, and for other values of `analysis_type` it returned `NULL`. `analysis_type` and `samp_num` are ignored with a warning. For `"glm.nb"`, bars are labelled `"<group> vs <reference>"`, from a `"reference"` attribute that `exposure_differential_analysis()` now sets; NA FDRs are left out with a message (#174).
* `k_select()` lowers `n` to the number of samples minus one, with a message; with the default `n = 10` it failed on ten or fewer samples. It stops with an error for fewer than three samples, or when `n` isn't a number of at least 2 (#185).

# musicatk 2.7.12 (2026-11-09)
* `predict_exposure()` and `auto_predict_grid()` take a `seed` argument, default `1`, so `algorithm = "lda"` gives the same exposures each time, as `discover_signatures()` does; the caller's random number state is left unchanged, and `seed = NULL` uses it as before. LDA predictions change from earlier, unseeded runs. withr moves from Suggests to Imports for this (ADR 0002) (#152).
* `generate_result_grid()` uses its `annotation` argument. Before, it always read a `Tumor_Type` column: with any other column it stored no models, and with `Tumor_Type` every group's model overwrote the previous group's. Models are now named `<algorithm><k>_<group>` (for example `lda3_Lung`), samples with a missing annotation are left out with a message, and an error in one group names the group (#140).
* `generate_result_grid()` keeps the new fit when a model of the same name already exists, storing it as `lda2.1` and so on. Before, it silently kept the old model and discarded the new fit. `auto_predict_grid()` had the same problem with its per-group models (#141).
* `predict_exposure()` matches signatures to the count table by motif name, and stores the signatures in the table's motif order, as does `combine_predict_grid()` for the models it combines. Before, signatures whose rows were in a different order from the table gave wrong exposures without an error, and signatures of another type (for example DBS78 on an SBS96 table) ran with only a recycling warning; these now stop with an error. `algorithm = "decompTumor2Sig"` stops for tables other than SBS96 (#142).
* Samples with no mutations in the count table are left out of the exposures with a message naming them, by `discover_signatures()`, `predict_exposure()`, `generate_result_grid()` and `auto_predict_grid()`; before, they were left out silently. `generate_result_grid()` and `compare_k_vals()` no longer fail on them ("subscript out of bounds", "NA in probability vector"). `discover_signatures()` needs at least two samples with mutations (#143).
* `auto_predict_grid()` works with a `result_name` other than `"result"` (#145) and with a factor annotation column (#146), and says which thresholds to lower when no signature passes them (#147).
* Repeated model ids are numbered `.1`, `.2`, … `.10`, `.11`; the 12th copy was `..11` (#148).
* `compare_k_vals()` stops with an error for an invalid `error_type` before running any fit; before, it ran every fit and returned NA (#149).
* `predict_exposure(algorithm = "lda")` no longer gives 100 warnings for a single sample. `generate_result_grid()` documents its reconstruction error metric (#152).
* `compare_cosmic_v2()` and the other comparison functions work for a model with one signature (#162).

# musicatk 2.7.11 (2026-11-09)
* The example files in `inst/extdata` are replaced by 12 TCGA samples chosen for strong mutational signatures, so the examples and vignette show clean signatures (#217, #238):
  * four lung tumors with a strong smoking signature (COSMIC SBS4) and four with a strong APOBEC signature (SBS2 and SBS13), and four melanomas with a strong UV signature (SBS7a and SBS7b), each with 150 to 500 SBS, from GDC open-access aliquot ensemble masked MAFs (hg38). `data-raw/choose_example_samples.R` chose them and `data-raw/example_files.R` writes the files.
  * The files are `public_LUSC_TCGA-22-5472-01A-01D-1632-08.vcf`, `public_TCGA_lung.maf` (the other seven lung tumors) and four `public_SKCM_*.vcf` files, 4,233 variants in all. They replace `public_LUAD_TCGA-97-7938.vcf`, `public_TCGA.LUSC.maf` and three `public_SKCM_*.vcf` files (911 variants from 7 samples), whose two lung signatures were weak and in which six insertions had a REF base that wasn't the hg38 base at POS. Code that reads the old file names must use the new ones.
  * VCF REF and ALT match hg38, and indels start with the base before the event, so building a musica object from the example files no longer warns "Reference bases for 6 out of 911 variants did not match".
  * `inst/extdata` is 659 KB instead of 912 KB.
* `inst/extdata/sample_annotations.txt` labels the lung tumors `"Lung"` and the melanomas `"Skin"`; the earlier file labelled its three melanomas `"Breast"`.
* The vignette names each discovered signature after the etiology of its most similar COSMIC V3 signature (Smoking, APOBEC or UV), instead of by position, and its comparison to COSMIC reports the matches it found. The vignette had named the signatures in the wrong order.
* The example datasets are rebuilt with the current import and table code, all from TCGA on hg38, and each has a `data-raw/` script and a `@source` (#217):
  * `musica`, `musica_annot`, `musica_sbs96` and `musica_sbs96_tiny` are built from the new example files by `data-raw/example_musica.R`: 12 samples, 4,014 SBS, 64 DBS, 73 deletions and 15 insertions (three complex deletion-insertions in the MAF are removed on import, as usual). `musica_sbs96_tiny` holds 40 SBS each of a smoking and an APOBEC lung tumor.
  * `res` and `res_annot` are a three-signature NMF fit on that SBS96 table: Signature1 is an APOBEC signature (0.99 cosine to a mix of SBS2 and SBS13), Signature2 a smoking signature (SBS4, 0.86) and Signature3 a UV signature (SBS7a, 0.95), each carrying at least 93% of the exposure of its group's samples.
  * These six objects take 137 KB instead of 45 KB.
  * `dbs_musica` holds all variants of the five TCGA-LUAD tumors with the most DBS in GDC open-access MAFs, leaving out whole-genome-amplified DNA (176 DBS), with SBS96 and DBS78 tables. Before, it mixed 400 SBS on hg38 with 400 DBS from other samples on hg19.
  * `indel_musica` holds the 211 insertions and deletions of the same five TCGA-UCEC tumors as before, from their current GDC MAFs. Its earlier 300 indels could not be traced to a source.

# musicatk 2.7.10 (2026-11-08)
* `extract_variants()` and the VCF readers keep records whose FILTER is `.` (no filters applied) when `filter = TRUE`, as well as `PASS`. Before, VCFs that mark passing calls with `.`, such as the PCAWG consensus VCFs, returned no variants. When no record passes, the warning lists the FILTER values found (#214).
* `extract_variants()` and `extract_variants_from_maf_file()` read MAFs whose records are all silent or non-coding; before, `maftools::read.maf()` stopped with "No non-synonymous mutations found". A MAF with no records gives a warning and no variants, as a VCF does, so one such file no longer stops `extract_variants()` from reading the others (#241).
* VCF alleles are trimmed to their minimal form. Splitting a multiallelic site can leave alleles that share more than one base: REF `ATG` with ALT `A,AT` gave `ATG>AT`, which was removed as a compound deletion, and REF `AT` with ALT `A,GT` gave the DBS `AT>GT`. They are now the deletion `TG>T` and the SBS `A>G`. Counts change only for VCFs with multiallelic sites or other alleles that share more than one base; the bundled example files have none (#213).
* ALT `*` (an allele missing because of an overlapping deletion) and symbolic alleles such as `<DEL>` are left out of VCF input with a message. Before, `*` was counted as an SBS (#215).
* IND83 tables count indels near the start or end of a chromosome, using the flank that is there; before, the whole table failed (#233).
* `build_standard_table(..., "IND83")` stops with an error when indels are not standardized (built with `standardize_indels = FALSE`); before, the VCF anchor base was counted as part of the indel and gave wrong bins (#234).
* IND83 tables warn when indels are left out because their motif isn't in the table, for example an insertion of `N` (#236).

# musicatk 2.7.9 (2026-11-08)
* IND83 tables put some indels in the wrong bin. On `indel_musica`, 16 of 300 indels move to another bin (12 of the 56 1 bp insertions, 3 longer insertions and 1 microhomology deletion); totals are unchanged. Rebuild IND83 tables made with an earlier version (#9):
  * Insertion flanks start next to the insertion. Before, the left flank skipped the base before the insertion, and for MAF insertions the right flank skipped the base after it, so for example chr1:153056355 +T in `indel_musica` was `INS_T_1_4` instead of `INS_T_1_5+` (#129).
  * Repeats of multi-base units are counted on the left of the indel too. Before, a unit like `AC` never matched on the left, so indels with copies of their unit on the left (for example, deleting any copy but the first in a repeat) counted too few copies (#128).
  * Microhomology on the left of a deletion is compared with the bases next to it, not the far end of the flank (#127).
  * Insertions of 2 bp or more look five unit lengths to each side, as deletions do. Before, they looked 10 bp, so `INS_repeats_4_5+` and `INS_repeats_5+_5+` were almost never counted (#130).
* `build_standard_table()` builds IND83 tables when some indel classes are missing, for example with only deletions or only insertions; before, it stopped with an error. With no indels at all, it stops with "No insertions or deletions found." The help page says how indels are binned (#126).

# musicatk 2.7.8 (2026-11-07)
* `annotate_transcript_strand()` leaves SBS in overlapping genes on opposite strands as `NA`, like SBS outside genes, so they are left out of SBS192_Trans tables. Before, they took the strand of whichever gene came last in the overlap list, so their `"T"` or `"U"` label was arbitrary. SBS192_Trans counts drop: on `data(musica)`, from 807 to 719 with hg38 and from 584 to 556 with hg19. Re-run `annotate_transcript_strand()` on objects annotated with an earlier version before building SBS192_Trans (#221).

# musicatk 2.7.7 (2026-11-06)
* `annotate_replication_strand()` gives A>* and G>* SBS the replication strand of their pyrimidine, the strand SBS192 counts them on. `rep_range` labels the plus strand, so these SBS now take the opposite label: a G>T in a lagging range is a C>A on the leading-strand template. Before, they kept the plus-strand label, so they were counted on the wrong strand in SBS192_Rep tables. SBS192_Rep counts and the `Replication_Strand` table change for every SBS with an A or G reference base; SBS with a C or T reference base are unchanged. Re-run `annotate_replication_strand()` on objects annotated with an earlier version before building SBS192_Rep. A GRanges passed as `rep_range` must label the plus strand, as `rep_range` does (#222).

# musicatk 2.7.6 (2026-11-06)
* New `rep_range_hg38`: `rep_range` lifted over to hg38, for annotating replication strand on hg38 variants. Labels are swapped in regions inverted between the builds. The examples and the Shiny app use the object that matches the genome build; for mm9 or mm10 the app now shows an error instead of annotating with human data (#220).
* `annotate_replication_strand()` accepts the genome build for `rep_range`: `"hg19"` (or `"GRCh37"`) uses `rep_range` and `"hg38"` (or `"GRCh38"`) uses `rep_range_hg38`. A GRanges still works as before (#220).
* `rep_range` is tagged with genome `"hg19"` and its chromosome lengths, and is documented with its source (GEO GSM4305465) and what `"leading"` and `"lagging"` mean. `data-raw/rep_range.R` rebuilds both objects (#220).
* `rep_range` bins are 1-based. Before, neighbouring bins shared an edge position, so a variant on an edge could match two bins and take whichever label came last. Replication strand changes only for variants on a bin edge (#220).

# musicatk 2.7.5 (2026-11-06)
* `annotate_transcript_strand()` leaves SBS outside genes as `NA`; before, they got the string `"NA"` and were labelled `"T"`. SBS192_Trans tables now count only SBS in genes, so their counts drop (on `data(musica)`, from 884 to 807) (#105).
* `add_flank_to_variants()` takes the flank of A>* and G>* SBS from the pyrimidine strand; before, a left flank returned the bases 3' of the mutated base and a right flank those 5' of it (#104).
* `annotate_transcript_strand()` accepts a `TxDb` object for `genome_build` (#117).
* `annotate_replication_strand()` no longer fails, or copies one strand to every row, when the variant table doesn't end in SBS. Variants without a replication strand are `NA` rather than the string `"NA"` (#118).
* `annotate_variant_length()` gives indels a length: the number of inserted bases, or minus the number of deleted bases. Before, every indel length was `NA` (#119).
* SBS96 and SBS192 tables warn when variants are left out because their motif isn't in the table. For SBS96 this happens when a reference base doesn't match `g`, for example when `g` is not the genome build the variants were called on; for both tables, when the context has N bases (#135).
* Examples and the vignette use the hg38 genome for the hg38 example data. The `rep_range` docs say it is on hg19 (#120).

# musicatk 2.7.4 (2026-11-05)
* `extract_variants_from_vcf()` with `multiallele = "expand"` gives each sample only the alternate alleles in its GT; before, a sample carrying any alternate allele at a multiallelic site got all of them (#106).
* `multiallele = "exclude"` returns `alt` as a character column and no longer drops every indel from an `ExpandedVCF` (#110).
* VCFs without a GT field (for example from Strelka) keep all variants passing the filter, with a message; before, they returned `NULL`. Such a VCF with several sample columns needs `id` or `rename`, and one without sample columns needs `rename` (#111).
* `extra_fields` keeps the column name of a single INFO field, gives per-allele (`Number=A`) fields the value for each allele, and joins other multi-value fields with commas (#121).
* `extract_variants_from_vcf_file()` with `fix_vcf_errors = TRUE` removes malformed lines and reads the rest, with a warning; before, it returned the function `stats::dt` (#112).
* `extract_variants()` reads `.vcf.gz` and `.maf.gz` files (#113).
* `create_musica_from_variants()` accepts a matrix and factor columns (#114), and reads lowercase bases as uppercase instead of removing them (#123).
* `create_musica_from_variants()` moves VCF-style deletions to start at the first deleted base, as in MAF files, so a deletion gets the same ID83 bin from either format. ID83 counts from VCF input can change (#108).
* Indels whose first base differs between `ref` and `alt` (for example `A > GGT`) are removed as complex variants with a message; before, they were rewritten as indels (#109).
* `convert_dbs = TRUE` merges adjacent SBS only within a sample; before, SBS from different samples could become one DBS, and unsorted neighbours were missed. A run of three or more adjacent SBS becomes one variant of type `"MNV"`, which no standard count table includes. DBS78 and SBS96 counts can change (#107).
* Variants on chromosomes missing from the genome are removed with a warning; before, they crashed the reference check (#115).
* `create_musica_from_counts()` matches rows to motifs by row name; before, it overwrote the row names, so tables in another order were mislabelled. COSMIC/SigProfiler row names (`A[C>A]A`, `AC>CA`) are converted to musicatk motif names; other row names give an error (#116).

# musicatk 2.7.3 (2026-11-05)
* `name_signatures()` returns the renamed `musica` object (invisibly), so `res <- name_signatures(res, ...)` works; before, it returned the vector of names. The object passed in is still renamed in place as before (#154).
* The setters `num_signatures<-`, `other_parameters<-`, `credible_intervals<-`, `model_id<-` and `modality<-` now work; they could never be called before. `model_id<-` and `modality<-` on a `musica` or `result_collection` object also rename or move the stored model, so `get_model()` finds it by the new id or modality (#156).
* `hyperparameter<-` on a `musica` object now assigns the value; before, it always failed (#157).
* `subset_musica_by_annotation()` gives its intended error when the requested values are missing (#163).
* `samp_annot<-` accepts Date and date-time (POSIXct, POSIXlt) vectors, and `variants<-` accepts a `data.frame` (#164).
* `subset_variants_by_samples()` returns all variants of every requested sample; with several names it returned only some of them (#48).

# musicatk 2.7.2 (2026-11-04)
* Samples are now matched by name everywhere. `sample_names()` (the `Samples` column of the sample annotations) sets the sample order; count tables and stored exposures follow it, and new validity methods for `musica`, `result_model` and `result_collection` check this (#153, #160). Objects saved with older versions whose tables are out of order give an error asking to rebuild the table with `overwrite = TRUE`; a `full_benchmark` saved with an older version must be recreated with `create_benchmark()`.
* `exposure_differential_analysis()`, `plot_heatmap()` and `plot_cluster()` pair annotations with exposures by sample name; before, results were wrong when the orders differed, as in the shipped `res_annot` (#153). `plot_cluster(group = "annotation")` now uses the `annotation` argument (#172).
* The shipped musica datasets and `full_benchmark_example` were repaired so their annotations, count tables and exposures share one sample order; no count or exposure values changed. The vignette now matches `sample_annotations.txt` to the samples by name; before, 6 of 7 samples got the wrong `Tumor_Subtypes` label.
* `sample_names()` returns a character vector and no longer partially matches a `Sample` column (#158).
* `create_musica_from_counts()` stores the column names as sample names, so its objects can be annotated; it stops on missing or duplicated column names (#159).
* `subset_musica_by_counts()` and `subset_musica_by_annotation()` also subset stored exposures (#155).
* Count tables no longer get all-zero columns for samples removed by subsetting (#134), IND83 columns are named correctly (#131), and `combine_count_tables()` aligns samples by name (#133).
* `add_result()` matches exposure columns to the samples by name and puts them in sample order; unknown or missing sample names give an error.
* `combine_predict_grid()` returns exposures in sample order and sets `num_signatures` and `model_id`; it stops when the models share samples (#144).
* `create_benchmark()` stores true exposures as signatures x samples and matches them to the count table by sample name; either orientation is accepted (#194).

# musicatk 2.2.1 (2025-07-15)
* Removed conclust dependency
* Fixed plotting issues with indels
* Cosmetic changes to k value and plotting code

# musicatk 2.2.0 (2025-04-15)
* Updated version number to match Bioconductor

# musicatk 2.0.0 (2024-10-30)
* Restructure of the musica class to contain mutational signature discovery and prediction results, thus eliminating the need for the musica_result class.
* Added new functionality to help determine the correct number of signatures to predict (k-value).
* Musica objects can be created from variants or directly from count tables.

# musicatk 1.14.1 (2024-07-13)
* Removed deconstructSigs

# musicatk 1.14.0 (2024-05-03)
* Updated version number to match Bioconductor 
* Fixed some minor bugs and dependencies

# musicatk 1.10.0 (2023-10-16)
* Updated version number to match Bioconductor

# musicatk 1.8.0 (2023-04-11)
* Updated version number to match Bioconductor

# musicatk 1.3.1 (2021-10-10)
* Updated version number to match Bioconductor

# musicatk 1.2.1 (2021-08-08)
* Fixed bug in fonts for some plots
* Added documentation site generated with pkgdown
* Added Shiny UI for interactive analysis of mutational signatures

# musicatk 1.2.0 (2021-08-08)
* Generated new release for Bioconductor 3.13

# musicatk 1.1.2 (2021-02-03)
* Added new functionality for analysis including clustering sample exposures using a variety of methods (e.g. k-means), differential exposure using a variety of methods (e.g., Wilcoxon rank-sum, Kruskal Wallis, GLM), and heatmap creation. Fixed a bug causing incorrect counting of mutations during creation of indel tables.
