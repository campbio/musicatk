# musicatk 2.7.7 (2026-11-06)
* `annotate_replication_strand()` gives A>* and G>* SBS the replication strand of their pyrimidine, the strand SBS192 counts them on. `rep_range` labels the plus strand, so these SBS now take the opposite label: a G>T in a lagging range is a C>A on the leading-strand template. Before, they kept the plus-strand label, so they were counted on the wrong strand in SBS192_Rep tables. SBS192_Rep counts and the `Replication_Strand` table change for every SBS with an A or G reference base; SBS with a C or T reference base are unchanged (#222).

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
