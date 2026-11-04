# musicatk 2.7.2 (2026-11-04)
* Samples are now matched by name everywhere. `sample_names()` (the `Samples` column of the sample annotations) sets the sample order; count tables and stored exposures follow it, and new validity methods for `musica`, `result_model` and `result_collection` check this (#153, #160). Objects saved with older versions whose tables are out of order give an error asking to rebuild the table with `overwrite = TRUE`.
* `exposure_differential_analysis()`, `plot_heatmap()` and `plot_cluster()` pair annotations with exposures by sample name; before, results were wrong when the orders differed, as in the shipped `res_annot` (#153). `plot_cluster(group = "annotation")` now uses the `annotation` argument (#172).
* The shipped musica datasets and `full_benchmark_example` were repaired so their annotations, count tables and exposures share one sample order; no count or exposure values changed. The vignette now matches `sample_annotations.txt` to the samples by name; before, 6 of 7 samples got the wrong `Tumor_Subtypes` label.
* `sample_names()` returns a character vector and no longer partially matches a `Sample` column (#158).
* `create_musica_from_counts()` stores the column names as sample names, so its objects can be annotated; it stops on missing or duplicated column names (#159).
* `subset_musica_by_counts()` and `subset_musica_by_annotation()` also subset stored exposures (#155).
* Count tables no longer get all-zero columns for samples removed by subsetting (#134), IND83 columns are named correctly (#131), and `combine_count_tables()` aligns samples by name (#133).
* `combine_predict_grid()` returns exposures in sample order and sets `num_signatures` and `model_id` (#144).
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
