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
