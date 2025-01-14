<!-- badges: start -->
[![R-CMD-check](https://github.com/campbio/musicatk/actions/workflows/R-CMD-check.yaml/badge.svg?branch=master)](https://github.com/campbio/musicatk/actions/workflows/R-CMD-check.yaml) <!-- badges: end -->


Mutational Signature Comprehensive Analysis Toolkit
================

# Introduction

A variety of exogenous exposures or endogenous biological processes can
contribute to the overall mutational load observed in human tumors. Many
different mutational patterns, or “mutational signatures”, have been
identified across different tumor types. These signatures can provide a
record of environmental exposure and can give clues about the etiology
of carcinogenesis. The Mutational Signature Comprehensive Analysis
Toolkit (musicatk) has utilities for extracting variants from a variety
of file formats, contains multiple methods for discovery of novel
signatures or prediction of known signatures, as well as many types of
downstream visualizations for exploratory analysis. This package has the
ability to parse and combine multiple motif classes in the mutational
signature discovery or prediction processes. Mutation motifs include
single base substitutions (SBS), double base substitutions (DBS),
insertions (INS) and deletions (DEL).

<img width="901" alt="image" src="https://github.com/user-attachments/assets/7c2d430e-0924-474a-b800-94b119c6d1b1">
<img width="903" alt="image" src="https://github.com/user-attachments/assets/5cc5ec7c-ee69-4625-aae6-dc75d49c1cf4">


# Installation

Currently musicatk can be installed from on Bioconductor using the
following code:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE)){
    install.packages("BiocManager")}
BiocManager::install("musicatk")
```

To install the latest version from Github, use the following code:

``` r
if (!requireNamespace("devtools", quietly=TRUE)){
    install.packages("devtools")}

library(devtools)
install_github("campbio/musicatk")
```

The package can be loaded using the `library` command.

``` r
library(musicatk)
```

# Tutorials

Detailed usage tutorials are available on the [Camplab website](https://camplab.net/musicatk/current/index.html) under the 'Vignettes' drop-down.

