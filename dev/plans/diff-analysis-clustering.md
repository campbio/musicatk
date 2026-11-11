# Spec: differential analysis and clustering statistics (PR G of #208)

Decisions by the maintainer on 2026-11-10:

- `cluster_exposure(method = "kmeans")` runs real k-means; a new
  `method = "fanny"` keeps the current fuzzy clustering.
- `plot_differential_analysis()` is redesigned, not deprecated.
- glm.nb columns for the reference group are labelled as the intercept.
- Devel only, no release port.
- wilcox: BH within each comparison; the p-value and FDR columns are
  renamed `p-value` and `fdr`, matching kruskal.
- `dis.method = "cosine"` uses 1 − cosine similarity (found while
  testing #173).

Issues: #169, #170, #171, #173, #174, #185, #34. All reproduced on
`devel` at e07c07a.

## Problem

- `exposure_differential_analysis(method = "glm.nb")` adjusts the wrong
  p-values. Every `(fdr)` column holds the next group's adjusted p, and
  the last holds the ANOVA's. The first group's `(coef)`, `(Std. Error)`,
  `(z)` and `(Pr(>|z|))` columns are the intercept, but they are labelled
  as a group effect (#169).
- `method = "wilcox"`:
  - three or more comparisons crash with "only 0's may be mixed with
    negative subscripts" (#170);
  - p-values and FDRs come back as character strings (#171).
- `plot_differential_analysis()` (#174):
  - for wilcox, it plots the melted text columns;
  - for glm.nb, it plots the intercept as an effect;
  - for anything else, it returns NULL.
- `cluster_exposure(method = "hclust")` passes a distance matrix to
  `hcut(isdiss = FALSE)`, which computes Euclidean distances between its
  rows, so `dis.method` acts only indirectly (#173).
- `cluster_exposure(method = "kmeans")` runs `cluster::fanny` with
  `maxit = 10`, which warns that it hasn't converged. `k_select()` uses
  `stats::kmeans` (#34).
- `cluster_exposure(dis.method = "cosine")` uses
  `philentropy::distance(method = "cosine")`, which returns cosine
  similarity, as a dissimilarity. pam and fanny treat the most similar
  samples as the farthest apart; hclust would too once it uses the matrix.
  The other `dis.method` options are distances.
- `k_select()` with the default `n = 10` fails for every method on ten or
  fewer samples (#185).

## Changes

`R/differential_analysis.R`, `exposure_differential_analysis()`:

- **wilcox** (#170, #171):
  - One row per signature and comparison, with columns `Signature`,
    `Group1`, `Group2`, `p-value` and `fdr`, all numeric where they hold
    numbers. `Pr(>|z|)` and `(fdr)` become `p-value` and `fdr`, matching
    kruskal; the Wilcoxon test has no z statistic.
  - The FDR is BH within each comparison, across signatures. For two
    groups this gives the same numbers as before.
- **glm.nb** (#169):
  - Coefficients are read by name from each fit (`(Intercept)`,
    `annotations<group>`), so an aliased coefficient gives NA instead of
    shifting the columns.
  - The reference (first) group's columns are named `Intercept(coef)`,
    `Intercept(Std. Error)`, `Intercept(z)` and `Intercept(Pr(>|z|))`.
    The other groups keep `<group>(...)`, and their coefficients are log
    fold changes against the reference.
  - `<group>(fdr)` is BH per group, across signatures, for each
    non-reference group. `ANOVA(fdr)` is unchanged.
  - Docs say which group is the reference.
- **Groups** come from the factor levels of the annotation, so the
  header and the model use the same order.

`plot_differential_analysis()` (#174):

- Plots −log10(FDR) per signature as a dodged bar chart, with a dashed
  line at `threshold` (new argument, default 0.05).
- The method is read from the table's columns:
  - wilcox: one bar per comparison, `Group1 vs Group2`;
  - kruskal: one bar per signature;
  - glm.nb: one bar per non-reference group against the reference, plus
    the ANOVA.
- Stops with a clear error for any other table.
- `analysis_type` and `samp_num` are no longer needed. They default to
  `NULL`, and giving them warns that they are ignored.

`R/kmeans.R`:

- `cluster_exposure()`:
  - `method = "kmeans"` runs `stats::kmeans()` on the exposures, using
    `iter.max`. It stops when `dis.method` isn't `"euclidean"`, since
    k-means uses Euclidean distance (#34).
  - New `method = "fanny"` runs the existing `cluster::fanny()` call,
    with `iter.max` and `tol`.
  - `method = "hclust"` uses `hcut(isdiss = TRUE)`, so `dis.method` is
    used (#173).
  - `dis.method = "cosine"` is 1 − cosine similarity, for every method
    that uses the dissimilarity matrix (pam, fanny, hclust).
  - Docs for `method`, `iter.max` and `tol` say which methods use them.
- `k_select()` lowers `n` to the number of samples minus one, with a
  message, and stops when there are fewer than three samples (#185).

Shiny: the cluster tab doesn't add `"fanny"`. It offers only Euclidean
distance for k-means, since `"kmeans"` now stops for other distances
(found while implementing). Its other calls still work.

## Tests

One regression test per issue, each failing on its assertion on `devel`
first. They go in `tests/testthat/test-differentialanalysis.R` and
`tests/testthat/test-kmeans.R`, using `res_annot` (12 samples) with a
three-group annotation where needed.

- **#169:** with three groups, each `<group>(fdr)` equals
  `p.adjust(<group>(Pr(>|z|)), "BH")`, and those p-values equal the
  `annotations<group>` p-values of a direct `MASS::glm.nb` fit. The
  intercept columns are named `Intercept(...)`.
- **#170:** three comparisons give 9 rows, each comparison's FDR is BH of
  its own p-values, and the comparison labels are right.
- **#171:** `p-value` and `fdr` are numeric, and `p-value` equals the
  `matrixTests::row_wilcoxon_twosample()` p-values.
- **#174:**
  - each method's table gives a ggplot whose values are −log10(FDR) of
    that table;
  - a table with other columns errors;
  - `samp_num` warns.
- **#173:** `hclust` clusters equal
  `cutree(hclust(<dis.method distance>, hc.method), k)` for `"cosine"`
  and `"euclidean"`, which differ.
- **cosine:** `pam` with `dis.method = "cosine"` equals `cluster::pam()`
  on 1 − cosine similarity.
- **#34:**
  - `"kmeans"` gives the same clusters as `stats::kmeans()` with the
    same seed, and doesn't warn;
  - `"fanny"` gives the same clusters as `cluster::fanny()`;
  - `"kmeans"` with `dis.method = "cosine"` errors.
- **#185:** `k_select()` on 10 samples with the default `n` gives a plot
  and a message; on 2 samples it errors.

## Docs, NEWS, version

- Roxygen for the changed arguments, with an example for the new plot.
- NEWS for 2.7.13, results-changing entries first.
- DESCRIPTION 2.7.12 → 2.7.13. `make docs`.
- No vignette uses these functions. The TCGA article uses
  `k_select(..., clust.method = "pam")` and
  `cluster_exposure(method = "pam")`, which are unchanged.

## Results change

Yes:

- glm.nb FDRs, which were shifted, and its column names;
- wilcox now returns numbers instead of strings, and its p-value and FDR
  columns are renamed;
- `cluster_exposure()` assignments for `"kmeans"`, now real k-means;
- `cluster_exposure()` assignments for `"hclust"` whenever `dis.method`
  isn't Euclidean;
- `cluster_exposure()` assignments for `"pam"` with
  `dis.method = "cosine"`;
- `plot_differential_analysis()` output.

Shipped data are unaffected. The PR's Scientific correctness section
needs a person.

## Checks

- `make test` and `make coverage` (must not drop).
- `make lint` on changed lines.
- `make check-full` and `make bioccheck` before the PR.
