# Spec: Plotting fixes (PR H of #208)

Changes results: plots only. No counts, exposures or statistics change.
Devel only; no release port proposed (none of these is a wrong numeric
result, and `subset_tumor` has never worked, so no release user relies
on it).

## Problem (all reproduced on devel a2aa073 unless noted)

- **#175 `plot_heatmap(subset_tumor =)`.** Samples are picked with
  `grepl()` over every annotation column, the heatmap annotation with an
  exact `==` on the first matching column, then passed to
  `HeatmapAnnotation(df =)` as a bare vector. It fails even in the
  simplest case: `plot_heatmap(res_annot, "res_annot", subset_tumor =
  "Lung")` stops with "number of observations in top annotation should be
  as same as ncol of the matrix". It also ignores `annotation`, uses only
  the first of several values (the Shiny app passes several), and calls
  `select()` without `all_of()`.
- **#176 `plot_heatmap(scale = TRUE)`** (the default) runs `scale()` on
  the signatures x samples matrix, which z-scores each sample (column
  means are all 0), not each signature.
- **#177 `plot_signatures()` header colours** are given by position
  (`fill = factor(color_mapping)`), while the bars are coloured by name.
  With `color_mapping = rev(default)`, the C>A header is pink and its bars
  blue. The white label for SBS96 is also positional (always the 2nd box).
- **#178 `.pivot_signatures()`**: a vector `color_variable` (one value
  per motif) is recycled over the motif-major long table, so only 22% of
  rows get the right value with 3 signatures; a column name is ignored
  (`annot[, tab@color_variable]` is used instead).
- **#179 `plot_exposures(plot_type = "box" or "violin")`** with the default
  `group_by = "none"` falls through to the bar plot.
- **#180 `plot_umap()`** colour limits are `c(0, max(signif(breaks, 2)))`;
  when rounding goes down (max 4549 gives 4500), the top samples fall
  outside the limits and are drawn grey. Not hit by `res_annot`
  (max = 1); shown by the regression test.
- **#181 `plot_umap(add_annotation_labels = TRUE)`**: the invisible
  centroid layer maps `y = UMAP_1`, stretching the y-axis.
- **#182 `plot_umap(color_by = "cluster", clust = )`** stops asking for
  `annotation`, which the cluster colouring doesn't use.
- **#183 Unclear errors.** `plot_umap()` before `create_umap()` fails in
  `cbind()` ("number of rows of matrices must match"). A sample whose
  exposures are all 0 gives 0/0: `create_umap()` fails in uwot with
  "Missing values found in 'X'"; `cluster_exposure()`, `plot_heatmap()`,
  `plot_exposures()` and `plot_umap()` with proportions get NaN. Since
  #143 the package's own discovery/prediction leaves such samples out,
  so this needs a hand-built or edited result.
- **#184 `plot_exposures(threshold =)`** drops rows at or below the
  threshold instead of setting values below it to 0 as documented, which
  removes points from box/violin distributions.
- **#186 `plot_signatures()`** tests `exists("annotation_text")`, which
  finds a global variable of that name and adds it as a layer.
- **#95** `aes_string()` (deprecated) in `R/kmeans.R` (`plot_cluster`), and
  `is.ggplot()` (deprecated) in `test-plotting.R` and `test-kmeans.R`.

## Approach

- **#175.** `subset_tumor` is one or more values. Keep the samples whose
  value in any sample annotation column (other than `Samples`) exactly
  equals one of them. Stop, naming them, if some values are in no
  column. Subset the exposures by sample name (`all_of()`). The top
  annotation shows the `annotation` columns if given, otherwise the
  column(s) where the values were found, built with the existing
  `.annot_for_samples()` so it is matched by sample name. Scaling and
  proportions are still computed on all samples before subsetting, as now.
- **#176.** `t(scale(t(exp)))`: each signature is z-scored across samples.
  A signature with zero variance gets 0 instead of NaN. The `scale` doc
  says this.
- **#177.** Header boxes are filled with `color_mapping[mutation_color]`,
  by name. Label text is white on dark fills (brightness below 0.25 by
  `col2rgb`), black otherwise: same look for the defaults (only SBS96 C>G
  is dark). An unnamed `color_mapping` is named in the order the values
  first appear along the motifs (before, an unnamed mapping stopped with
  an error, found in spec review); a `color_mapping` missing a value
  stops, naming the values. When
  `color_mapping` is NULL and the table's default doesn't cover the
  chosen `color_variable`, colours are generated with `.discrete_colors()`.
- **#178.** A column name uses that column; a vector is matched to the
  motifs by name (`match(df$motif, rownames(signatures))`).
- **#179.** Box/violin with `group_by = "none"` puts signatures on the
  x-axis, as `group_by = "annotation"` does, without facets.
- **#180.** Limits are `c(0, max(breaks, max(exposure)))` in both
  `same_scale` branches; breaks unchanged.
- **#181.** `y = UMAP_2`.
- **#182.** `color_by = "cluster"` needs `clust` (with a message naming
  `cluster_exposure()`), not `annotation`. Fix the grammar of the
  annotation message; drop the two dead `annot <-` lines in that block.
- **#183.** `plot_umap()` stops with "No UMAP for model '<name>'. Run
  create_umap() first." A new helper `.proportional_exposures(exposures)`
  in `R/plotting.R` divides by sample totals and stops, naming the
  samples, when a total is 0. It replaces the five `sweep(..., "/")`
  calls: `create_umap()`, `plot_umap()`, `cluster_exposure()`,
  `plot_heatmap()`, `plot_exposures()`, and `k_select()` (added after
  spec review: it failed with "NA/NaN/Inf in foreign function call").
  (`plot_cluster()` has a seventh;
  it only runs after `create_umap()`, which now stops first, so it is
  left alone.)
- **#184.** `exposure[exposure < threshold] <- 0` instead of `filter()`;
  error message says "a single non-negative number".
- **#186.** `if (!is.null(annotation))`.
- **#95.** `aes(x = .data$UMAP_1, y = .data$UMAP_2, colour = .data$cluster)`
  in `plot_cluster()`; tests use `expect_s3_class(p, "ggplot")` (works
  with ggplot2 4.0.3's S7 objects, checked).

Out of scope (logged, not fixed): the Shiny heatmap tab offers values
only from a hard-coded `Tumor_Subtypes` column; `plot_cluster()`'s
"Run create_umap(result)" message and other copied validation (#187,
PR K).

## Tests

One regression test per issue, each failing on its assertion first:

- #175 `subset_tumor = "Lung"` gives a heatmap with the 8 Lung samples and
  a `Tumor_Subtypes` top annotation in sample order; two values; an
  unknown value errors with its name.
- #176 with `scale = TRUE`, row means are 0 and row sds 1.
- #177/#186 mock `ggpubr::ggarrange` with `local_mocked_bindings()` to
  capture the two panels: with `rev()` of the default mapping, each
  header box's fill equals `color_mapping[mutation_color]`; with a global
  `annotation_text` defined, the bar panel has one `geom_text` layer.
- #178 `.pivot_signatures()` with a vector and with `"context"`.
- #179 box and violin with `group_by = "none"` use `GeomBoxplot` and
  `GeomViolin`, x = signature.
- #180 counts with max 4549: no point is coloured `grey50`.
- #181 the centroid layer maps y to `UMAP_2`.
- #182 `color_by = "cluster"` with `clust` and no `annotation` works;
  without `clust` errors with "clust".
- #183 the two messages.
- #184 with a threshold, all sample x signature rows are kept and values
  below it are 0.
- #95 `make test` shows no `aes_string`/`is.ggplot` deprecation warnings.

## Other

- Docs: `scale`, `subset_tumor`, `threshold` (unchanged wording, now true),
  `color_variable`, `color_mapping`, `plot_umap` `clust`. `make docs`.
- NEWS under 2.7.15; version 2.7.14 -> 2.7.15.
- Vignette `musicatk.Rmd` doesn't call `plot_heatmap()`; the article
  `tutorial_tcga_cnsl.Rmd` does with the default `scale = TRUE`, so its
  heatmaps change. Maintainer decision (2026-11-11): update the article's
  text to match the new heatmaps and rebuild it with pkgdown. The agent
  edits the `.Rmd` and previews it with `make article
  FILTER=tutorial_tcga_cnsl`; the maintainer rebuilds `docs/` with
  `make site` (people-only).
- Checks: `make test`, `make coverage`, `make lint` on changed files;
  `check-full` and `bioccheck` before the PR.
