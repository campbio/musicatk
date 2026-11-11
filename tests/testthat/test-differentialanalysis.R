context("Test Differential Analysis Functions")
library("musicatk")

test_that(desc = "Diff Anal Input validity", {
  # expect_error(exposure_differential_analysis(data.frame(0)), # error here
  #             regexp = "must be a")
  data("res_annot")
  samp_annot(res_annot, "test") <- rep(c("C", "B", "A"), c(6, 5, 1))
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "DNE",
      method = "wilcox"
    ),
    regexp = "does not exist"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "dummy"
    ),
    regexp = "should be one of"
  )
  expect_message(
    exposure_differential_analysis(res_annot, "res_annot", "Tumor_Subtypes",
      method = "wilcox",
      group1 = "test1",
      group2 = "test2"
    ),
    regexp = "'annotations' is of length 2."
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox"
    ),
    regexp = "are required for annotations"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = "test1",
      group2 = "a"
    ),
    regexp = "does not exist in annotations"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = "a",
      group2 = "test1"
    ),
    regexp = "does not exist in annotations"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = 1,
      group2 = "B"
    ),
    regexp = "must be character vectors"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = c("A", "C"),
      group2 = "B"
    ),
    regexp = "must be the same length"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = c("A", "C"),
      group2 = c("B", "C")
    ),
    regexp = "must be unique"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = c("A", "t"),
      group2 = c("B", "C")
    ),
    regexp = "does not exist in annotations"
  )
  expect_error(
    exposure_differential_analysis(res_annot, "res_annot", "test",
      method = "wilcox",
      group1 = c("A", "B"),
      group2 = c("t", "C")
    ),
    regexp = "does not exist in annotations"
  )
})

test_that(desc = "Diff Anal functionality", {
  data("res_annot")
  expect_type(
    exposure_differential_analysis(res_annot, "res_annot", "Tumor_Subtypes",
      method = "wilcox"
    ),
    "list"
  )
  expect_type(
    exposure_differential_analysis(res_annot, "res_annot", "Tumor_Subtypes",
      method = "kruskal"
    ),
    "list"
  )
  expect_warning(
    exposure_differential_analysis(res_annot, "res_annot", "Tumor_Subtypes",
      method = "glm.nb"
    ),
    regexp = "NaNs produced"
  )
})

test_that(desc = "Diff Anal pairs annotations with exposures by name", {
  data("res_annot")
  samp_annot(res_annot, "test") <- rep(c("A", "B", "C"), 4)
  # A model holding a subset of the samples, as discovery returns when it
  # drops zero-count samples, so exposure columns and annotation rows differ
  results <- result_list(res_annot)
  model <- results[["result"]]@modality[["SBS96"]][["res_annot"]]
  model@exposures <- model@exposures[, -2]
  results[["result"]]@modality[["SBS96"]][["res_annot"]] <- model
  result_list(res_annot) <- results

  exposure_matrix <- exposures(res_annot, "result", "SBS96", "res_annot")
  annot <- samp_annot(res_annot)
  groups <- annot$test[
    match(colnames(exposure_matrix), as.character(annot$Samples))
  ]
  expected <- matrixTests::row_kruskalwallis(exposure_matrix, groups)$pvalue

  observed <- exposure_differential_analysis(res_annot, "res_annot", "test",
    method = "kruskal"
  )
  expect_equal(observed[["p-value"]], expected)
})

# Exposures of res_annot with one annotation value per exposure column
.diff_fixture <- function() {
  data("res_annot", envir = environment())
  samp_annot(res_annot, "test") <- rep(c("A", "B", "C"), 4)
  e <- exposures(res_annot, "result", "SBS96", "res_annot")
  annot <- samp_annot(res_annot)
  index <- match(colnames(e), as.character(annot$Samples))
  list(
    musica = res_annot, exposures = e, test = annot$test[index],
    subtype = annot$Tumor_Subtypes[index]
  )
}

test_that(desc = "glm.nb labels the intercept and adjusts each group (#169)", {
  fx <- .diff_fixture()
  res <- suppressWarnings(exposure_differential_analysis(
    fx$musica, "res_annot", "test",
    method = "glm.nb"
  ))
  expect_equal(colnames(res), c(
    "Intercept(coef)", "B(coef)", "C(coef)",
    "Intercept(Std. Error)", "B(Std. Error)", "C(Std. Error)",
    "Intercept(z)", "B(z)", "C(z)",
    "Intercept(Pr(>|z|))", "B(Pr(>|z|))", "C(Pr(>|z|))",
    "B(fdr)", "C(fdr)", "ANOVA(Pr(>Chi))", "ANOVA(fdr)"
  ))
  groups <- factor(fx$test)
  for (g in c("B", "C")) {
    expected <- vapply(seq_len(nrow(fx$exposures)), function(i) {
      fit <- suppressWarnings(MASS::glm.nb(round(fx$exposures[i, ]) ~ groups))
      summary(fit)$coefficients[paste0("groups", g), "Pr(>|z|)"]
    }, numeric(1))
    expect_equal(res[[paste0(g, "(Pr(>|z|))")]], expected)
    expect_equal(res[[paste0(g, "(fdr)")]], p.adjust(expected, "BH"))
  }
})

test_that(desc = "glm.nb uses the first factor level as the reference", {
  fx <- .diff_fixture()
  samp_annot(fx$musica, "test") <- factor(
    rep(c("A", "B", "C"), 4),
    levels = c("C", "A", "B")
  )
  res <- suppressWarnings(exposure_differential_analysis(
    fx$musica, "res_annot", "test",
    method = "glm.nb"
  ))
  expect_equal(
    grep("\\(fdr\\)$", colnames(res), value = TRUE),
    c("A(fdr)", "B(fdr)", "ANOVA(fdr)")
  )
})

test_that(desc = "wilcox returns numeric p-values and FDRs (#171)", {
  fx <- .diff_fixture()
  res <- exposure_differential_analysis(fx$musica, "res_annot",
    "Tumor_Subtypes",
    method = "wilcox"
  )
  expect_equal(
    colnames(res),
    c("Signature", "Group1", "Group2", "p-value", "fdr")
  )
  expect_type(res[["p-value"]], "double")
  expect_type(res$fdr, "double")
  expected <- matrixTests::row_wilcoxon_twosample(
    fx$exposures[, fx$subtype == "Lung"],
    fx$exposures[, fx$subtype == "Skin"]
  )$pvalue
  expect_equal(res[["p-value"]], expected)
  expect_equal(res$fdr, p.adjust(expected, "BH"))
})

test_that(desc = "wilcox runs three comparisons, BH within each (#170)", {
  fx <- .diff_fixture()
  res <- exposure_differential_analysis(fx$musica, "res_annot", "test",
    method = "wilcox",
    group1 = c("A", "A", "B"), group2 = c("B", "C", "C")
  )
  n_sigs <- nrow(fx$exposures)
  expect_equal(nrow(res), 3 * n_sigs)
  expect_equal(res$Signature, rep(rownames(fx$exposures), 3))
  expect_equal(res$Group1, rep(c("A", "A", "B"), each = n_sigs))
  expect_equal(res$Group2, rep(c("B", "C", "C"), each = n_sigs))
  for (i in seq_len(3)) {
    rows <- seq_len(n_sigs) + (i - 1) * n_sigs
    x <- fx$exposures[, fx$test == res$Group1[rows[1]]]
    y <- fx$exposures[, fx$test == res$Group2[rows[1]]]
    expected <- matrixTests::row_wilcoxon_twosample(x, y)$pvalue
    expect_equal(res[["p-value"]][rows], expected)
    expect_equal(res$fdr[rows], p.adjust(expected, "BH"))
  }
})

test_that(desc = "plot_differential_analysis plots -log10(FDR) (#174)", {
  fx <- .diff_fixture()
  sigs <- rownames(fx$exposures)

  wil <- exposure_differential_analysis(fx$musica, "res_annot", "test",
    method = "wilcox", group1 = c("A", "A"), group2 = c("B", "C")
  )
  p <- plot_differential_analysis(wil)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$value, -log10(wil$fdr))
  expect_equal(p$data$signature, wil$Signature)
  expect_equal(
    as.character(p$data$comparison),
    paste(wil$Group1, "vs", wil$Group2)
  )

  kw <- exposure_differential_analysis(fx$musica, "res_annot", "test",
    method = "kruskal"
  )
  p <- plot_differential_analysis(kw)
  expect_s3_class(p, "ggplot")
  expect_equal(p$data$value, -log10(kw$fdr))
  expect_equal(p$data$signature, sigs)

  nb <- suppressWarnings(exposure_differential_analysis(
    fx$musica, "res_annot", "test",
    method = "glm.nb"
  ))
  p <- plot_differential_analysis(nb)
  expect_s3_class(p, "ggplot")
  expect_equal(
    p$data$value,
    -log10(c(nb[["B(fdr)"]], nb[["C(fdr)"]], nb[["ANOVA(fdr)"]]))
  )
  expect_equal(p$data$signature, rep(sigs, 3))
  expect_equal(
    as.character(p$data$comparison),
    rep(c("B", "C", "ANOVA"), each = length(sigs))
  )
})

test_that(desc = "plot_differential_analysis checks its input (#174)", {
  fx <- .diff_fixture()
  kw <- exposure_differential_analysis(fx$musica, "res_annot", "test",
    method = "kruskal"
  )
  expect_error(
    plot_differential_analysis(data.frame(a = 1)),
    "exposure_differential_analysis"
  )
  expect_warning(
    plot_differential_analysis(kw, "glm", 2),
    "ignored"
  )
})

test_that(desc = "glm.nb rejects groups named like its other columns", {
  fx <- .diff_fixture()
  samp_annot(fx$musica, "test") <- rep(c("A", "ANOVA", "C"), 4)
  expect_error(
    exposure_differential_analysis(fx$musica, "res_annot", "test",
      method = "glm.nb"
    ),
    "can't be named 'ANOVA'"
  )
})
