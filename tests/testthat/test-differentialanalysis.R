context("Test Differential Analysis Functions")
library("musicatk")

test_that(desc = "Diff Anal Input validity", {
  expect_error(exposure_differential_analysis(data.frame(0)),
               regexp = "must be a")
  data("res_annot")
  samp_annot(res_annot, "test") <- c("C", "C", "C", "B", "B", "B", "A")
  expect_error(exposure_differential_analysis(res_annot, "DNE",
                                              method = "wilcox"),
               regexp = "does not exist")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "dummy"),
               regexp = "should be one of")
  expect_message(exposure_differential_analysis(res_annot, "Tumor_Subtypes",
                                              method = "wilcox",
                                              group1 = "test1",
                                              group2 = "test2"),
               regexp = "'annotations' is of length 2.")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox"),
               regexp = "are required for annotations")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = "test1",
                                              group2 = "a"),
               regexp = "does not exist in annotations")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = "a",
                                              group2 = "test1"),
               regexp = "does not exist in annotations")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = 1,
                                              group2 = "B"),
               regexp = "must be character vectors")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = c("A", "C"),
                                              group2 = "B"),
               regexp = "must be the same length")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = c("A", "C"),
                                              group2 = c("B", "C")),
               regexp = "must be unique")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = c("A", "t"),
                                              group2 = c("B", "C")),
               regexp = "does not exist in annotations")
  expect_error(exposure_differential_analysis(res_annot, "test",
                                              method = "wilcox",
                                              group1 = c("A", "B"),
                                              group2 = c("t", "C")),
               regexp = "does not exist in annotations")

})

test_that(desc = "Diff Anal functionality",  {
  data("res_annot")
  expect_type(exposure_differential_analysis(res_annot, "Tumor_Subtypes",
                                             method = "wilcox"),
                  "list")
  expect_type(exposure_differential_analysis(res_annot, "Tumor_Subtypes",
                                             method = "kruskal"),
              "list")
  expect_warning(exposure_differential_analysis(res_annot, "Tumor_Subtypes",
                                                method = "glm.nb"),
                 regexp = "NaNs produced")
})
