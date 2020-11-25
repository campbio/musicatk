context("Test Differential Analysis Functions")
library("musicatk")

test_that(desc = "Diff Anal Input validity", {
  expect_error(compare_samples(data.frame(0)), regexp = "must be a")
  data("res_annot")
  expect_error(compare_samples(res_annot, "dummy", method="wilcox"), 
               regexp="does not exist")
  expect_error(compare_samples(res_annot, "Tumor_Subtypes", method="dummy"), 
               regexp="Method is not supported")
})

test_that(desc = "Diff Anal functionality",  {
  data("res_annot")
  expect_type(compare_samples(res_annot, "Tumor_Subtypes", method="wilcox"), 
                  "double")
  expect_type(compare_samples(res_annot, "Tumor_Subtypes", method="kruskal"), 
              "double")
  expect_warning(compare_samples(res_annot, "Tumor_Subtypes", method="glm.nb"), 
                 regexp="NaNs produced")
  suppressWarnings(expect_type(compare_samples(res_annot, "Tumor_Subtypes", method="glm.nb"), 
              "double"))
})
