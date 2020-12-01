context("Test Differential Analysis Functions")
library("musicatk")

test_that(desc = "Diff Anal Input validity", {
  expect_error(differential_exposure(data.frame(0)), regexp = "must be a")
  data("res_annot")
  expect_error(differential_exposure(res_annot, "dummy", method="wilcox"), 
               regexp="does not exist")
  expect_error(differential_exposure(res_annot, "Tumor_Subtypes", method="dummy"), 
               regexp="should be one of")
})

test_that(desc = "Diff Anal functionality",  {
  data("res_annot")
  expect_type(differential_exposure(res_annot, "Tumor_Subtypes", method="wilcox"), 
                  "double")
  expect_type(differential_exposure(res_annot, "Tumor_Subtypes", method="kruskal"), 
              "double")
  expect_warning(differential_exposure(res_annot, "Tumor_Subtypes", method="glm.nb"), 
                 regexp="NaNs produced")
  suppressWarnings(expect_type(differential_exposure(res_annot, "Tumor_Subtypes", method="glm.nb"), 
              "double"))
})
