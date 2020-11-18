context("Test Differential Analysis Functions")
library("musicatk")

test_that(desc = "Diff Anal Input validity", {
  expect_error(compare_samples(data.frame(0)), regexp = "must be a")
  result <- readRDS(system.file("testdata", "res_annot.rds", package = "musicatk"))
  expect_error(compare_samples(result, "dummy", method="wilcox"), 
               regexp="does not exist")
  expect_error(compare_samples(result, "Tumor_Subtypes", method="dummy"), 
               regexp="Method is not supported")
})

# test_that(desc = "Diff Anal functionality") {
#   result <- readRDS(system.file("testdata", "res_annot.rds", package = "musicatk"))
#   
# 
# }
