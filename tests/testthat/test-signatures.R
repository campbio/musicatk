context("Test Find Signature Functions")
library("BAGEL")

test_that(desc = "Inputs are correct", {
  #incomplete_bagel <- readRDS(system.file("testdata", "bagel.rds",
  #                        package = "BAGEL"))

  #expect_error(discover_signatures(incomplete_bagel, "SBS96"),
  #             regexp = "malformed")
  #expect_error(predict_exposure(incomplete_bagel, "SNV96", cosmic_v2_sigs),
  #             regexp = "malformed")

  expect_error(discover_signatures(data.frame(0)), regexp = "trying to get")
  #expect_error(predict_exposure(methods::new("bagel", variants =
  #                                             data.table::data.table(0)),
  #                              "SNV96", cosmic_v2_sigs), regexp = "malformed")
})
