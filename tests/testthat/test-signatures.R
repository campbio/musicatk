context("Test Find Signature Functions")
library("musicatk")

test_that(desc = "Inputs are correct", {
  #incomplete_musicatk <- readRDS(system.file("testdata", "musicatk",
  #                        package = "musicatk"))

  #expect_error(discover_signatures(incomplete_musicatk, "SBS96"),
  #             regexp = "malformed")
  #expect_error(predict_exposure(incomplete_musicatk, "SNV96", cosmic_v2_sigs),
  #             regexp = "malformed")

  expect_error(discover_signatures(data.frame(0)), regexp = "must be a")
  #expect_error(predict_exposure(methods::new("bagel", variants =
  #                                             data.table::data.table(0)),
  #                              "SNV96", cosmic_v2_sigs), regexp = "malformed")
})
