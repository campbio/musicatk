context("Test Find Signature Functions")
library("musicatk")

test_that(desc = "Inputs are correct", {
  expect_error(discover_signatures(data.frame(0)), regexp = "must be a")
  expect_error(predict_exposure(methods::new("musica", variants =
                                               data.table::data.table(0)),
                                "SBS96", cosmic_v2_sigs), regexp = "malformed")
})
