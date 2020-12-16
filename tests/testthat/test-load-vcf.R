context("Test Loading VCF File")
library("musicatk")

test_that(desc = "Testing VCF Input", {
  vcf_file <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
                          package = "musicatk")
  vcf <- musicatk::extract_variants_from_vcf_file(vcf_file)
  expect_s3_class(vcf, "data.table")
  expect_equal(nrow(vcf), 121)
})
