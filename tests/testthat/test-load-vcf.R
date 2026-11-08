context("Test Loading VCF File")
library("musicatk")

test_that(desc = "Testing VCF Input", {
  vcf_file <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
    package = "musicatk"
  )
  vcf <- musicatk::extract_variants_from_vcf_file(vcf_file)
  expect_s3_class(vcf, "data.table")
  expect_equal(nrow(vcf), 121)
})

test_that("multiallele = 'expand' keeps only the alleles each sample has", {
  # Issue #106: S1 carries A, S2 carries T, S3 carries both
  path <- write_test_vcf(
    list(c("chr1", "100", ".", "C", "A,T", ".", "PASS", ".", "GT",
           "0/1", "0/2", "1/2")),
    samples = c("S1", "S2", "S3")
  )
  s1 <- extract_variants_from_vcf_file(path, id = "S1")
  s2 <- extract_variants_from_vcf_file(path, id = "S2")
  s3 <- extract_variants_from_vcf_file(path, id = "S3")
  expect_identical(as.character(s1$alt), "A")
  expect_identical(as.character(s2$alt), "T")
  expect_identical(as.character(s3$alt), c("A", "T"))

  # Haploid and partly missing genotypes
  path <- write_test_vcf(
    list(c("chr1", "100", ".", "C", "A,T", ".", "PASS", ".", "GT",
           "2", "./1", "./.")),
    samples = c("S1", "S2", "S3")
  )
  expect_identical(extract_variants_from_vcf_file(path, id = "S1")$alt, "T")
  expect_identical(extract_variants_from_vcf_file(path, id = "S2")$alt, "A")
  expect_warning(
    expect_null(extract_variants_from_vcf_file(path, id = "S3")),
    regexp = "S3"
  )
})

test_that("multiallele = 'exclude' drops only multiallelic sites", {
  # Issue #110
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "C", "A,T", ".", "PASS", ".", "GT", "1/2"),
      c("chr1", "200", ".", "G", "A", ".", "PASS", ".", "GT", "0/1"),
      c("chr1", "300", ".", "G", "GTT", ".", "PASS", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  collapsed <- extract_variants_from_vcf_file(path, multiallele = "exclude")
  expect_type(collapsed$alt, "character")
  expect_identical(collapsed$alt, c("A", "GTT"))

  vcf <- VariantAnnotation::readVcf(path)
  expanded <- extract_variants_from_vcf(
    VariantAnnotation::expand(vcf),
    multiallele = "exclude"
  )
  expect_type(expanded$alt, "character")
  expect_identical(expanded$alt, c("A", "GTT"))
})

test_that("VCFs without a GT field keep all variants", {
  # Issue #111
  rows <- list(
    c("chr1", "100", ".", "C", "A", ".", "PASS", ".", "DP", "10"),
    c("chr1", "200", ".", "G", "T", ".", "PASS", ".", "DP", "12")
  )
  path <- write_test_vcf(rows,
    samples = "TUMOR",
    format = paste0(
      "##FORMAT=<ID=DP,Number=1,Type=Integer,",
      "Description=\"Depth\">"
    )
  )
  expect_message(
    dt <- extract_variants_from_vcf_file(path),
    regexp = "no GT field"
  )
  expect_identical(dt$alt, c("A", "T"))
  expect_identical(as.character(dt$sample), c("TUMOR", "TUMOR"))

  # With several sample columns, the sample must be chosen
  two_samples <- write_test_vcf(
    lapply(rows, c, "8"),
    samples = c("NORMAL", "TUMOR"),
    format = paste0(
      "##FORMAT=<ID=DP,Number=1,Type=Integer,",
      "Description=\"Depth\">"
    )
  )
  expect_error(
    extract_variants_from_vcf_file(two_samples),
    regexp = "'id' or 'rename'"
  )
  dt <- suppressMessages(
    extract_variants_from_vcf_file(two_samples, id = "TUMOR")
  )
  expect_identical(as.character(dt$sample), c("TUMOR", "TUMOR"))

  # A VCF with no sample columns needs a name for the sample
  sites_only <- write_test_vcf(
    lapply(rows, "[", seq_len(8)),
    samples = character(0), format = character(0)
  )
  expect_error(
    extract_variants_from_vcf_file(sites_only),
    regexp = "rename"
  )
  expect_error(
    extract_variants_from_vcf_file(sites_only, id = "TUMOR", rename = "T1"),
    regexp = "no sample columns"
  )
  dt <- suppressMessages(
    extract_variants_from_vcf_file(sites_only, rename = "T1")
  )
  expect_identical(as.character(dt$sample), c("T1", "T1"))

  empty <- write_test_vcf(list(),
    samples = character(0), format = character(0)
  )
  expect_warning(
    suppressMessages(
      extract_variants_from_vcf_file(empty, rename = "T1", filter = FALSE)
    ),
    regexp = "sample: T1\\. "
  )
})

test_that("fix_vcf_errors = TRUE recovers the variants of a malformed VCF", {
  # Issue #112: a REF of "0" makes readVcf() fail
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "C", "A", ".", "PASS", ".", "GT", "0/1"),
      c("chr1", "150", ".", "0", "A", ".", "PASS", ".", "GT", "0/1"),
      c("chr1", "200", ".", "G", "T", ".", "PASS", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  expect_warning(
    dt <- extract_variants_from_vcf_file(path),
    regexp = "1 malformed"
  )
  expect_s3_class(dt, "data.table")
  expect_identical(dt$start, c(100L, 200L))
  expect_error(
    extract_variants_from_vcf_file(path, fix_vcf_errors = FALSE),
    regexp = "malformed"
  )
})

test_that("extract_variants() reads gzipped VCF files", {
  # Issue #113
  path <- write_test_vcf(
    list(c("chr1", "100", ".", "C", "A", ".", "PASS", ".", "GT", "0/1")),
    samples = "S1"
  )
  # writeVcf(index = TRUE) writes a bgzipped "<file>.bgz"
  bgz_path <- VariantAnnotation::writeVcf(
    VariantAnnotation::readVcf(path),
    tempfile(fileext = ".vcf"),
    index = TRUE
  )
  gz_path <- sub("\\.bgz$", ".gz", bgz_path)
  file.rename(bgz_path, gz_path)
  dt <- extract_variants(gz_path, verbose = FALSE)
  expect_identical(dt$alt, "A")

  bgz_path <- VariantAnnotation::writeVcf(
    VariantAnnotation::readVcf(path),
    tempfile(fileext = ".vcf"),
    index = TRUE
  )
  dt <- extract_variants(bgz_path, verbose = FALSE)
  expect_identical(dt$alt, "A")
})

test_that("extra_fields keeps names and handles per-allele INFO fields", {
  # Issue #121
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "C", "A,T", ".", "PASS",
        "DP=10;AF=0.1,0.2;TAG=x,y,z", "GT", "1/2"),
      c("chr1", "200", ".", "G", "A", ".", "PASS",
        "DP=12;AF=0.3;TAG=w", "GT", "0/1")
    ),
    samples = "S1",
    info = c(
      "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
      "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
      "##INFO=<ID=TAG,Number=.,Type=String,Description=\"Tags\">"
    )
  )
  one <- extract_variants_from_vcf_file(path, extra_fields = "DP")
  expect_true("DP" %in% colnames(one))
  expect_identical(one$DP, c(10L, 10L, 12L))

  all <- extract_variants_from_vcf_file(path,
    extra_fields = c("DP", "AF", "TAG")
  )
  expect_equal(all$AF, c(0.1, 0.2, 0.3))
  expect_identical(all$TAG, c("x,y,z", "x,y,z", "w"))
})

test_that("filter = TRUE keeps FILTER '.' and names the values seen", {
  # Issue #214
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "C", "A", ".", "PASS", ".", "GT", "0/1"),
      c("chr1", "200", ".", "G", "T", ".", ".", ".", "GT", "0/1"),
      c("chr1", "300", ".", "G", "C", ".", "LowQual", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  dt <- extract_variants_from_vcf_file(path)
  expect_identical(dt$start, c(100L, 200L))

  failed <- write_test_vcf(
    list(
      c("chr1", "300", ".", "G", "C", ".", "LowQual", ".", "GT", "0/1"),
      c("chr1", "400", ".", "G", "C", ".", "NORMALPANEL", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  expect_warning(
    expect_null(extract_variants_from_vcf_file(failed)),
    regexp = "LowQual, NORMALPANEL"
  )
})

test_that("ALT '*' and symbolic alleles are left out", {
  # Issue #215
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "C", "A,*", ".", "PASS", ".", "GT", "1/2"),
      c("chr1", "200", ".", "G", "<DEL>", ".", "PASS", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  expect_message(
    dt <- extract_variants_from_vcf_file(path),
    regexp = "2 alternate alleles"
  )
  expect_identical(dt$start, 100L)
  expect_identical(dt$alt, "A")
})

test_that("VCF alleles are trimmed to their minimal form", {
  # Issue #213: multiallelic sites split into alleles that share bases
  path <- write_test_vcf(
    list(
      c("chr1", "100", ".", "ATG", "A,AT", ".", "PASS", ".", "GT", "1/2"),
      c("chr1", "200", ".", "AT", "A,GT", ".", "PASS", ".", "GT", "1/2"),
      c("chr1", "300", ".", "ATG", "AG", ".", "PASS", ".", "GT", "0/1"),
      c("chr1", "400", ".", "A", "GTT", ".", "PASS", ".", "GT", "0/1")
    ),
    samples = "S1"
  )
  dt <- extract_variants_from_vcf_file(path)
  expect_identical(dt$start, c(100L, 101L, 200L, 200L, 300L, 400L))
  expect_identical(dt$end, c(102L, 102L, 201L, 200L, 301L, 400L))
  expect_identical(dt$ref, c("ATG", "TG", "AT", "A", "AT", "A"))
  expect_identical(dt$alt, c("A", "T", "A", "G", "A", "GTT"))

  # None is removed as compound; only the complex A>GTT is removed
  g <- select_genome("hg38")
  expect_no_message(
    musica <- create_musica_from_variants(dt, g,
      check_ref_bases = FALSE, verbose = FALSE
    ),
    message = "compound"
  )
  v <- variants(musica)
  expect_identical(v$Variant_Type, c("DEL", "DEL", "DEL", "SBS", "DEL"))
  expect_identical(v$start, c(101L, 102L, 201L, 200L, 301L))
})
