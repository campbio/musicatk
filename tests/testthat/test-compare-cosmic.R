library(musicatk)

# compare_cosmic_v3() also draws a plot; keep it out of the test directory
.quiet_compare <- function(...) {
  path <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  compare_cosmic_v3(...)
}

# The comparison compare_cosmic_v3() should return against a given catalogue
.expected <- function(musica, modality, model_id, cosmic, threshold) {
  sig_compare(
    sig1 = signatures(musica, "result", modality, model_id),
    sig2 = signatures(cosmic), threshold = threshold
  )
}

test_that("the COSMIC v3.6 objects have the published signature counts", {
  expect_identical(ncol(signatures(cosmic_v3_sbs_sigs)), 101L)
  expect_identical(ncol(signatures(cosmic_v3_dbs_sigs)), 22L)
  expect_identical(ncol(signatures(cosmic_v3_indel_sigs)), 25L)
  expect_identical(ncol(signatures(cosmic_v3_sbs_sigs_exome)), 65L)
  expect_identical(
    colnames(signatures(cosmic_v3_sbs_sigs_grch38)),
    colnames(signatures(cosmic_v3_sbs_sigs))
  )
  expect_identical(
    colnames(signatures(cosmic_v3_dbs_sigs_grch38)),
    colnames(signatures(cosmic_v3_dbs_sigs))
  )
  # GRCh37 and GRCh38 catalogues are normalized differently
  expect_false(isTRUE(all.equal(
    signatures(cosmic_v3_sbs_sigs_grch38), signatures(cosmic_v3_sbs_sigs)
  )))
})

test_that("compare_cosmic_v3 picks the SBS catalogue for genome_build", {
  data(res)
  hg38 <- .quiet_compare(res, "res", "genome", genome_build = "hg38",
    threshold = 0.5
  )
  expect_equal(hg38, .expected(res, "SBS96", "res",
    cosmic_v3_sbs_sigs_grch38, 0.5
  ))
  expect_equal(
    .quiet_compare(res, "res", "genome", genome_build = "GRCh38",
      threshold = 0.5
    ),
    hg38
  )
  # The default stays hg19 (GRCh37)
  expect_equal(
    .quiet_compare(res, "res", "genome", threshold = 0.5),
    .expected(res, "SBS96", "res", cosmic_v3_sbs_sigs, 0.5)
  )
})

test_that("compare_cosmic_v3 picks the DBS catalogue for genome_build", {
  data(dbs_musica)
  m <- withr::with_seed(1, predict_exposure(dbs_musica, "DBS78",
    cosmic_v3_dbs_sigs, "lda",
    model_id = "dbs", make_copy = TRUE
  ))
  expect_equal(
    .quiet_compare(m, "dbs", "genome", modality = "DBS78",
      genome_build = "hg38", threshold = 0.5
    ),
    .expected(m, "DBS78", "dbs", cosmic_v3_dbs_sigs_grch38, 0.5)
  )
  expect_equal(
    .quiet_compare(m, "dbs", "genome", modality = "DBS78",
      genome_build = "hg19", threshold = 0.5
    ),
    .expected(m, "DBS78", "dbs", cosmic_v3_dbs_sigs, 0.5)
  )
})

test_that("compare_cosmic_v3 has one Indel catalogue for both builds", {
  data(indel_musica)
  g <- select_genome("hg38")
  m <- indel_musica
  build_standard_table(m, g, "IND83")
  m <- withr::with_seed(1, predict_exposure(m, "IND83",
    cosmic_v3_indel_sigs, "lda",
    model_id = "ind", make_copy = TRUE
  ))
  expect_message(
    got <- .quiet_compare(m, "ind", "genome", modality = "IND83",
      genome_build = "hg38", threshold = 0.5
    ),
    "genome_build is ignored"
  )
  expect_equal(got, .expected(m, "IND83", "ind", cosmic_v3_indel_sigs, 0.5))
})

test_that("compare_cosmic_v3 checks genome_build and warns for exome", {
  data(res)
  expect_error(
    .quiet_compare(res, "res", "genome", genome_build = "mm10"),
    "genome_build must be one of"
  )
  # e.g. a Shiny input that has not rendered yet
  expect_error(
    .quiet_compare(res, "res", "genome", genome_build = NULL),
    "genome_build must be one of"
  )
  expect_warning(
    got <- .quiet_compare(res, "res", "exome", threshold = 0.5),
    "COSMIC v3.0"
  )
  expect_equal(got, .expected(res, "SBS96", "res",
    cosmic_v3_sbs_sigs_exome, 0.5
  ))
})

test_that("COSMIC objects carry no exposures, like the objects they replace", {
  for (obj in list(
    cosmic_v3_sbs_sigs, cosmic_v3_sbs_sigs_grch38, cosmic_v3_dbs_sigs,
    cosmic_v3_dbs_sigs_grch38, cosmic_v3_indel_sigs, cosmic_v3_sbs_sigs_exome
  )) {
    expect_identical(dim(exposures(obj)), c(0L, 0L))
    expect_type(exposures(obj), "double")
  }
})

test_that("only the SBS and DBS model ids name a genome build", {
  expect_identical(
    model_id(cosmic_v3_indel_sigs), "Cosmic_v3.6_Indel_Genome"
  )
  expect_identical(
    model_id(cosmic_v3_sbs_sigs), "Cosmic_v3.6_SBS_Genome_GRCh37"
  )
  expect_identical(
    model_id(cosmic_v3_dbs_sigs_grch38), "Cosmic_v3.6_DBS_Genome_GRCh38"
  )
})

test_that("compare_cosmic_v3 ignores genome_build for exome", {
  data(res)
  expect_warning(
    got <- .quiet_compare(res, "res", "exome", genome_build = "mm10",
      threshold = 0.5
    ),
    "COSMIC v3.0"
  )
  expect_equal(got, .expected(res, "SBS96", "res",
    cosmic_v3_sbs_sigs_exome, 0.5
  ))
})
