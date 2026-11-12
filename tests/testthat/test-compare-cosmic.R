library(musicatk)

# compare_cosmic_v3() also draws a plot; keep it out of the test directory
.quiet_compare <- function(...) {
  path <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  compare_cosmic_v3(...)
}

# The same, without the message that names the catalogue used
.silent_compare <- function(...) suppressMessages(.quiet_compare(...))

# The comparison compare_cosmic_v3() should return against a given catalogue
.expected <- function(musica, modality, model_id, cosmic, threshold) {
  sig_compare(
    sig1 = signatures(musica, "result", modality, model_id),
    sig2 = signatures(cosmic), threshold = threshold
  )
}

# A small DBS78 model to compare, predicted from the v3.0 DBS signatures
.dbs_model <- function() {
  data(dbs_musica)
  withr::with_seed(1, predict_exposure(dbs_musica, "DBS78",
    cosmic_v3_dbs_sigs, "lda",
    model_id = "dbs", make_copy = TRUE
  ))
}

test_that("the v3.0 objects are still COSMIC v3.0 (ADR 0003)", {
  expect_identical(ncol(signatures(cosmic_v3_sbs_sigs)), 67L)
  expect_identical(ncol(signatures(cosmic_v3_dbs_sigs)), 11L)
  expect_identical(ncol(signatures(cosmic_v3_indel_sigs)), 17L)
  expect_identical(ncol(signatures(cosmic_v3_sbs_sigs_exome)), 65L)
  expect_identical(model_id(cosmic_v3_sbs_sigs), "Cosmic_v3_SBS_Genome")
  expect_identical(model_id(cosmic_v3_dbs_sigs), "Cosmic_v3_DBS_Genome")
  expect_identical(model_id(cosmic_v3_indel_sigs), "Cosmic_v3_Indel_Genome")
  expect_identical(model_id(cosmic_v3_sbs_sigs_exome), "Cosmic_v3_SNV_Exome")
})

test_that("the v3.6 objects have the published signature counts", {
  expect_identical(ncol(signatures(cosmic_v3_6_sbs_sigs_grch37)), 101L)
  expect_identical(ncol(signatures(cosmic_v3_6_dbs_sigs_grch37)), 22L)
  expect_identical(ncol(signatures(cosmic_v3_6_indel_sigs)), 25L)
  expect_identical(
    colnames(signatures(cosmic_v3_6_sbs_sigs_grch38)),
    colnames(signatures(cosmic_v3_6_sbs_sigs_grch37))
  )
  expect_identical(
    colnames(signatures(cosmic_v3_6_dbs_sigs_grch38)),
    colnames(signatures(cosmic_v3_6_dbs_sigs_grch37))
  )
  # GRCh37 and GRCh38 catalogues are normalized differently
  expect_false(isTRUE(all.equal(
    signatures(cosmic_v3_6_sbs_sigs_grch38),
    signatures(cosmic_v3_6_sbs_sigs_grch37)
  )))
  # Motifs in the same order as the v3.0 objects
  expect_identical(
    rownames(signatures(cosmic_v3_6_sbs_sigs_grch37)),
    rownames(signatures(cosmic_v3_sbs_sigs))
  )
  expect_identical(
    rownames(signatures(cosmic_v3_6_indel_sigs)),
    rownames(signatures(cosmic_v3_indel_sigs))
  )
})

test_that("the v3.6 objects carry no exposures and name version and build", {
  objects <- list(
    cosmic_v3_6_sbs_sigs_grch37 = "Cosmic_v3.6_SBS_Genome_GRCh37",
    cosmic_v3_6_sbs_sigs_grch38 = "Cosmic_v3.6_SBS_Genome_GRCh38",
    cosmic_v3_6_dbs_sigs_grch37 = "Cosmic_v3.6_DBS_Genome_GRCh37",
    cosmic_v3_6_dbs_sigs_grch38 = "Cosmic_v3.6_DBS_Genome_GRCh38",
    cosmic_v3_6_indel_sigs = "Cosmic_v3.6_Indel_Genome"
  )
  for (name in names(objects)) {
    obj <- get(name)
    expect_identical(model_id(obj), objects[[name]])
    expect_identical(dim(exposures(obj)), c(0L, 0L))
    expect_type(exposures(obj), "double")
  }
})

test_that("compare_cosmic_v3 defaults to the latest SBS catalogue per build", {
  data(res)
  expect_message(
    hg19 <- .quiet_compare(res, "res", "genome", threshold = 0.5),
    "COSMIC v3.6 SBS96 genome signatures \\(GRCh37\\)"
  )
  expect_equal(hg19, .expected(res, "SBS96", "res",
    cosmic_v3_6_sbs_sigs_grch37, 0.5
  ))
  expect_message(
    hg38 <- .quiet_compare(res, "res", "genome", genome_build = "hg38",
      threshold = 0.5
    ),
    "COSMIC v3.6 SBS96 genome signatures \\(GRCh38\\)"
  )
  expect_equal(hg38, .expected(res, "SBS96", "res",
    cosmic_v3_6_sbs_sigs_grch38, 0.5
  ))
  expect_equal(
    .silent_compare(res, "res", "genome", genome_build = "GRCh38",
      cosmic_version = "3.6", threshold = 0.5
    ),
    hg38
  )
})

test_that("compare_cosmic_v3(cosmic_version = '3.0') gives the v3.0 result", {
  data(res)
  expect_message(
    got <- .quiet_compare(res, "res", "genome", cosmic_version = "3.0",
      threshold = 0.5
    ),
    "COSMIC v3.0 SBS96 genome signatures \\(GRCh37\\)"
  )
  expect_equal(got, .expected(res, "SBS96", "res", cosmic_v3_sbs_sigs, 0.5))
  expect_equal(
    .silent_compare(res, "res", "genome", cosmic_version = 3.0,
      threshold = 0.5
    ),
    got
  )
  expect_equal(
    .silent_compare(res, "res", "genome", cosmic_version = "v3.0",
      threshold = 0.5
    ),
    got
  )
})

test_that("compare_cosmic_v3 picks the DBS catalogue by version and build", {
  m <- .dbs_model()
  expect_equal(
    .silent_compare(m, "dbs", "genome", modality = "DBS78",
      genome_build = "hg38", threshold = 0.5
    ),
    .expected(m, "DBS78", "dbs", cosmic_v3_6_dbs_sigs_grch38, 0.5)
  )
  expect_equal(
    .silent_compare(m, "dbs", "genome", modality = "DBS78",
      threshold = 0.5
    ),
    .expected(m, "DBS78", "dbs", cosmic_v3_6_dbs_sigs_grch37, 0.5)
  )
  expect_equal(
    .silent_compare(m, "dbs", "genome", modality = "DBS78",
      cosmic_version = "3.0", threshold = 0.5
    ),
    .expected(m, "DBS78", "dbs", cosmic_v3_dbs_sigs, 0.5)
  )
})

test_that("compare_cosmic_v3 has one Indel catalogue per version", {
  data(indel_musica)
  g <- select_genome("hg38")
  m <- indel_musica
  build_standard_table(m, g, "IND83")
  m <- withr::with_seed(1, predict_exposure(m, "IND83",
    cosmic_v3_indel_sigs, "lda",
    model_id = "ind", make_copy = TRUE
  ))
  latest <- .expected(m, "IND83", "ind", cosmic_v3_6_indel_sigs, 0.5)
  for (build in c("hg19", "hg38")) {
    expect_message(
      got <- .quiet_compare(m, "ind", "genome", modality = "IND83",
        genome_build = build, threshold = 0.5
      ),
      "COSMIC v3.6 IND83 genome signatures \\(one catalogue for all builds\\)"
    )
    expect_equal(got, latest)
  }
  expect_equal(
    .silent_compare(m, "ind", "genome", modality = "IND83",
      genome_build = "hg38", cosmic_version = "3.0", threshold = 0.5
    ),
    .expected(m, "IND83", "ind", cosmic_v3_indel_sigs, 0.5)
  )
})

test_that("compare_cosmic_v3 uses the v3.0 exome catalogue and says why", {
  data(res)
  expect_message(
    got <- .quiet_compare(res, "res", "exome", threshold = 0.5),
    "COSMIC v3.0 SBS96 exome signatures \\(GRCh37\\).*after v3.0"
  )
  expect_equal(got, .expected(res, "SBS96", "res",
    cosmic_v3_sbs_sigs_exome, 0.5
  ))
  expect_no_warning(
    .silent_compare(res, "res", "exome", threshold = 0.5)
  )
})

test_that("compare_cosmic_v3 errors on a catalogue musicatk doesn't ship", {
  data(res)
  expect_error(
    .quiet_compare(res, "res", "exome", cosmic_version = "3.6"),
    "no COSMIC v3.6 SBS96 exome catalogue.*v3.0 GRCh37"
  )
  expect_error(
    .quiet_compare(res, "res", "exome", genome_build = "hg38"),
    "no COSMIC SBS96 exome catalogue for GRCh38.*v3.0 GRCh37"
  )
  expect_error(
    .quiet_compare(res, "res", "genome", genome_build = "hg38",
      cosmic_version = "3.0"
    ),
    "no COSMIC v3.0 SBS96 genome catalogue for GRCh38.*v3.6 GRCh38"
  )
})

test_that("compare_cosmic_v3 checks cosmic_version and genome_build", {
  data(res)
  expect_error(
    .quiet_compare(res, "res", "genome", cosmic_version = "2.0"),
    "cosmic_version must be \"latest\" or one of: 3.0, 3.6"
  )
  expect_error(
    .quiet_compare(res, "res", "genome", cosmic_version = NULL),
    "cosmic_version must be"
  )
  expect_error(
    .quiet_compare(res, "res", "genome", genome_build = "mm10"),
    "genome_build must be one of"
  )
  # e.g. a Shiny input that has not rendered yet
  expect_error(
    .quiet_compare(res, "res", "genome", genome_build = NULL),
    "genome_build must be one of"
  )
})
