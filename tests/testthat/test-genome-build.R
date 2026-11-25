library(musicatk)

test_that("supported_genomes lists the builds for each table (#227)", {
  expect_identical(supported_genomes(), c("hg19", "hg38", "mm9", "mm10"))
  expect_identical(supported_genomes("Transcript_Strand"), c("hg19", "hg38"))
  expect_identical(supported_genomes("Replication_Strand"), c("hg19", "hg38"))
})

test_that("supported_genomes rejects other strand types", {
  expect_error(supported_genomes("Strand"), regexp = "strand_type must be")
  expect_error(
    supported_genomes(c("Transcript_Strand", "Replication_Strand")),
    regexp = "strand_type must be"
  )
})

test_that("select_genome stops clearly for more than one build", {
  expect_error(
    select_genome(c("hg19", "hg38")),
    regexp = "That genome build is not currently selectable"
  )
})

test_that("select_genome takes a build number", {
  expect_identical(select_genome(38), select_genome("hg38"))
})

test_that("annotate_transcript_strand takes build names like select_genome", {
  # Issue #223
  data(musica)
  by_name <- musica
  annotate_transcript_strand(by_name, "hg38", build_table = FALSE)
  for (build in c("HG38", "GRCh38", "38")) {
    other <- musica
    annotate_transcript_strand(other, build, build_table = FALSE)
    expect_identical(
      variants(other)$Transcript_Strand,
      variants(by_name)$Transcript_Strand
    )
  }
})

test_that("annotate_transcript_strand rejects mouse builds clearly", {
  data(musica)
  expect_error(
    annotate_transcript_strand(musica, "mm10"),
    regexp = "Please select hg19, hg38, or provide a TxDb object"
  )
})
