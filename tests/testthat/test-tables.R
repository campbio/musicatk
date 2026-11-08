library(musicatk)

test_that("IND83 columns follow sample_names(), matched by name", {
  data(indel_musica)
  g <- select_genome("hg38")
  samples <- as.character(unique(variants(indel_musica)$sample))
  forward <- new("musica",
    variants = variants(indel_musica),
    sample_annotations = data.frame(Samples = samples)
  )
  reverse <- new("musica",
    variants = variants(indel_musica),
    sample_annotations = data.frame(Samples = rev(samples))
  )
  build_standard_table(forward, g, "IND83")
  build_standard_table(reverse, g, "IND83")
  forward_counts <- get_count_table(tables(forward)$IND83)
  reverse_counts <- get_count_table(tables(reverse)$IND83)

  expect_identical(colnames(forward_counts), samples)
  expect_identical(colnames(reverse_counts), rev(samples))
  expect_equal(reverse_counts[, samples], forward_counts)
})

test_that("Count tables have no columns for stale sample factor levels", {
  g <- select_genome("hg38")
  # Variants whose sample factor keeps levels for samples that are gone,
  # as left behind by subsetting in older versions
  add_stale_levels <- function(musica) {
    variants <- variants(musica)
    present <- sample_names(musica)
    variants$sample <- factor(as.character(variants$sample),
      levels = c("stale_1", present, "stale_2")
    )
    new("musica",
      variants = variants,
      sample_annotations = data.frame(Samples = present)
    )
  }

  data(musica_sbs96_tiny)
  tiny <- add_stale_levels(musica_sbs96_tiny)
  build_standard_table(tiny, g, "SBS96")
  expect_identical(
    colnames(get_count_table(tables(tiny)$SBS96)),
    sample_names(tiny)
  )

  data(dbs_musica)
  dbs <- add_stale_levels(dbs_musica)
  build_standard_table(dbs, g, "DBS78")
  expect_identical(
    colnames(get_count_table(tables(dbs)$DBS78)),
    sample_names(dbs)
  )
})

test_that("combine_count_tables() and custom tables align samples by name", {
  data(musica_sbs96)
  sbs96 <- tables(musica_sbs96)$SBS96
  counts <- get_count_table(sbs96)
  samples <- colnames(counts)

  musica <- new("musica",
    variants = variants(musica_sbs96),
    count_tables = list(SBS96 = sbs96),
    sample_annotations = data.frame(Samples = samples)
  )

  # Columns in another order are put back in sample order by name
  reversed_counts <- counts[, rev(samples)]
  expect_identical(
    .align_count_columns(musica, reversed_counts, "reversed"),
    counts
  )

  # A table in another order, as saved by older versions, is not bound by
  # position: combining stops and asks for the table to be rebuilt
  reversed <- sbs96
  reversed@count_table <- reversed_counts
  legacy <- musica
  legacy@count_tables$reversed <- reversed
  expect_error(
    combine_count_tables(legacy, c("SBS96", "reversed"), "combo"),
    regexp = "count table 'reversed'"
  )

  # Custom tables use sample_names() order
  custom <- new("musica",
    variants = variants(musica_sbs96),
    sample_annotations = data.frame(Samples = rev(samples))
  )
  build_custom_table(custom, "Variant_Type", "type_table")
  type_counts <- get_count_table(tables(custom)$type_table)
  expect_identical(colnames(type_counts), rev(samples))
  per_sample <- table(as.character(variants(musica_sbs96)$sample))
  expect_equal(colSums(type_counts), c(per_sample[rev(samples)]),
    ignore_attr = TRUE
  )

  # Tables whose samples differ from sample_names() are refused
  extra <- sbs96
  extra@count_table <- counts[, -1]
  wrong <- musica
  wrong@count_tables$extra <- extra
  expect_error(
    combine_count_tables(wrong, c("SBS96", "extra"), "combo"),
    regexp = "do not match the sample names"
  )
})

test_that(".subset_count_tables() selects samples by name, not factor code", {
  data(musica_sbs96)
  samples <- sample_names(musica_sbs96)
  # A factor whose integer codes point at other columns
  wanted <- factor(samples[c(6, 7)], levels = rev(samples))
  subset <- .subset_count_tables(musica_sbs96, wanted)
  expect_identical(
    colnames(get_count_table(subset$SBS96)),
    samples[c(6, 7)]
  )
})

test_that("SBS tables warn about variants whose motif isn't in the schema", {
  # Issue #135
  g <- select_genome("hg38")
  # hg38 chr1:43430030 is G, so the first ref doesn't match the genome;
  # chr1:5000 lies in the N gap at the start of chr1, so its flanks are N
  tab <- data.frame(
    chr = "chr1", start = c(43430030L, 43430032L, 5000L),
    end = c(43430030L, 43430032L, 5000L), ref = "C", alt = "A",
    sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    check_ref_bases = FALSE, verbose = FALSE
  )
  expect_warning(
    build_standard_table(musica, g, "SBS96"),
    regexp = "2 of 3 variants have motifs not in the SBS96 schema"
  )
  expect_equal(sum(get_count_table(tables(musica)$SBS96)), 1)

  dat <- variants(musica)
  dat$Transcript_Strand <- "T"
  variants(musica) <- dat
  expect_warning(
    build_standard_table(musica, g, "SBS192", "Transcript_Strand"),
    regexp = "1 of 3 variants have motifs not in the SBS192_Trans schema"
  )
})

# IND83 bin of each indel in 'tab' (one sample per row), named by sample
.ind83_bins <- function(tab, g) {
  musica <- suppressMessages(
    create_musica_from_variants(tab, g, verbose = FALSE)
  )
  suppressMessages(build_standard_table(musica, g, "IND83"))
  counts <- get_count_table(tables(musica)$IND83)
  bins <- apply(counts, 2, function(x) {
    paste(rownames(counts)[x > 0], collapse = ",")
  })
  bins[tab$sample]
}

# IND83 counts of a subset of indel_musica's variants, with all its samples
.ind83_subset <- function(keep) {
  data(indel_musica)
  var <- variants(indel_musica)
  samples <- as.character(unique(var$sample))
  musica <- new("musica",
    variants = var[keep, ],
    sample_annotations = data.frame(Samples = samples)
  )
  g <- select_genome("hg38")
  suppressMessages(build_standard_table(musica, g, "IND83"))
  get_count_table(tables(musica)$IND83)
}

test_that("IND83 tables build when some indel classes are missing", {
  # Issue #126
  data(indel_musica)
  var <- variants(indel_musica)
  ins <- var$Variant_Type == "INS"
  len <- ifelse(ins, nchar(var$alt), nchar(var$ref))
  subsets <- list(
    deletions_only = !ins,
    insertions_only = ins,
    no_1bp_deletions = ins | len > 1,
    one_long_deletion = seq_along(ins) == which(!ins & len > 1)[1]
  )
  for (keep in subsets) {
    expect_no_error(counts <- .ind83_subset(keep))
    expect_identical(dim(counts), c(83L, 5L))
    expect_equal(sum(counts), sum(keep))
  }
})

test_that("IND83 left microhomology uses the bases next to the deletion", {
  # Issue #127
  expect_equal(
    musicatk:::.micro_left("ACGTT", "TTTTTTTTTTTTTTTTTTTTTTGTT"), 3
  )
  # hg38 chr9:128695506-128695521 is TTAC|AATA|AAACTAAAT: 2 bp of
  # microhomology on the right (AA), none on the left
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr9", start = 128695510L, end = 128695513L, ref = "AATA",
    alt = "-", sample = "S1"
  )
  expect_equal(unname(.ind83_bins(tab, g)), "DEL_MH_4_2")
})

test_that("IND83 counts repeats of multi-base units on the left", {
  # Issue #128
  expect_equal(musicatk:::.count_repeat_left("AC", "TTACACAC"), 3)
  # hg38 chr1:2253479-2253490 is GG ACACACAC TT. Deleting any of the four
  # AC copies leaves three.
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1", start = c(2253481L, 2253483L, 2253487L),
    end = c(2253482L, 2253484L, 2253488L), ref = "AC", alt = "-",
    sample = c("first", "second", "last")
  )
  expect_equal(unname(.ind83_bins(tab, g)), rep("DEL_repeats_2_3", 3))
})

test_that("IND83 insertion flanks include the bases next to the insertion", {
  # Issue #129
  # hg38 chr1:153056348-153056363 is AAAGAAAC TTTTTGCT. Inserting a T
  # anywhere next to the run of five T gives INS_T_1_5+.
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1",
    start = c(153056355L, 153056355L, 153056356L),
    end = c(153056356L, 153056355L, 153056356L),
    ref = c("-", "C", "T"), alt = c("T", "CT", "TT"),
    sample = c("maf", "vcf", "vcf_shifted")
  )
  expect_equal(unname(.ind83_bins(tab, g)), rep("INS_T_1_5+", 3))
})

test_that("IND83 counts five or more copies of longer inserted units", {
  # Issue #130
  # hg38 chr1:2629380 starts seven copies of CAG after a C
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1", start = 2629379L, end = 2629379L, ref = "C",
    alt = "CCAG", sample = "S1"
  )
  expect_equal(unname(.ind83_bins(tab, g)), "INS_repeats_3_5+")
})

test_that("IND83 says when there are no indels or no variant types", {
  data(musica)
  g <- select_genome("hg38")
  var <- variants(musica)
  sbs_only <- new("musica",
    variants = var[var$Variant_Type == "SBS", ],
    sample_annotations = samp_annot(musica)
  )
  expect_error(
    build_standard_table(sbs_only, g, "IND83"),
    regexp = "No insertions or deletions found"
  )
  var$Variant_Type <- NULL
  untyped <- new("musica",
    variants = var,
    sample_annotations = samp_annot(musica)
  )
  expect_error(
    build_standard_table(untyped, g, "IND83"),
    regexp = "No Variant_Type annotation found"
  )
})

test_that("IND83 counts indels near chromosome ends", {
  # Issue #233: flanks reach past the start and end of chrM. hg38 chrM
  # (the rCRS) starts GATCACAGG.
  g <- select_genome("hg38")
  chrm_length <- GenomeInfoDb::seqlengths(g)[["chrM"]]
  last_base <- as.character(BSgenome::getSeq(g, "chrM", chrm_length,
    chrm_length
  ))
  tab <- data.frame(
    chr = "chrM", start = c(3L, 2L, chrm_length),
    end = c(4L, 3L, chrm_length), ref = c("-", "-", last_base),
    alt = c("A", "T", "-"), sample = c("first", "second", "last")
  )
  bins <- .ind83_bins(tab, g)
  # An A after GAT is a T on the other strand, between ATC and ACCTGTG
  expect_identical(unname(bins[c("first", "second")]),
    c("INS_T_1_0", "INS_T_1_1")
  )
  expect_match(bins[["last"]], "^DEL_[CT]_1_")

  # Chromosome lengths are looked up by name when 'chr' is a factor. Here
  # chr1 has code 25, the position of chrM in seqnames(g).
  insertion <- data.frame(
    chr = "chr1", start = 153056355L, end = 153056356L, ref = "-",
    alt = "T", sample = "S1"
  )
  musica <- create_musica_from_variants(insertion, g, verbose = FALSE)
  var <- variants(musica)
  var$chr <- factor(var$chr,
    levels = c(paste0("chr", 2:22), "chrX", "chrY", "chrM", "chr1")
  )
  variants(musica) <- var
  suppressMessages(build_standard_table(musica, g, "IND83"))
  counts <- get_count_table(tables(musica)$IND83)
  expect_identical(counts["INS_T_1_5+", "S1"], 1L)

  # A chromosome missing from 'g' gives a clear error
  var$chr <- "chrFoo"
  variants(musica) <- var
  expect_error(
    suppressWarnings(
      build_standard_table(musica, g, "IND83", overwrite = TRUE)
    ),
    regexp = "chrFoo"
  )
})

test_that("IND83 stops on indels that aren't standardized", {
  # Issue #234
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1", start = 153056355L, end = 153056355L, ref = "C",
    alt = "CT", sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g,
    standardize_indels = FALSE, verbose = FALSE
  )
  expect_error(
    build_standard_table(musica, g, "IND83"),
    regexp = "standardize_indels = TRUE"
  )
})

test_that("IND83 warns about indels whose motif isn't in the schema", {
  # Issue #236
  g <- select_genome("hg38")
  tab <- data.frame(
    chr = "chr1", start = c(153056355L, 153056355L),
    end = c(153056356L, 153056356L), ref = "-", alt = c("T", "N"),
    sample = "S1"
  )
  musica <- create_musica_from_variants(tab, g, verbose = FALSE)
  expect_warning(
    build_standard_table(musica, g, "IND83"),
    regexp = "1 of 2 variants have motifs not in the IND83 schema"
  )
  expect_equal(sum(get_count_table(tables(musica)$IND83)), 1)
})
