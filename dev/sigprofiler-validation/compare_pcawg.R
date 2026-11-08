# Compare musicatk IND83 with SigProfilerMatrixGenerator (SPMG) ID83 on a
# folder of PCAWG consensus indel VCFs (GRCh37). Part of the validation for
# issue #218; not part of the package. See README.md.
#
# Usage, from the package root:
#   Rscript dev/sigprofiler-validation/compare_pcawg.R <vcf_dir> <out_dir> \
#     <python>
#
# <vcf_dir> holds *.vcf.gz files without sample columns; each file is one
# sample, named after the file up to its first ".". <python> is a Python
# with SPMG and its GRCh37 reference.
#
# PCAWG data is controlled access. This script writes only aggregate tables
# (totals, agreement counts and bin-by-bin cross-tables, with no sample names
# or positions) to <out_dir>. Its working copies of the PASS records are
# deleted at the end.
#
# Both tools read the same working copies of the VCFs, which keep only the
# passing records: FILTER "PASS" or "." (PCAWG marks passing records with
# ".", which musicatk's default filter = TRUE drops). Records on MT are left
# out, as hg19's chrM is not the GRCh37 MT sequence.
# - musicatk reads the copies with extract_variants(filter = FALSE), adds
#   "chr" to the chromosome names and builds IND83 on hg19.
# - SPMG runs on the copies with GRCh37.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript compare_pcawg.R <vcf_dir> <out_dir> <python>")
}
vcf_dir <- args[1]
out_dir <- args[2]
python <- args[3]
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("dev", "sigprofiler-validation", "helpers.R"))
message("musicatk ", as.character(utils::packageVersion("musicatk")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
work <- file.path(out_dir, "work")

# Everything in a function, so the working copies are deleted even on an
# error
main <- function() {
  files <- list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
  ids <- sub("\\..*$", "", basename(files))
  g <- select_genome("hg19")

  # Working copies: passing records, without MT
  input_dir <- file.path(work, "input_vcf")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  records <- 0
  passing <- 0
  mt <- 0
  for (i in seq_along(files)) {
    lines <- readLines(gzfile(files[i]))
    header <- startsWith(lines, "#")
    fields <- strsplit(lines[!header], "\t", fixed = TRUE)
    chrom <- vapply(fields, `[`, character(1), 1)
    filter <- vapply(fields, `[`, character(1), 7)
    pass <- filter %in% c("PASS", ".")
    keep <- pass & chrom != "MT"
    records <- records + length(fields)
    passing <- passing + sum(pass)
    mt <- mt + sum(pass & chrom == "MT")
    writeLines(c(lines[header], lines[!header][keep]),
      file.path(input_dir, paste0(ids[i], ".vcf"))
    )
  }
  .run_spmg(python, input_dir, "pcawg", "GRCh37", file.path(work, "spmg"))

  # musicatk
  copies <- file.path(input_dir, paste0(ids, ".vcf"))
  ext <- .quietly(extract_variants(copies, rename = ids, filter = FALSE,
    verbose = FALSE
  ))
  v <- ext$value
  v$chr <- paste0("chr", v$chr)
  m <- .quietly(create_musica_from_variants(v, g))
  musica <- m$value
  tab <- .quietly(build_standard_table(musica, g, "IND83"))
  mus <- get_count_table(tables(musica)[["IND83"]])
  spmg <- .read_spmg_matrix(file.path(work, "spmg"), "pcawg", "ID83",
    rownames(mus)
  )

  # Counts, with sample names removed
  cmp <- .compare_counts(mus, spmg, "IND83")
  by_motif <- data.frame(
    motif = rownames(mus),
    musicatk = rowSums(mus),
    spmg = rowSums(spmg[rownames(mus), , drop = FALSE])
  )
  by_motif <- by_motif[by_motif$musicatk != by_motif$spmg, ]

  # Per-variant bins, reduced to counts
  mus_id <- .musicatk_indels(musica, g)
  spmg_id <- .read_spmg_indels(file.path(work, "spmg"))
  names(spmg_id)[names(spmg_id) == "spmg"] <- "spmg_bin"
  j <- .join_indels(mus_id, spmg_id, "musicatk", "spmg_bin")
  jj <- j$joined
  indel_counts <- do.call(rbind, lapply(c("INS", "DEL"), function(t) {
    k <- jj$type == t
    both <- k & !is.na(jj$musicatk) & !is.na(jj$spmg_bin)
    data.frame(
      type = t, matched = sum(both),
      agree = sum(both & jj$musicatk == jj$spmg_bin),
      only_musicatk = sum(k & !is.na(jj$musicatk) & is.na(jj$spmg_bin)),
      only_spmg = sum(k & is.na(jj$musicatk) & !is.na(jj$spmg_bin))
    )
  }))
  info <- data.frame(
    files = length(files), musicatk_samples = ncol(mus),
    spmg_samples = ncol(spmg), records = records, passing = passing,
    mt_left_out = mt, musicatk_indels = nrow(mus_id),
    spmg_indels = nrow(spmg_id)
  )
  messages <- unique(c(ext$messages, m$messages, tab$messages))

  write.csv(cmp$summary, file.path(out_dir, "pcawg_summary.csv"),
    row.names = FALSE
  )
  write.csv(by_motif, file.path(out_dir, "pcawg_motif_totals_different.csv"),
    row.names = FALSE
  )
  write.csv(indel_counts, file.path(out_dir, "pcawg_indel_counts.csv"),
    row.names = FALSE
  )
  write.csv(j$disagree, file.path(out_dir, "pcawg_indel_disagree.csv"),
    row.names = FALSE
  )
  write.csv(info, file.path(out_dir, "pcawg_info.csv"), row.names = FALSE)
  # Messages without sample names (the "no GT field" note names each file)
  messages <- unique(gsub("[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12}",
    "<sample>", messages
  ))
  writeLines(messages, file.path(out_dir, "pcawg_messages.txt"))
}
tryCatch(main(), finally = unlink(work, recursive = TRUE))
