# Compare musicatk SBS96, DBS78 and IND83 counts with SigProfilerMatrixGenerator
# (SPMG) on the bundled example files. One-off validation for issue #218; not
# part of the package. See README.md in this folder.
#
# Usage, from the package root (run_spmg.py first):
#   Rscript dev/sigprofiler-validation/compare.R inst/extdata <out_dir>
#
# <out_dir> is the folder run_spmg.py wrote to. The comparison tables are
# written there too. The package is loaded from the working tree with
# pkgload, so the comparison checks this checkout, not an installed version.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compare.R <extdata_dir> <out_dir>")
}
extdata_dir <- args[1]
out_dir <- args[2]
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
message("musicatk ", as.character(utils::packageVersion("musicatk")))

# musicatk with default arguments, as in the vignette
vcf_files <- list.files(extdata_dir, pattern = "\\.vcf$", full.names = TRUE)
maf_file <- file.path(extdata_dir, "public_TCGA.LUSC.maf")
variants <- extract_variants(c(vcf_files, maf_file), verbose = FALSE)
g <- select_genome("hg38")
musica <- create_musica_from_variants(variants, g)
build_standard_table(musica, g, "SBS96")
build_standard_table(musica, g, "DBS78")
build_standard_table(musica, g, "IND83")
var <- as.data.frame(variants(musica))

# SPMG names VCF samples after the file ("public_LUAD_TCGA-97-7938"); musicatk
# uses the tumour column of the VCF ("TCGA-97-7938-01A-11D-2167-08")
.musicatk_sample <- function(spmg_names) {
  vapply(spmg_names, function(s) {
    barcode <- sub("^public_[A-Z]+_", "", s)
    hit <- sample_names(musica)[startsWith(sample_names(musica), barcode)]
    if (length(hit) != 1) stop("No unique musicatk sample for ", s)
    hit
  }, character(1), USE.NAMES = FALSE)
}

# SPMG ID83 name to musicatk IND83 name: 1:Del:C:5 is DEL_C_1_5+,
# 2:Ins:R:0 is INS_repeats_2_0 and 5:Del:M:5 is DEL_MH_5+_5+
.id83_to_musicatk <- function(x) {
  if (length(x) == 0) return(character(0))
  p <- do.call(rbind, strsplit(x, ":"))
  size <- ifelse(p[, 1] == "5", "5+", p[, 1])
  n <- ifelse(p[, 4] == "5", "5+", p[, 4])
  type <- toupper(p[, 2])
  kind <- c(R = "repeats", M = "MH")[p[, 3]]
  ifelse(p[, 1] == "1", paste(type, p[, 3], "1", n, sep = "_"),
    paste(type, kind, size, n, sep = "_")
  )
}

.read_spmg <- function(modality, ext, motifs) {
  mats <- lapply(c("vcf", "maf"), function(project) {
    dir <- sub("[0-9]+$", "", ext)
    path <- file.path(
      out_dir, paste0("spmg_", project), dir,
      paste0(project, ".", ext, ".all")
    )
    m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
    colnames(m) <- if (project == "vcf") .musicatk_sample(colnames(m)) else
      colnames(m)
    m
  })
  m <- do.call(cbind, mats)
  if (modality == "IND83") {
    rownames(m) <- .id83_to_musicatk(rownames(m))
  } else {
    m <- .match_count_rows(m, motifs, modality)
  }
  m
}

.compare <- function(modality, ext) {
  mus <- get_count_table(tables(musica)[[modality]])
  spmg <- .read_spmg(modality, ext, rownames(mus))
  if (!setequal(rownames(spmg), rownames(mus))) {
    stop(modality, ": motif names don't match after conversion")
  }
  mus_raw <- mus
  spmg_raw <- spmg
  samples <- union(colnames(mus), colnames(spmg))
  full <- function(m) {
    out <- matrix(0, nrow(mus), length(samples),
      dimnames = list(rownames(mus), samples)
    )
    out[rownames(m), colnames(m)] <- m[rownames(m), ]
    out
  }
  mus <- full(mus)
  spmg <- full(spmg)
  diff <- which(mus != spmg, arr.ind = TRUE)
  cells <- data.frame(
    table = rep(modality, nrow(diff)),
    sample = samples[diff[, "col"]],
    motif = rownames(mus)[diff[, "row"]],
    musicatk = mus[diff], spmg = spmg[diff]
  )
  summary <- data.frame(
    table = modality, samples = length(samples),
    musicatk_total = sum(mus), spmg_total = sum(spmg),
    samples_identical = sum(colSums(mus != spmg) == 0),
    samples_only_musicatk = length(setdiff(samples, colnames(spmg_raw))),
    samples_only_spmg = length(setdiff(samples, colnames(mus_raw))),
    cells_different = nrow(diff),
    abs_difference = sum(abs(mus - spmg))
  )
  list(summary = summary, cells = cells)
}

res <- list(
  .compare("SBS96", "SBS96"),
  .compare("DBS78", "DBS78"),
  .compare("IND83", "ID83")
)
summary <- do.call(rbind, lapply(res, `[[`, "summary"))
cells <- do.call(rbind, lapply(res, `[[`, "cells"))
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(cells, file.path(out_dir, "cells_different.csv"), row.names = FALSE)

# Per-variant IND83 labels: musicatk's classifier against SPMG's seqInfo
indels <- var[var$Variant_Type %in% c("INS", "DEL"), ]
indels$musicatk <- .classify_indels(indels, g)
indels$chrom <- sub("^chr", "", indels$chr)
# SPMG's per-variant indel output (seqInfo): sample, chrom, position of
# the base before the indel, ID83 label, ref, alt, strand
.read_spmg_indels <- function(project) {
  files <- list.files(
    file.path(out_dir, paste0("spmg_", project), "vcf_files", "ID"),
    full.names = TRUE
  )
  rows <- lapply(files[file.size(files) > 0], read.delim,
    header = FALSE, colClasses = "character"
  )
  d <- do.call(rbind, rows)
  if (is.null(d)) {
    return(data.frame(
      sample = character(0), chrom = character(0), pos = integer(0),
      type = character(0), spmg = character(0), spmg_ref = character(0)
    ))
  }
  label <- sub("^[A-Z]:", "", d$V4)
  data.frame(
    sample = d$V1, chrom = d$V2, pos = as.integer(d$V3),
    type = ifelse(grepl(":Ins:", label), "INS", "DEL"),
    spmg = .id83_to_musicatk(label), spmg_ref = d$V5
  )
}
spmg_vcf <- .read_spmg_indels("vcf")
spmg_vcf$sample <- .musicatk_sample(spmg_vcf$sample)
spmg_maf <- .read_spmg_indels("maf")
# SPMG reads every MAF indel as starting one base before Start_Position.
# For a GDC MAF insertion Start_Position is already the base before it, so
# SPMG puts the insertion one base to the left; shift it back to join
spmg_maf$pos <- spmg_maf$pos + (spmg_maf$type == "INS")
spmg_id <- rbind(spmg_vcf, spmg_maf)
# SPMG reports deletions and insertions at the base before the event;
# musicatk deletions start at the first deleted base
indels$pos <- ifelse(indels$Variant_Type == "DEL",
  indels$start - 1L, indels$start
)
indels$type <- indels$Variant_Type
indel_cmp <- merge(indels[, c("sample", "chrom", "pos", "type", "ref", "alt",
                              "musicatk")],
  spmg_id, by = c("sample", "chrom", "pos", "type"), all = TRUE
)
indel_cmp$agree <- !is.na(indel_cmp$musicatk) & !is.na(indel_cmp$spmg) &
  indel_cmp$musicatk == indel_cmp$spmg
write.csv(indel_cmp, file.path(out_dir, "indels_per_variant.csv"),
  row.names = FALSE
)

# Variants musicatk did not count in a standard table
other <- var[!var$Variant_Type %in% c("SBS", "DBS", "INS", "DEL"), ]
write.csv(other, file.path(out_dir, "variants_not_counted.csv"),
  row.names = FALSE
)

# The MAF insertions given to SPMG as VCF records anchored on the base at
# Start_Position ("insfix" in run_spmg.py), against musicatk's MAF bins
insfix <- .read_spmg_indels("insfix")
insfix <- merge(insfix[, c("sample", "chrom", "pos", "spmg")],
  transform(indels[indels$type == "INS", c("sample", "chrom", "pos",
                                           "musicatk")],
    sample = substr(sample, 1, 12)
  ),
  by = c("sample", "chrom", "pos"), all.x = TRUE
)

# The synthetic run of three SBS and complex indel in musicatk
synthetic_vcf <- file.path(out_dir, "input_synthetic", "synthetic.vcf")
synthetic <- create_musica_from_variants(
  extract_variants(synthetic_vcf, rename = "synthetic", verbose = FALSE), g
)

print(summary, row.names = FALSE)
cat("\nDifferent cells:\n")
print(cells, row.names = FALSE)
cat("\nIndels where the tools disagree or only one has the variant:\n")
print(indel_cmp[!indel_cmp$agree, ], row.names = FALSE)
cat("\nmusicatk variant types:\n")
print(table(var$Variant_Type))
cat("\nMAF insertions anchored on Start_Position, SPMG against musicatk:\n")
print(insfix, row.names = FALSE)
cat("\nSynthetic VCF in musicatk:\n")
print(as.data.frame(variants(synthetic))[, c("chr", "start", "end", "ref",
                                             "alt", "Variant_Type")])
