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
source(file.path("dev", "sigprofiler-validation", "helpers.R"))
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

.compare <- function(modality, ext) {
  mus <- get_count_table(tables(musica)[[modality]])
  spmg <- lapply(c("vcf", "maf"), function(project) {
    m <- .read_spmg_matrix(file.path(out_dir, paste0("spmg_", project)),
      project, ext, rownames(mus)
    )
    if (project == "vcf") colnames(m) <- .musicatk_sample(colnames(m))
    m
  })
  .compare_counts(mus, do.call(cbind, spmg), modality)
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
spmg_vcf <- .read_spmg_indels(file.path(out_dir, "spmg_vcf"))
spmg_vcf$sample <- .musicatk_sample(spmg_vcf$sample)
spmg_maf <- .read_spmg_indels(file.path(out_dir, "spmg_maf"))
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
insfix <- .read_spmg_indels(file.path(out_dir, "spmg_insfix"))
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
