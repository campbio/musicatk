# Build the example musica objects from the files in inst/extdata (#217,
# #238). Run from the package root: Rscript data-raw/example_musica.R
#
# Source: the bundled example files, 12 TCGA samples on hg38 from GDC
# open-access MAFs, written by data-raw/example_files.R: four lung tumors
# each with a strong smoking or APOBEC signature (one lung VCF and a
# seven-sample lung MAF) and four melanomas with a strong UV signature (one
# VCF each), plus sample_annotations.txt with the tumor type of each sample.
#
# Builds, with the current extract_variants() and
# create_musica_from_variants() defaults:
# - musica: all variants of the 12 samples, with an SBS96 table
# - musica_annot: musica with the Tumor_Subtypes annotation
# - musica_sbs96: the SBS of musica, with an SBS96 table
# - musica_sbs96_tiny: the first 40 SBS by position of one smoking and one
#   APOBEC lung sample, with an SBS96 table
# - res, res_annot: musica and musica_annot with a 3-signature NMF model
#   on SBS96, with model_id "res" and "res_annot"
#
# The script stops unless the model shows the expected pattern: the closest
# COSMIC v3 exome signature of each signature is an APOBEC (SBS2, SBS13),
# smoking (SBS4) or UV (SBS7a, SBS7b) signature, one of each; the smoking
# and UV signatures have a cosine of at least 0.8 to it, and the APOBEC
# signature, which combines SBS2 and SBS13, a cosine of at least 0.9 to the
# best mix of the two; and each sample's largest exposure is its group's
# signature. Signatures are ordered APOBEC, smoking, UV.

pkgload::load_all(".", quiet = TRUE)
seed <- 1
set.seed(seed)

g <- select_genome("hg38")
extdata <- file.path("inst", "extdata")
files <- list.files(extdata, pattern = "^public_.*(vcf|maf)$",
  full.names = TRUE
)
variants <- extract_variants(files, verbose = FALSE)
if (nrow(variants) != 4233) {
  stop("Expected 4233 variants, got ", nrow(variants))
}

# The signature group of each sample, from data-raw/choose_example_samples.R
groups <- c(
  "TCGA-22-5472-01A-01D-1632-08" = "smoking",
  "TCGA-95-7562-01A-11D-2238-08" = "smoking",
  "TCGA-44-2655-01A-01D-1553-08" = "smoking",
  "TCGA-44-6774-01A-21D-1855-08" = "smoking",
  "TCGA-92-8065-01A-11D-2244-08" = "apobec",
  "TCGA-55-A48Y-01A-11D-A24D-08" = "apobec",
  "TCGA-91-6830-01A-11D-1945-08" = "apobec",
  "TCGA-95-7944-01A-11D-2184-08" = "apobec",
  "TCGA-ER-A3PL-06A-11D-A23B-08" = "uv",
  "TCGA-EB-A3XD-01A-22D-A23B-08" = "uv",
  "TCGA-BF-A1PX-01A-12D-A19A-08" = "uv",
  "TCGA-EB-A41B-01A-11D-A24R-08" = "uv"
)

# musica and musica_annot
musica <- create_musica_from_variants(variants, g, verbose = FALSE)
build_standard_table(musica, g, "SBS96")
musica_annot <- musica
annot <- read.table(file.path(extdata, "sample_annotations.txt"),
  sep = "\t", header = TRUE
)
annot_rows <- match(sample_names(musica_annot), annot$Sample_Names)
if (anyNA(annot_rows)) stop("Samples missing from sample_annotations.txt")
samp_annot(musica_annot, "Tumor_Subtypes") <- annot$Tumor_Subtypes[annot_rows]

# musica_sbs96: the SBS only
all_variants <- variants(musica)
sbs <- all_variants[all_variants$Variant_Type == "SBS", ]
musica_sbs96 <- create_musica_from_variants(sbs, g, verbose = FALSE)
build_standard_table(musica_sbs96, g, "SBS96")

# musica_sbs96_tiny: 40 SBS of a smoking and an APOBEC sample
tiny_samples <- c(
  "TCGA-95-7562-01A-11D-2238-08", "TCGA-92-8065-01A-11D-2244-08"
)
tiny <- sbs[sbs$sample %in% tiny_samples, ]
tiny <- tiny[order(
  tiny$sample, match(tiny$chr, GenomeInfoDb::seqnames(g)), tiny$start
), ]
rank_in_sample <- stats::ave(seq_len(nrow(tiny)), tiny$sample,
  FUN = seq_along
)
tiny <- tiny[rank_in_sample <= 40, ]
tiny$sample <- as.character(tiny$sample)
musica_sbs96_tiny <- create_musica_from_variants(tiny, g, verbose = FALSE)
build_standard_table(musica_sbs96_tiny, g, "SBS96")

# res and res_annot: NMF with three signatures, ordered APOBEC, smoking, UV
etiology <- c(
  SBS2 = "apobec", SBS13 = "apobec", SBS4 = "smoking", SBS7a = "uv",
  SBS7b = "uv"
)
group_order <- c("apobec", "smoking", "uv")

order_signatures <- function(musica, model_id) {
  sigs <- signatures(musica, "result", "SBS96", model_id)
  expos <- exposures(musica, "result", "SBS96", model_id)
  cosmic <- signatures(cosmic_v3_sbs_sigs_exome)
  cos <- .cosine(sigs, cosmic)
  closest <- colnames(cosmic)[apply(cos, 1, which.max)]
  best <- apply(cos, 1, max)
  message(
    "Closest COSMIC v3 exome signatures: ",
    toString(paste0(closest, " ", round(best, 2)))
  )
  sig_group <- unname(etiology[closest])
  if (!setequal(sig_group, group_order) || length(sig_group) != 3) {
    stop("Expected one APOBEC, one smoking and one UV signature")
  }
  new_order <- match(group_order, sig_group)
  # SBS2 and SBS13 occur together, so the APOBEC signature combines them:
  # compare it to the best mix of the two instead of to either alone
  apobec <- sigs[, new_order[1], drop = FALSE]
  mix_cos <- vapply(seq(0, 1, by = 0.05), function(w) {
    mix <- w * cosmic[, "SBS2"] + (1 - w) * cosmic[, "SBS13"]
    as.numeric(.cosine(apobec, matrix(mix, ncol = 1)))
  }, 0)
  message(
    "APOBEC signature vs the best SBS2/SBS13 mix: ", round(max(mix_cos), 2)
  )
  if (max(mix_cos) < 0.9 || any(best[new_order[2:3]] < 0.8)) {
    stop(
      "Expected an APOBEC signature with a cosine of at least 0.9 to a ",
      "mix of SBS2 and SBS13, and smoking and UV signatures with a cosine ",
      "of at least 0.8 to their closest COSMIC signature"
    )
  }
  sigs <- sigs[, new_order]
  expos <- expos[new_order, ]
  colnames(sigs) <- rownames(expos) <- paste0("Signature", seq_len(3))
  signatures(musica, "result", "SBS96", model_id) <- sigs
  exposures(musica, "result", "SBS96", model_id) <- expos
  musica
}

# Each sample's largest exposure must be its group's signature
check_pattern <- function(musica, model_id) {
  expos <- exposures(musica, "result", "SBS96", model_id)
  largest <- group_order[apply(expos, 2, which.max)]
  expected <- groups[colnames(expos)]
  share <- round(apply(expos, 2, max) / colSums(expos), 2)
  message(
    "Largest signature by sample: ",
    toString(paste0(expected, " ", largest, " ", share))
  )
  if (any(largest != expected)) {
    stop("A sample's largest exposure isn't its group's signature")
  }
}

res <- musica
discover_signatures(res, "SBS96", num_signatures = 3, algorithm = "nmf",
  model_id = "res", seed = seed, nstart = 10
)
res <- order_signatures(res, "res")
check_pattern(res, "res")

res_annot <- musica_annot
discover_signatures(res_annot, "SBS96", num_signatures = 3,
  algorithm = "nmf", model_id = "res_annot", seed = seed, nstart = 10
)
res_annot <- order_signatures(res_annot, "res_annot")
check_pattern(res_annot, "res_annot")

for (name in c(
  "musica", "musica_annot", "musica_sbs96", "musica_sbs96_tiny", "res",
  "res_annot"
)) {
  validObject(get(name))
  save(list = name, file = file.path("data", paste0(name, ".rda")),
    compress = "xz"
  )
}
