# Build the example musica objects from the files in inst/extdata (#217,
# #238). Run from the package root: Rscript data-raw/example_musica.R
#
# Source: the bundled example files, all from TCGA on hg38 (GDC):
# - public_LUAD_TCGA-97-7938.vcf and the three public_SKCM_*.vcf files,
#   converted from GDC MuTect2 MAFs with maf2vcf
# - public_TCGA.LUSC.maf, three LUSC samples from a GDC MuTect2 MAF
# - sample_annotations.txt, the tumor type of each sample
#
# Builds, with the current extract_variants() and
# create_musica_from_variants() defaults:
# - musica: all variants of the seven samples, with an SBS96 table
# - musica_annot: musica with the Tumor_Subtypes annotation
# - musica_sbs96: the SBS of musica, with an SBS96 table
# - musica_sbs96_tiny: the first 40 SBS by position of two LUSC samples,
#   with an SBS96 table
# - res, res_annot: musica and musica_annot with a 3-signature NMF model
#   on SBS96, with model_id "res" and "res_annot"
#
# The script stops unless the model shows the expected pattern: one UV
# signature (closest COSMIC v3 signature SBS7*) carrying most of the
# melanoma (Skin) exposure, and two other signatures carrying most of the
# lung exposure. Signatures are put in a fixed order: the lung signature
# closest to COSMIC SBS13 first, the other lung signature second, the UV
# signature last.

pkgload::load_all(".", quiet = TRUE)
seed <- 1

g <- select_genome("hg38")
extdata <- file.path("inst", "extdata")
files <- c(
  list.files(extdata, pattern = "vcf$", full.names = TRUE),
  file.path(extdata, "public_TCGA.LUSC.maf")
)
variants <- extract_variants(files, verbose = FALSE)
if (nrow(variants) != 911) stop("Expected 911 variants, got ", nrow(variants))

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

# musica_sbs96_tiny: 40 SBS of each of two samples
tiny_samples <- c(
  "TCGA-56-7582-01A-11D-2042-08", "TCGA-94-7557-01A-11D-2122-08"
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

# res and res_annot: NMF with three signatures, in a fixed order
order_signatures <- function(musica, model_id) {
  sigs <- signatures(musica, "result", "SBS96", model_id)
  expos <- exposures(musica, "result", "SBS96", model_id)
  cosmic <- signatures(cosmic_v3_sbs_sigs)
  cos <- .cosine(sigs, cosmic)
  closest <- colnames(cosmic)[apply(cos, 1, which.max)]
  uv <- which(startsWith(closest, "SBS7"))
  if (length(uv) != 1) {
    stop("Expected one UV signature; closest COSMIC: ", toString(closest))
  }
  lung <- setdiff(seq_len(3), uv)
  lung <- lung[order(-cos[lung, "SBS13"])]
  new_order <- c(lung, uv)
  sigs <- sigs[, new_order]
  expos <- expos[new_order, ]
  colnames(sigs) <- rownames(expos) <- paste0("Signature", seq_len(3))
  signatures(musica, "result", "SBS96", model_id) <- sigs
  exposures(musica, "result", "SBS96", model_id) <- expos
  musica
}

check_pattern <- function(musica, model_id, annot) {
  expos <- exposures(musica, "result", "SBS96", model_id)
  skin <- annot$Sample_Names[annot$Tumor_Subtypes == "Skin"]
  lung <- annot$Sample_Names[annot$Tumor_Subtypes == "Lung"]
  uv_share_skin <- sum(expos[3, skin]) / sum(expos[, skin])
  uv_share_lung <- sum(expos[3, lung]) / sum(expos[, lung])
  message(
    "UV signature share: ", round(uv_share_skin, 2), " of Skin, ",
    round(uv_share_lung, 2), " of Lung exposure"
  )
  if (uv_share_skin < 0.8 || uv_share_lung > 0.2) {
    stop("The UV signature doesn't separate Skin from Lung samples")
  }
}

data(cosmic_v3_sbs_sigs)
res <- musica
discover_signatures(res, "SBS96", num_signatures = 3, algorithm = "nmf",
  model_id = "res", seed = seed, nstart = 10
)
res <- order_signatures(res, "res")
check_pattern(res, "res", annot)

res_annot <- musica_annot
discover_signatures(res_annot, "SBS96", num_signatures = 3,
  algorithm = "nmf", model_id = "res_annot", seed = seed, nstart = 10
)
res_annot <- order_signatures(res_annot, "res_annot")
check_pattern(res_annot, "res_annot", annot)

for (name in c(
  "musica", "musica_annot", "musica_sbs96", "musica_sbs96_tiny", "res",
  "res_annot"
)) {
  validObject(get(name))
  save(list = name, file = file.path("data", paste0(name, ".rda")),
    compress = "xz"
  )
}
