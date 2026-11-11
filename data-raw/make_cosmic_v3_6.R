# Regenerate the COSMIC v3.6 signature objects shipped in data/.
#
# Source matrices were downloaded from https://cancer.sanger.ac.uk/signatures/
# (Human Cancer v3.6, May 2026) and live in data-raw/cosmic_v3.6/. COSMIC no
# longer publishes a separate exome-normalized SBS matrix, so
# cosmic_v3_sbs_sigs_exome is NOT regenerated here; it remains the COSMIC v3.0
# exome catalogue and is documented as such in R/cosmic_data.R.
#
# Run from the package root:  Rscript data-raw/make_cosmic_v3_6.R

library(musicatk)

raw_dir <- "data-raw/cosmic_v3.6"

# ---- label translation ------------------------------------------------------
# COSMIC and musicatk name the mutation motifs differently. musicatk row order
# is taken from the existing objects; we map by name and assert full coverage
# rather than relying on the two files happening to agree on order.

# A[C>A]A -> C>A_ACA
sbs_cosmic_to_musicatk <- function(x) {
  five <- substr(x, 1, 1)
  mut <- substr(x, 3, 5)
  three <- substr(x, 7, 7)
  paste0(mut, "_", five, substr(mut, 1, 1), three)
}

# AC>CA -> AC>NN_CA
dbs_cosmic_to_musicatk <- function(x) {
  parts <- strsplit(x, ">", fixed = TRUE)
  vapply(parts, function(p) paste0(p[1], ">NN_", p[2]), character(1))
}

# 1:Del:C:0 -> DEL_C_1_0 ; 5:Del:M:5 -> DEL_MH_5+_5+
id_cosmic_to_musicatk <- function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)
  vapply(parts, function(p) {
    size <- if (p[1] == "5") "5+" else p[1]
    type <- toupper(p[2])
    subtype <- switch(p[3], "R" = "repeats", "M" = "MH", p[3])
    count <- if (p[4] == "5") "5+" else p[4]
    paste(type, subtype, size, count, sep = "_")
  }, character(1))
}

# ---- build ------------------------------------------------------------------

read_cosmic <- function(file, translate, template) {
  tbl <- read.delim(file.path(raw_dir, file), check.names = FALSE)
  mat <- as.matrix(tbl[, -1, drop = FALSE])
  rownames(mat) <- translate(tbl$Type)

  # The shipped objects define the canonical motif order; reorder to match and
  # fail loudly if the two label sets disagree.
  want <- rownames(signatures(template))
  missing <- setdiff(want, rownames(mat))
  extra <- setdiff(rownames(mat), want)
  if (length(missing) || length(extra)) {
    stop(file, ": motif mismatch.\n  missing: ",
         paste(missing, collapse = ", "), "\n  unexpected: ",
         paste(extra, collapse = ", "))
  }
  mat <- mat[want, , drop = FALSE]

  # COSMIC columns are probabilities per signature; confirm before shipping.
  sums <- colSums(mat)
  if (any(abs(sums - 1) > 1e-4)) {
    stop(file, ": signature(s) do not sum to 1: ",
         paste(names(sums)[abs(sums - 1) > 1e-4], collapse = ", "))
  }
  mat
}

make_result <- function(mat, model_id, modality) {
  methods::new(
    "result_model",
    signatures = mat,
    exposures = matrix(),
    num_signatures = ncol(mat),
    modality = modality,
    model_id = model_id
  )
}

# Templates supply the canonical motif ordering for each modality.
data(cosmic_v3_sbs_sigs, cosmic_v3_dbs_sigs, cosmic_v3_indel_sigs,
     package = "musicatk")

specs <- list(
  list(obj = "cosmic_v3_sbs_sigs_grch38", file = "COSMIC_v3.6_SBS_GRCh38.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Genome_GRCh38"),
  list(obj = "cosmic_v3_dbs_sigs_grch38", file = "COSMIC_v3.6_DBS_GRCh38.txt",
       fn = dbs_cosmic_to_musicatk, tmpl = cosmic_v3_dbs_sigs,
       modality = "DBS78", model_id = "Cosmic_v3.6_DBS_Genome_GRCh38"),
  # GRCh37 objects keep the unsuffixed historical names: COSMIC itself defaults
  # to GRCh37, as do the objects these replace.
  list(obj = "cosmic_v3_sbs_sigs", file = "COSMIC_v3.6_SBS_GRCh37.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Genome_GRCh37"),
  list(obj = "cosmic_v3_dbs_sigs", file = "COSMIC_v3.6_DBS_GRCh37.txt",
       fn = dbs_cosmic_to_musicatk, tmpl = cosmic_v3_dbs_sigs,
       modality = "DBS78", model_id = "Cosmic_v3.6_DBS_Genome_GRCh37"),
  list(obj = "cosmic_v3_indel_sigs", file = "COSMIC_v3.6_ID_GRCh37.txt",
       fn = id_cosmic_to_musicatk, tmpl = cosmic_v3_indel_sigs,
       modality = "IND83", model_id = "Cosmic_v3.6_Indel_Genome_GRCh37")
)

for (s in specs) {
  mat <- read_cosmic(s$file, s$fn, s$tmpl)
  assign(s$obj, make_result(mat, s$model_id, s$modality))
  save(list = s$obj, file = file.path("data", paste0(s$obj, ".rda")),
       compress = "xz")
  message(sprintf("%-28s %d motifs x %d signatures", s$obj, nrow(mat),
                  ncol(mat)))
}
