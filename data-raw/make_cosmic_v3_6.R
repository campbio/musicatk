# Regenerate the COSMIC v3.6 signature objects shipped in data/.
#
# Source matrices live in data-raw/cosmic_v3.6/:
# - The genome SBS, DBS and ID files were downloaded from
#   https://cancer.sanger.ac.uk/signatures/ (Human Cancer v3.6, May 2026).
# - The exome SBS files (COSMIC_v3.6_SBS_GRCh3[78]_exome.txt) come from
#   SigProfilerAssignment, which the COSMIC signature authors maintain:
#   https://github.com/SigProfilerSuite/SigProfilerAssignment, commit
#   ff61b0f56d43c916582b7c505ce72364c883dfbd, under
#   SigProfilerAssignment/data/Reference_Signatures/GRCh3[78]/. That commit's
#   v3.6 genome SBS files are identical to the COSMIC downloads above.
# Following ADR 0003 (dev/adr/0003-versioned-cosmic-signatures.md), each
# version gets its own objects named
# cosmic_v<major>_<minor>_<modality>_sigs[_<build>][_exome]. The unversioned
# cosmic_v3_* objects stay COSMIC v3.0 and are only read here, for their motif
# order.
#
# Run from the package root:  Rscript data-raw/make_cosmic_v3_6.R

# Use the package source, not an installed musicatk, for the class definitions
pkgload::load_all(".", export_all = FALSE)

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
    # Signature-only objects carry an empty exposure matrix
    exposures = matrix(numeric(0), 0, 0),
    num_signatures = ncol(mat),
    modality = modality,
    model_id = model_id
  )
}

# Templates supply the canonical motif ordering for each modality: the frozen
# v3.0 objects in data/ of this source tree, which this script never writes.
tmpl <- new.env()
for (f in c("cosmic_v3_sbs_sigs", "cosmic_v3_dbs_sigs",
            "cosmic_v3_indel_sigs")) {
  load(file.path("data", paste0(f, ".rda")), envir = tmpl)
}

specs <- list(
  list(obj = "cosmic_v3_6_sbs_sigs_grch37", file = "COSMIC_v3.6_SBS_GRCh37.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Genome_GRCh37"),
  list(obj = "cosmic_v3_6_sbs_sigs_grch38", file = "COSMIC_v3.6_SBS_GRCh38.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Genome_GRCh38"),
  list(obj = "cosmic_v3_6_sbs_sigs_grch37_exome",
       file = "COSMIC_v3.6_SBS_GRCh37_exome.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Exome_GRCh37"),
  list(obj = "cosmic_v3_6_sbs_sigs_grch38_exome",
       file = "COSMIC_v3.6_SBS_GRCh38_exome.txt",
       fn = sbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_sbs_sigs,
       modality = "SBS96", model_id = "Cosmic_v3.6_SBS_Exome_GRCh38"),
  list(obj = "cosmic_v3_6_dbs_sigs_grch37", file = "COSMIC_v3.6_DBS_GRCh37.txt",
       fn = dbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_dbs_sigs,
       modality = "DBS78", model_id = "Cosmic_v3.6_DBS_Genome_GRCh37"),
  list(obj = "cosmic_v3_6_dbs_sigs_grch38", file = "COSMIC_v3.6_DBS_GRCh38.txt",
       fn = dbs_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_dbs_sigs,
       modality = "DBS78", model_id = "Cosmic_v3.6_DBS_Genome_GRCh38"),
  # COSMIC publishes one Indel catalogue (labelled GRCh37) for all builds, so
  # neither its object name nor its model id names a build.
  list(obj = "cosmic_v3_6_indel_sigs", file = "COSMIC_v3.6_ID_GRCh37.txt",
       fn = id_cosmic_to_musicatk, tmpl = tmpl$cosmic_v3_indel_sigs,
       modality = "IND83", model_id = "Cosmic_v3.6_Indel_Genome")
)

for (s in specs) {
  mat <- read_cosmic(s$file, s$fn, s$tmpl)
  assign(s$obj, make_result(mat, s$model_id, s$modality))
  save(list = s$obj, file = file.path("data", paste0(s$obj, ".rda")),
       compress = "xz")
  message(sprintf("%-34s %d motifs x %d signatures", s$obj, nrow(mat),
                  ncol(mat)))
}
