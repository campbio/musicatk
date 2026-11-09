# Choose the TCGA samples of the bundled example files (#217). Run from the
# package root: Rscript data-raw/choose_example_samples.R [download_dir]
#
# Not needed to rebuild the example files: data-raw/example_files.R lists
# the chosen file IDs. This script shows how they were chosen and checks
# the choice against the current GDC release.
#
# It lists the open-access TCGA-LUAD, TCGA-LUSC and TCGA-SKCM "Aliquot
# Ensemble Somatic Variant Merging and Masking" MAFs (WXS, hg38), downloads
# them (about 1,500 files; pass download_dir to keep them) and picks four
# samples each with a strong smoking (SBS4), APOBEC (SBS2 + SBS13) or UV
# (SBS7a + SBS7b) signature:
#
# 1. Keep aliquots that aren't whole-genome-amplified DNA (analyte W or
#    X), are primary tumors (sample type 01; for SKCM also metastases, 06,
#    since most TCGA melanomas are metastases), and are the first by
#    barcode when a patient has several.
# 2. Keep those with 150 to 500 SBS (SNPs not adjacent to another SNP of
#    the sample, as create_musica_from_variants() merges adjacent ones):
#    enough for a clean profile, few enough to keep the files small.
# 3. Estimate exposures to eight COSMIC v3 exome signatures with
#    predict_exposure() and turn them into shares of each sample's total.
# 4. Rank by the group's share (ties by barcode) and take four: smoking
#    from LUAD/LUSC with an APOBEC share below 0.1, APOBEC from LUAD/LUSC
#    with a smoking share below 0.1, UV from SKCM. Stop if a pick has a
#    group share below 0.5.

pkgload::load_all(".", quiet = TRUE)
set.seed(1)
args <- commandArgs(trailingOnly = TRUE)
download_dir <- if (length(args) > 0) args[1] else tempfile("example_mafs")
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)

# List the files of each project
list_mafs <- function(project) {
  filters <- list(op = "and", content = list(
    list(op = "in", content = list(
      field = "cases.project.project_id", value = list(project)
    )),
    list(op = "in", content = list(
      field = "data_type", value = list("Masked Somatic Mutation")
    )),
    list(op = "in", content = list(
      field = "analysis.workflow_type",
      value = list("Aliquot Ensemble Somatic Variant Merging and Masking")
    )),
    list(op = "in", content = list(field = "access", value = list("open")))
  ))
  query <- paste0(
    "https://api.gdc.cancer.gov/files?size=5000&format=json&fields=file_id",
    "&filters=", utils::URLencode(
      jsonlite::toJSON(filters, auto_unbox = TRUE),
      reserved = TRUE
    )
  )
  ids <- jsonlite::fromJSON(query)$data$hits$file_id
  data.frame(file_id = ids, project = rep(project, length(ids)))
}
files <- do.call(rbind, lapply(
  c("TCGA-LUAD", "TCGA-LUSC", "TCGA-SKCM"), list_mafs
))
message(nrow(files), " open MAFs")

# Download the ones not already in download_dir
files$path <- file.path(download_dir, paste0(files$file_id, ".maf.gz"))
for (i in which(!file.exists(files$path))) {
  utils::download.file(
    paste0("https://api.gdc.cancer.gov/data/", files$file_id[i]),
    files$path[i],
    mode = "wb", quiet = TRUE
  )
}

# Aliquot and SBS count of each file
count_sbs <- function(path) {
  maf <- data.table::fread(path,
    skip = "Hugo_Symbol", quote = "", showProgress = FALSE,
    select = c(
      "Chromosome", "Start_Position", "Variant_Type", "Tumor_Sample_Barcode"
    )
  )
  snp <- maf[maf$Variant_Type == "SNP", ]
  snp <- snp[order(snp$Chromosome, snp$Start_Position), ]
  next_to <- diff(snp$Start_Position) == 1 &
    utils::head(snp$Chromosome, -1) == utils::tail(snp$Chromosome, -1)
  adjacent <- c(next_to, FALSE) | c(FALSE, next_to)
  data.frame(aliquot = maf$Tumor_Sample_Barcode[1], sbs = sum(!adjacent))
}
files <- cbind(files, do.call(rbind, lapply(files$path, count_sbs)))

# Aliquot filters: analyte, sample type, one per patient, SBS window. A MAF
# with no records has no aliquot.
message(sum(is.na(files$aliquot)), " MAFs without records")
files <- files[!is.na(files$aliquot), ]
fields <- strsplit(files$aliquot, "-")
analyte <- substr(vapply(fields, `[`, "", 5), 3, 3)
sample_type <- substr(vapply(fields, `[`, "", 4), 1, 2)
keep_type <- sample_type == "01" |
  (files$project == "TCGA-SKCM" & sample_type == "06")
files <- files[!analyte %in% c("W", "X") & keep_type, ]
files <- files[order(files$aliquot), ]
files <- files[!duplicated(substr(files$aliquot, 1, 12)), ]
files <- files[files$sbs >= 150 & files$sbs <= 500, ]
message(nrow(files), " aliquots with 150 to 500 SBS")

# Exposures to the eight COSMIC v3 exome signatures
g <- select_genome("hg38")
variants <- extract_variants(files$path, verbose = FALSE)
candidates <- create_musica_from_variants(variants, g, verbose = FALSE)
build_standard_table(candidates, g, "SBS96")
use <- c("SBS1", "SBS2", "SBS4", "SBS5", "SBS7a", "SBS7b", "SBS13", "SBS40")
predict_exposure(candidates, "SBS96",
  signature_res = cosmic_v3_sbs_sigs_exome, signatures_to_use = use,
  algorithm = "lda", model_id = "cosmic"
)
expos <- exposures(candidates, "result", "SBS96", "cosmic")
share <- t(expos) / colSums(expos)
files$smoking <- share[files$aliquot, "SBS4"]
files$apobec <- rowSums(share[files$aliquot, c("SBS2", "SBS13")])
files$uv <- rowSums(share[files$aliquot, c("SBS7a", "SBS7b")])

# Four per group
lung <- files$project != "TCGA-SKCM"
pick <- function(rows, share) {
  chosen <- files[rows, ]
  chosen <- chosen[order(-chosen[[share]], chosen$aliquot), ]
  chosen <- utils::head(chosen, 4)
  if (nrow(chosen) < 4 || any(chosen[[share]] < 0.5)) {
    stop("Fewer than four samples with a ", share, " share of 0.5 or more")
  }
  cbind(group = share, chosen)
}
chosen <- rbind(
  pick(lung & files$apobec < 0.1, "smoking"),
  pick(lung & files$smoking < 0.1, "apobec"),
  pick(files$project == "TCGA-SKCM", "uv")
)
chosen[, c("smoking", "apobec", "uv")] <-
  round(chosen[, c("smoking", "apobec", "uv")], 2)
print(chosen[, c(
  "group", "aliquot", "project", "file_id", "sbs", "smoking", "apobec", "uv"
)], row.names = FALSE)
