# Choose the TCGA-LUAD aliquots for dbs_musica (#217). Run from the package
# root: Rscript data-raw/choose_dbs_aliquots.R [download_dir]
#
# Not needed to rebuild dbs_musica: data-raw/gdc_example_musica.R lists the
# chosen file IDs. This script shows how they were chosen and checks the
# choice against the current GDC release.
#
# It lists the open-access TCGA-LUAD "Aliquot Ensemble Somatic Variant
# Merging and Masking" MAFs (WXS, hg38), downloads all of them (about
# 120 MB; pass download_dir to keep them), counts the DBS of each aliquot
# (GDC DNP records plus adjacent SNPs, as create_musica_from_variants()
# merges them) and prints the five with the most DBS. Aliquots of
# whole-genome-amplified DNA (analyte code W or X in the barcode) are left
# out, since amplification can create artifacts. Ties are broken by
# aliquot barcode. No random numbers are used; the seed is set only so
# that later edits stay reproducible.

set.seed(1)
args <- commandArgs(trailingOnly = TRUE)
download_dir <- if (length(args) > 0) args[1] else tempfile("luad_mafs")
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)

# List the files
filters <- list(op = "and", content = list(
  list(op = "in", content = list(
    field = "cases.project.project_id", value = list("TCGA-LUAD")
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
file_ids <- jsonlite::fromJSON(query)$data$hits$file_id
message(length(file_ids), " open TCGA-LUAD MAFs")

# Download the ones not already in download_dir
paths <- file.path(download_dir, paste0(file_ids, ".maf.gz"))
for (i in which(!file.exists(paths))) {
  utils::download.file(
    paste0("https://api.gdc.cancer.gov/data/", file_ids[i]), paths[i],
    mode = "wb", quiet = TRUE
  )
}

# Count the DBS of each aliquot
count_dbs <- function(path) {
  maf <- data.table::fread(path,
    skip = "Hugo_Symbol", quote = "", showProgress = FALSE,
    select = c(
      "Chromosome", "Start_Position", "Variant_Type", "Tumor_Sample_Barcode"
    )
  )
  snp <- maf[maf$Variant_Type == "SNP", ]
  snp <- snp[order(snp$Chromosome, snp$Start_Position), ]
  adjacent <- sum(
    diff(snp$Start_Position) == 1 &
      utils::head(snp$Chromosome, -1) == utils::tail(snp$Chromosome, -1)
  )
  data.frame(
    aliquot = maf$Tumor_Sample_Barcode[1],
    dbs = sum(maf$Variant_Type == "DNP") + adjacent
  )
}
counts <- do.call(rbind, lapply(paths, count_dbs))
counts$file_id <- file_ids

# Leave out whole-genome-amplified DNA: the analyte code is the third
# character of the fifth barcode field (W in TCGA-17-Z031-01A-01W-...)
analyte <- substr(vapply(strsplit(counts$aliquot, "-"), `[`, "", 5), 3, 3)
counts <- counts[!analyte %in% c("W", "X"), ]
counts <- counts[order(-counts$dbs, counts$aliquot), ]
print(utils::head(counts, 5), row.names = FALSE)
