# Build dbs_musica and indel_musica from GDC open-access TCGA MAFs (#217).
# Run from the package root: Rscript data-raw/gdc_example_musica.R
#
# Source: the "Aliquot Ensemble Somatic Variant Merging and Masking" MAFs
# (WXS, hg38) of the NCI Genomic Data Commons, GDC Data Release 46.0
# (August 10, 2026), downloaded by file ID from the GDC API. No random
# numbers are used; the seed is set only so that later edits stay
# reproducible.
#
# - dbs_musica: all variants of five TCGA-LUAD aliquots, with SBS96 and
#   DBS78 tables. The five are the open TCGA-LUAD aliquot MAFs with the
#   most DBS (GDC DNP records plus adjacent SNPs in the same sample) among
#   the 618 of Data Release 46.0, leaving out whole-genome-amplified DNA.
#   data-raw/choose_dbs_aliquots.R makes this choice. The object is about
#   50 KB, larger than the earlier dbs_musica (24 KB); five samples give
#   176 DBS, enough for DBS78 examples.
# - indel_musica: the insertions and deletions of five TCGA-UCEC aliquots,
#   the same aliquots as the earlier indel_musica, with no count tables.

pkgload::load_all(".", quiet = TRUE)
set.seed(1)

luad_files <- c(
  "TCGA-05-4382-01A-01D-1931-08" = "879297cc-5afe-48aa-9a38-e8b1c6fe0a32",
  "TCGA-86-8073-01A-11D-2238-08" = "0a250294-e760-418a-bed3-be259ae161e1",
  "TCGA-69-7979-01A-11D-2184-08" = "c3721df1-653c-408e-b3fa-b7c11c1e3413",
  "TCGA-L9-A7SV-01A-11D-A397-08" = "44a776b7-0eb2-463c-ab1d-80c654278f56",
  "TCGA-55-7994-01A-11D-2184-08" = "03add666-b87c-414c-8695-5667f3b599c1"
)
ucec_files <- c(
  "TCGA-5B-A90C-01A-11D-A37N-09" = "e43a4d87-20f5-486b-92db-a7421a09943f",
  "TCGA-A5-A0R7-01A-31D-A10B-09" = "2c39f9cc-e7c1-423f-bd70-25147c8d72f9",
  "TCGA-A5-A2K7-01A-11D-A17W-09" = "577f5df7-8d80-4f63-8f48-2369b19b6342",
  "TCGA-BK-A4ZD-01A-11D-A27P-09" = "aed050de-572e-4979-9e5b-0b4361b049d9",
  "TCGA-DF-A2KU-01A-11D-A17W-09" = "11b46b47-6409-4e21-96dc-6ba3bac52ea6"
)

# Download MAFs by GDC file ID into a temporary folder and read their
# variants; each MAF must hold exactly its own aliquot
read_gdc_mafs <- function(file_ids) {
  folder <- tempfile("gdc_mafs")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  paths <- file.path(folder, paste0(names(file_ids), ".maf.gz"))
  for (i in seq_along(file_ids)) {
    utils::download.file(
      paste0("https://api.gdc.cancer.gov/data/", file_ids[[i]]), paths[i],
      mode = "wb", quiet = TRUE
    )
  }
  variants <- extract_variants(paths, verbose = FALSE)
  if (!setequal(as.character(variants$sample), names(file_ids))) {
    stop("The MAFs don't hold the expected aliquots")
  }
  variants
}

g <- select_genome("hg38")

luad <- read_gdc_mafs(luad_files)
dbs_musica <- create_musica_from_variants(luad, g, verbose = FALSE)
build_standard_table(dbs_musica, g, "SBS96")
build_standard_table(dbs_musica, g, "DBS78")

ucec <- read_gdc_mafs(ucec_files)
ucec <- ucec[ucec$ref == "-" | ucec$alt == "-", ]
indel_musica <- create_musica_from_variants(ucec, g, verbose = FALSE)

for (name in c("dbs_musica", "indel_musica")) {
  validObject(get(name))
  save(list = name, file = file.path("data", paste0(name, ".rda")),
    compress = "xz"
  )
}
