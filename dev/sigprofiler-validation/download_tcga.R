# Download the open-access GDC masked somatic mutation MAFs (hg38) for every
# TCGA project with TCGAbiolinks. Part of the SigProfilerMatrixGenerator
# validation for issue #218; not part of the package. See README.md.
#
# Usage: Rscript dev/sigprofiler-validation/download_tcga.R <out_dir>
#
# Writes <out_dir>/GDCdata/TCGA-<project>/.../<file>.maf.gz, one MAF per
# aliquot, as TCGAbiolinks::GDCdownload() lays them out, and
# <out_dir>/tcga_files.csv listing the project of each file.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript download_tcga.R <out_dir>")
}
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir)
gdc_dir <- file.path(out_dir, "GDCdata")
# GDCdownload() writes MANIFEST.txt to the working directory
old_wd <- setwd(out_dir)

projects <- TCGAbiolinks::getGDCprojects()$project_id
projects <- sort(projects[startsWith(projects, "TCGA-")])
message(length(projects), " TCGA projects")

files <- NULL
for (project in projects) {
  query <- TCGAbiolinks::GDCquery(
    project = project,
    data.category = "Simple Nucleotide Variation",
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
    experimental.strategy = "WXS",
    data.format = "maf",
    access = "open"
  )
  TCGAbiolinks::GDCdownload(query, directory = gdc_dir,
    files.per.chunk = 200
  )
  res <- TCGAbiolinks::getResults(query)
  paths <- list.files(file.path(gdc_dir, project), pattern = "\\.maf\\.gz$",
    recursive = TRUE, full.names = TRUE
  )
  message(project, ": ", nrow(res), " files in query, ", length(paths),
    " downloaded"
  )
  files <- rbind(files, data.frame(project = project, file = paths))
}
write.csv(files, file.path(out_dir, "tcga_files.csv"), row.names = FALSE)
setwd(old_wd)
