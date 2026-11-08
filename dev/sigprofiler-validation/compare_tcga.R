# Compare musicatk with SigProfilerMatrixGenerator (SPMG) on the GDC TCGA MAFs
# from download_tcga.R, one project at a time. Part of the validation for
# issue #218; not part of the package. See README.md.
#
# Usage, from the package root:
#   Rscript dev/sigprofiler-validation/compare_tcga.R <tcga_dir> <out_dir> \
#     <python> [project ...]
#
# <tcga_dir> is the out_dir of download_tcga.R; <python> is a Python with
# SPMG and its GRCh38 reference. Without projects, all projects are run;
# projects already in <out_dir>/projects/ are skipped.
#
# For each project:
# - musicatk gets the MAF records with extract_variants_from_matrix(), as in
#   the TCGA tutorial article, and builds SBS96, DBS78 and IND83 with default
#   arguments, and SBS96 again with convert_dbs = FALSE, which counts every
#   base as an SBS as SPMG does. (extract_variants() on the files stops at
#   the first MAF with only silent mutations, because maftools::read.maf()
#   errors on it.)
# - SPMG runs on the same MAFs ("maf"), and on the same records written as
#   VCF ("vcf"): substitutions as they are, deletions anchored on the base
#   before Start_Position and insertions on the base at Start_Position (the
#   GDC convention), with the anchor taken from hg38.
# Writes per-project tables under <out_dir>/projects/<project>/ and, at the
# end, the combined tables in <out_dir>.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compare_tcga.R <tcga_dir> <out_dir> <python> ",
       "[project ...]")
}
tcga_dir <- args[1]
out_dir <- args[2]
python <- args[3]
suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
source(file.path("dev", "sigprofiler-validation", "helpers.R"))
message("musicatk ", as.character(utils::packageVersion("musicatk")))

gdc_dir <- file.path(tcga_dir, "GDCdata")
projects <- if (length(args) > 3) args[-(1:3)] else list.files(gdc_dir)
.project_files <- function(project) {
  list.files(file.path(gdc_dir, project), pattern = "\\.maf\\.gz$",
    recursive = TRUE, full.names = TRUE
  )
}
g <- select_genome("hg38")
maf_cols <- c(
  "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
  "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"
)

# The MAF records as VCF records, one data.frame per sample
.maf_to_vcf <- function(maf) {
  ref <- maf$Reference_Allele
  alt <- maf$Tumor_Seq_Allele2
  start <- as.integer(maf$Start_Position)
  ins <- ref == "-"
  del <- alt == "-"
  pos <- ifelse(del, start - 1L, start)
  anchor <- rep("", nrow(maf))
  indel <- ins | del
  anchor[indel] <- BSgenome::getSeq(g, maf$Chromosome[indel], pos[indel],
    pos[indel],
    as.character = TRUE
  )
  data.frame(
    sample = maf$Tumor_Sample_Barcode, chrom = maf$Chromosome, pos = pos,
    ref = ifelse(ins, anchor, ifelse(del, paste0(anchor, ref), ref)),
    alt = ifelse(ins, paste0(anchor, alt), ifelse(del, anchor, alt))
  )
}

.write_vcfs <- function(vcf, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (s in unique(vcf$sample)) {
    v <- vcf[vcf$sample == s, ]
    lines <- paste(v$chrom, v$pos, ".", v$ref, v$alt, ".", "PASS", ".",
      sep = "\t"
    )
    writeLines(
      c("##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", lines),
      file.path(dir, paste0(s, ".vcf"))
    )
  }
}

.run_project <- function(project, pdir) {
  work <- file.path(out_dir, "work", project)
  maf_files <- .project_files(project)

  # Raw MAF records, read as text so no record is dropped
  maf <- do.call(rbind, lapply(maf_files, function(f) {
    data.table::fread(f, skip = "Hugo_Symbol", select = maf_cols,
      colClasses = "character", showProgress = FALSE
    )
  }))
  maf <- as.data.frame(maf)

  # SPMG inputs: the MAFs uncompressed, and the anchored VCFs
  maf_dir <- file.path(work, "input_maf")
  dir.create(maf_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in maf_files) {
    writeLines(readLines(gzfile(f)),
      file.path(maf_dir, sub("\\.gz$", "", basename(f)))
    )
  }
  .write_vcfs(.maf_to_vcf(maf), file.path(work, "input_vcf"))
  .run_spmg(python, maf_dir, "maf", "GRCh38", file.path(work, "spmg_maf"))
  .run_spmg(python, file.path(work, "input_vcf"), "vcf", "GRCh38",
    file.path(work, "spmg_vcf")
  )

  # musicatk with default arguments, and with convert_dbs = FALSE
  ext <- .quietly(extract_variants_from_matrix(maf))
  m1 <- .quietly(create_musica_from_variants(ext$value, g))
  m2 <- .quietly(create_musica_from_variants(ext$value, g,
    convert_dbs = FALSE
  ))
  musica <- m1$value
  musica_all_sbs <- m2$value
  tabs <- .quietly({
    build_standard_table(musica, g, "SBS96")
    build_standard_table(musica, g, "DBS78")
    build_standard_table(musica, g, "IND83")
    build_standard_table(musica_all_sbs, g, "SBS96")
  })
  ct <- function(m, name) get_count_table(tables(m)[[name]])
  motifs <- list(
    SBS96 = rownames(ct(musica, "SBS96")),
    DBS78 = rownames(ct(musica, "DBS78")),
    IND83 = rownames(ct(musica, "IND83"))
  )
  spmg <- function(p, ext, name) {
    .read_spmg_matrix(file.path(work, paste0("spmg_", p)), p, ext,
      motifs[[name]]
    )
  }

  cmp <- list(
    .compare_counts(
      ct(musica, "SBS96"), spmg("maf", "SBS96", "SBS96"),
      "musicatk_vs_spmg_maf:SBS96"
    ),
    .compare_counts(
      ct(musica_all_sbs, "SBS96"), spmg("maf", "SBS96", "SBS96"),
      "musicatk_no_dbs_vs_spmg_maf:SBS96"
    ),
    .compare_counts(
      ct(musica, "DBS78"), spmg("maf", "DBS78", "DBS78"),
      "musicatk_vs_spmg_maf:DBS78"
    ),
    .compare_counts(
      ct(musica, "IND83"), spmg("maf", "ID83", "IND83"),
      "musicatk_vs_spmg_maf:IND83"
    ),
    .compare_counts(
      ct(musica_all_sbs, "SBS96"), spmg("vcf", "SBS96", "SBS96"),
      "musicatk_no_dbs_vs_spmg_vcf:SBS96"
    ),
    .compare_counts(
      ct(musica, "DBS78"), spmg("vcf", "DBS78", "DBS78"),
      "musicatk_vs_spmg_vcf:DBS78"
    ),
    .compare_counts(
      ct(musica, "IND83"), spmg("vcf", "ID83", "IND83"),
      "musicatk_vs_spmg_vcf:IND83"
    ),
    .compare_counts(
      spmg("maf", "ID83", "IND83"), spmg("vcf", "ID83", "IND83"),
      "spmg_maf_vs_spmg_vcf:IND83", "spmg_maf", "spmg_vcf"
    )
  )
  # Same column names for every comparison
  summary <- do.call(rbind, lapply(cmp, function(x) {
    s <- x$summary
    names(s) <- c("table", "samples", "x_total", "y_total",
                  "samples_identical", "samples_only_x", "samples_only_y",
                  "cells_different", "abs_difference")
    s
  }))
  cells <- do.call(rbind, lapply(cmp, function(x) {
    setNames(x$cells, c("table", "sample", "motif", "x", "y"))
  }))

  # Per-variant indel bins. SPMG reads MAF insertions one base left of the
  # GDC position; shift them back so each indel joins to itself
  # Records that change length but are not plain indels (e.g. GC>T, C>TT);
  # musicatk removes them, SPMG bins them as indels at Start_Position
  ref <- maf$Reference_Allele
  alt <- maf$Tumor_Seq_Allele2
  complex <- ref != "-" & alt != "-" & nchar(ref) != nchar(alt)
  complex_key <- paste(maf$Tumor_Sample_Barcode[complex],
    sub("^chr", "", maf$Chromosome[complex]), maf$Start_Position[complex]
  )
  mus_id <- .musicatk_indels(musica, g)
  maf_id <- .read_spmg_indels(file.path(work, "spmg_maf"))
  maf_id$pos <- maf_id$pos + (maf_id$type == "INS")
  vcf_id <- .read_spmg_indels(file.path(work, "spmg_vcf"))
  names(maf_id)[names(maf_id) == "spmg"] <- "spmg_maf"
  names(vcf_id)[names(vcf_id) == "spmg"] <- "spmg_vcf"
  joins <- list(
    musicatk_vs_spmg_vcf = .join_indels(
      mus_id, vcf_id, "musicatk", "spmg_vcf"
    ),
    spmg_maf_vs_spmg_vcf = .join_indels(
      maf_id, vcf_id, "spmg_maf", "spmg_vcf"
    )
  )
  indel_counts <- do.call(rbind, lapply(names(joins), function(n) {
    j <- joins[[n]]$joined
    x <- names(j)[6]
    y <- names(j)[7]
    both <- !is.na(j[[x]]) & !is.na(j[[y]])
    do.call(rbind, lapply(c("INS", "DEL"), function(t) {
      k <- j$type == t
      data.frame(
        comparison = n, type = t,
        matched = sum(both & k), agree = sum(both & k & j[[x]] == j[[y]]),
        only_x = sum(k & !is.na(j[[x]]) & is.na(j[[y]])),
        only_y = sum(k & is.na(j[[x]]) & !is.na(j[[y]])),
        only_y_complex = sum(k & is.na(j[[x]]) & !is.na(j[[y]]) &
                               paste(j$sample, j$chrom, j$pos) %in%
                                 complex_key)
      )
    }))
  }))
  disagree <- do.call(rbind, lapply(names(joins), function(n) {
    d <- joins[[n]]$disagree
    if (nrow(d) == 0) return(NULL)
    data.frame(comparison = n, x = d[[1]], y = d[[2]], n = d$n)
  }))
  disagree_variants <- do.call(rbind, lapply(names(joins), function(n) {
    j <- joins[[n]]$joined
    x <- names(j)[6]
    y <- names(j)[7]
    bad <- is.na(j[[x]]) | is.na(j[[y]]) | j[[x]] != j[[y]]
    if (!any(bad)) return(NULL)
    data.frame(comparison = n, j[bad, 1:5], x = j[[x]][bad], y = j[[y]][bad])
  }))

  # Records where Tumor_Seq_Allele1 (musicatk's ref column for MAFs) is not
  # Reference_Allele
  tsa1 <- maf$Tumor_Seq_Allele1 != maf$Reference_Allele
  messages <- c(ext$messages, m1$messages, tabs$messages)
  info <- data.frame(
    maf_records = nrow(maf),
    maf_samples = length(unique(maf$Tumor_Sample_Barcode)),
    tsa1_not_ref = sum(tsa1),
    tsa1_not_ref_is_alt = sum(tsa1 & maf$Tumor_Seq_Allele1 ==
                                maf$Tumor_Seq_Allele2),
    complex_ins = sum(complex & nchar(alt) > nchar(ref)),
    complex_del = sum(complex & nchar(alt) < nchar(ref)),
    musicatk_variants = nrow(variants(musica))
  )

  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  out <- list(
    summary = summary, cells = cells, indel_counts = indel_counts,
    indel_disagree = disagree, indel_disagree_variants = disagree_variants,
    info = info, messages = data.frame(message = messages)
  )
  for (n in names(out)) {
    if (!is.null(out[[n]]) && nrow(out[[n]]) > 0) {
      write.csv(cbind(project = project, out[[n]]),
        file.path(pdir, paste0(n, ".csv")),
        row.names = FALSE
      )
    }
  }
  unlink(work, recursive = TRUE)
  file.create(file.path(pdir, "done"))
}

for (project in projects) {
  pdir <- file.path(out_dir, "projects", project)
  if (file.exists(file.path(pdir, "done"))) next
  t0 <- Sys.time()
  .run_project(project, pdir)
  minutes <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  message(project, " done in ", minutes, " min")
}

# Combine the per-project tables
for (n in c("summary", "cells", "indel_counts", "indel_disagree",
            "indel_disagree_variants", "info", "messages")) {
  paths <- file.path(out_dir, "projects", projects, paste0(n, ".csv"))
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) next
  combined <- do.call(rbind, lapply(paths, read.csv))
  write.csv(combined, file.path(out_dir, paste0("tcga_", n, ".csv")),
    row.names = FALSE
  )
}
