# Write the bundled example files in inst/extdata (#217, #238). Run from
# the package root: Rscript data-raw/example_files.R
#
# Source: the "Aliquot Ensemble Somatic Variant Merging and Masking" MAFs
# (WXS, hg38) of the NCI Genomic Data Commons, GDC Data Release 46.0
# (August 10, 2026), downloaded by file ID from the GDC API. The 12 samples
# were chosen by data-raw/choose_example_samples.R: four lung tumors each
# with a strong smoking (SBS4) or APOBEC (SBS2 + SBS13) signature, and four
# melanomas with a strong UV (SBS7a + SBS7b) signature.
#
# Writes:
# - public_<project>_<aliquot>.vcf: the smoking lung sample with the
#   highest smoking share
# - public_TCGA_lung.maf: the other seven lung samples, with the 34
#   standard TCGA MAF columns and the read counts
# - public_SKCM_<aliquot>.vcf: one VCF per melanoma
# - sample_annotations.txt: the tumor type (Lung or Skin) of each sample
#
# VCF indels are written with the hg38 base before the event, as VCF
# requires. No random numbers are used; the seed is set only so that later
# edits stay reproducible. Data-derived objects are rebuilt from these files
# by data-raw/example_musica.R.

pkgload::load_all(".", quiet = TRUE)
set.seed(1)

# In the order of choose_example_samples.R: each group by its share
chosen <- data.frame(
  group = rep(c("smoking", "apobec", "uv"), each = 4),
  project = c(
    "TCGA-LUSC", "TCGA-LUAD",
    "TCGA-LUAD", "TCGA-LUAD",
    "TCGA-LUSC", "TCGA-LUAD",
    "TCGA-LUAD", "TCGA-LUAD",
    "TCGA-SKCM", "TCGA-SKCM",
    "TCGA-SKCM", "TCGA-SKCM"
  ),
  aliquot = c(
    "TCGA-22-5472-01A-01D-1632-08", "TCGA-95-7562-01A-11D-2238-08",
    "TCGA-44-2655-01A-01D-1553-08", "TCGA-44-6774-01A-21D-1855-08",
    "TCGA-92-8065-01A-11D-2244-08", "TCGA-55-A48Y-01A-11D-A24D-08",
    "TCGA-91-6830-01A-11D-1945-08", "TCGA-95-7944-01A-11D-2184-08",
    "TCGA-ER-A3PL-06A-11D-A23B-08", "TCGA-EB-A3XD-01A-22D-A23B-08",
    "TCGA-BF-A1PX-01A-12D-A19A-08", "TCGA-EB-A41B-01A-11D-A24R-08"
  ),
  file_id = c(
    "9acd8b46-e076-4a5f-8f40-752c11f7632b",
    "897fa266-afac-4ab7-90d0-5046444486fc",
    "3c34b137-098b-4e76-9ed1-52e1be916fd0",
    "80f7616d-b100-45a1-92c5-e804690ccbd6",
    "ef05797d-30c8-430f-904d-8cd40bc6f3d2",
    "ad79d2d6-8895-4d5a-9293-bd019d9ceed2",
    "09400a18-aff7-4b9d-ac61-c65239c28b11",
    "d26d1add-37bd-4a6b-83ab-97e23bc9b7a1",
    "3b576835-5dad-4868-84da-8ab00146489e",
    "e6d31646-d920-47ab-8001-c59808876ae8",
    "8d913ab1-ef3b-44c2-add6-b9cd011c7ff4",
    "4693c34f-36d3-4bb7-8671-f1f593cc5941"
  )
)
extdata <- file.path("inst", "extdata")
maf_columns <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
  "Start_Position", "End_Position", "Strand", "Variant_Classification",
  "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
  "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
  "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
  "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
  "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
  "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
  "Verification_Status", "Validation_Status", "Mutation_Status",
  "Sequencing_Phase", "Sequence_Source", "Validation_Method", "Score",
  "BAM_File", "Sequencer", "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID",
  "t_depth", "t_ref_count", "t_alt_count", "n_depth"
)
g <- select_genome("hg38")

# Download a MAF by GDC file ID and keep the standard columns; the MAF must
# hold exactly its aliquot
read_gdc_maf <- function(file_id, aliquot) {
  path <- tempfile(fileext = ".maf.gz")
  on.exit(unlink(path), add = TRUE)
  utils::download.file(
    paste0("https://api.gdc.cancer.gov/data/", file_id), path,
    mode = "wb", quiet = TRUE
  )
  maf <- data.table::fread(path,
    skip = "Hugo_Symbol", quote = "", colClasses = "character",
    showProgress = FALSE, na.strings = NULL
  )
  if (!all(maf$Tumor_Sample_Barcode == aliquot)) {
    stop("The MAF ", file_id, " doesn't hold only ", aliquot)
  }
  maf[, maf_columns, with = FALSE]
}

# Write one sample's MAF records as a VCF with a tumor and a normal column
write_vcf <- function(maf, path, file_id) {
  pos <- as.integer(maf$Start_Position)
  ref <- maf$Reference_Allele
  alt <- maf$Tumor_Seq_Allele2
  ins <- ref == "-"
  del <- alt == "-"
  # VCF indels start at the base before the event: the MAF start of an
  # insertion, the base before the MAF start of a deletion
  pos[del] <- pos[del] - 1L
  anchor <- as.character(BSgenome::getSeq(g, GenomicRanges::GRanges(
    maf$Chromosome, IRanges::IRanges(pos, width = 1)
  )))
  ref[ins] <- anchor[ins]
  alt[ins] <- paste0(anchor[ins], maf$Tumor_Seq_Allele2[ins])
  ref[del] <- paste0(anchor[del], maf$Reference_Allele[del])
  alt[del] <- anchor[del]

  contigs <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
  records <- data.frame(
    CHROM = maf$Chromosome, POS = pos, ID = ".", REF = ref, ALT = alt,
    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT:AD:DP",
    tumor = paste0(
      "0/1:", maf$t_ref_count, ",", maf$t_alt_count, ":", maf$t_depth
    ),
    normal = paste0("0/0:.,.:", maf$n_depth)
  )
  records <- records[order(match(records$CHROM, contigs), records$POS), ]
  colnames(records)[1] <- "#CHROM"
  colnames(records)[10:11] <- c(
    maf$Tumor_Sample_Barcode[1], maf$Matched_Norm_Sample_Barcode[1]
  )
  header <- c(
    "##fileformat=VCFv4.2",
    paste0(
      "##source=GDC aliquot ensemble masked MAF ", file_id,
      ", converted by data-raw/example_files.R"
    ),
    "##reference=GRCh38",
    paste0(
      "##contig=<ID=", contigs, ",length=",
      GenomeInfoDb::seqlengths(g)[contigs], ">"
    ),
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste0(
      "##FORMAT=<ID=AD,Number=R,Type=Integer,",
      "Description=\"Allelic depths for the ref and alt alleles\">"
    ),
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
  )
  writeLines(header, path)
  data.table::fwrite(records, path,
    sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE
  )
}

mafs <- Map(read_gdc_maf, chosen$file_id, chosen$aliquot)
names(mafs) <- chosen$aliquot
# Remove all earlier example variant files, including the ones these
# replaced (public_LUAD_TCGA-97-7938.vcf, public_TCGA.LUSC.maf and three
# public_SKCM_*.vcf files)
unlink(list.files(extdata, pattern = "^public_", full.names = TRUE))

# The single lung VCF: the first smoking sample
vcf_sample <- chosen$aliquot[1]
vcf_project <- sub("TCGA-", "", chosen$project[1])
write_vcf(mafs[[vcf_sample]],
  file.path(extdata, paste0("public_", vcf_project, "_", vcf_sample, ".vcf")),
  chosen$file_id[1]
)

# The lung MAF: the other lung samples
lung <- chosen$aliquot[chosen$group != "uv" & chosen$aliquot != vcf_sample]
maf_path <- file.path(extdata, "public_TCGA_lung.maf")
writeLines("#version gdc-1.0.0", maf_path)
data.table::fwrite(data.table::rbindlist(mafs[lung]), maf_path,
  sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE
)

# The melanoma VCFs
for (i in which(chosen$group == "uv")) {
  write_vcf(mafs[[chosen$aliquot[i]]],
    file.path(extdata, paste0("public_SKCM_", chosen$aliquot[i], ".vcf")),
    chosen$file_id[i]
  )
}

utils::write.table(
  data.frame(
    Sample_Names = chosen$aliquot,
    Tumor_Subtypes = ifelse(chosen$group == "uv", "Skin", "Lung")
  ),
  file.path(extdata, "sample_annotations.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Check: every record is read back, and every reference base matches hg38
files <- c(
  list.files(extdata, pattern = "^public_.*vcf$", full.names = TRUE),
  maf_path
)
variants <- extract_variants(files, verbose = FALSE)
if (nrow(variants) != sum(vapply(mafs, nrow, 0L))) {
  stop("Expected ", sum(vapply(mafs, nrow, 0L)), " variants, got ",
    nrow(variants)
  )
}
invisible(withCallingHandlers(
  create_musica_from_variants(variants, g, verbose = FALSE),
  warning = function(w) stop("Unexpected warning: ", conditionMessage(w))
))
# The earlier example files took 911,845 bytes
size <- sum(file.size(list.files(extdata, full.names = TRUE)))
message(nrow(variants), " variants; inst/extdata is ", size, " bytes")
if (size > 911845) stop("inst/extdata grew")
