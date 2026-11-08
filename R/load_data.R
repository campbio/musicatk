#' @importFrom methods is
NULL

#' Helper function to load common human or mouse genomes
#'
#' @param x Select the hg19 or hg38 human genome or the mm9 or mm10
#' mouse genome in UCSC format
#' @return Returns BSgenome of given version
#' @examples
#' g <- select_genome(x = "hg38")
#' @export
select_genome <- function(x) {
  # Choose genome build version
  # Keep for now as a helper function
  if (tolower(x) %in% c("hg19", "19", "grch37")) {
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (tolower(x) %in% c("hg38", "38", "grch38")) {
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (tolower(x) %in% c("mm10", "grcm38")) {
    g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else if (tolower(x) %in% c("mm9", "mgscv37")) {
    g <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
  } else {
    stop(
      "That genome build is not currently selectable with this function. ",
      "Run BSgenome::available.genomes() to get a full list of availble ",
      "genomes or use BSgenome::forgeBSgenomeDataPkg to create a ",
      "custome genome package."
    )
  }
  return(g)
}

#' Extract variants from mutliple objects
#'
#' Chooses the correct function to extract variants from input based on
#' the class of the object or the file extension. Different types of objects
#' can be mixed within the list. For example, the list can include VCF files
#' and maf objects. Certain parameters such as \code{id} and \code{rename}
#' only apply to VCF objects or files and need to be individually specified
#' for each VCF. Therefore, these parameters should be suppied as a vector
#' that is the same length as the number of inputs. If other types of
#' objects are in the input list, then the value of \code{id} and \code{rename}
#' will be ignored for these items.
#'
#' @param inputs A vector or list of objects or file names. Objects can be
#' \linkS4class{CollapsedVCF}, \linkS4class{ExpandedVCF}, \linkS4class{MAF},
#' an object that inherits from \code{matrix} or \code{data.frame}, or
#' character strings that denote the path to a vcf or maf file.
#' @param id A character vector the same length as \code{inputs} denoting
#' the sample to extract from a vcf.
#' See \code{\link{extract_variants_from_vcf}} for more details.
#' Only used if the input is a vcf object or file. Default \code{NULL}.
#' @param rename  A character vector the same length as \code{inputs} denoting
#' what the same will be renamed to.
#' See \code{\link{extract_variants_from_vcf}} for more details.
#' Only used if the input is a vcf object or file. Default \code{NULL}.
#' @param sample_field Some algoriths will save the name of the
#' sample in the ##SAMPLE portion of header in the VCF.
#' See \code{\link{extract_variants_from_vcf}} for more details.
#' Default \code{NULL}.
#' @param filter Exclude variants that do not have a \code{PASS} in the
#' \code{FILTER} column of VCF inputs.
#' @param multiallele Multialleles are when multiple alternative variants
#' are listed in the same row in the vcf.
#' See \code{\link{extract_variants_from_vcf}} for more details.
#' Only used if the input is a vcf object or file. Default \code{"expand"}.
#' @param filename_as_id If set to \code{TRUE}, the file name will be used
#' as the sample name.
#' See \code{\link{extract_variants_from_vcf_file}} for more details.
#' Only used if the input is a vcf file. Default \code{TRUE}.
#' @param strip_extension Only used if \code{filename_as_id} is set to
#' \code{TRUE}. If set to \code{TRUE}, the file extention will be stripped
#' from the filename before setting the sample name.
#' See \code{\link{extract_variants_from_vcf_file}} for more details.
#' Only used if the input is a vcf file.
#' Default \code{c(".vcf",".vcf.gz",".gz")}
#' @param fix_vcf_errors Attempt to automatically fix VCF file
#' formatting errors.
#' See \code{\link{extract_variants_from_vcf_file}} for more details.
#' Only used if the input is a vcf file. Default \code{TRUE}.
#' @param chromosome_col The name of the column that contains the chromosome
#' reference for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"Chromosome"}.
#' @param start_col The name of the column that contains the start
#' position for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"Start_Position"}.
#' @param end_col The name of the column that contains the end
#' position for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"End_Position"}.
#' @param ref_col The name of the column that contains the reference
#' base(s) for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"Tumor_Seq_Allele1"}.
#' @param alt_col The name of the column that contains the alternative
#' base(s) for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"Tumor_Seq_Allele2"}.
#' @param sample_col The name of the column that contains the sample
#' id for each variant. Only used if the input is a matrix or data.frame.
#' Default \code{"sample"}.
#' @param extra_fields Optionally extract additional fields from all input
#' objects. Default \code{NULL}.
#' @param verbose Show progress of variant extraction. Default \code{TRUE}.
#' @return Returns a data.table of variants from a vcf
#' @examples
#' # Get loations of two vcf files and a maf file
#' luad_vcf_file <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "musicatk"
#' )
#' lusc_maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#'   package = "musicatk"
#' )
#' melanoma_vcfs <- list.files(system.file("extdata", package = "musicatk"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE
#' )
#'
#' # Read all files in at once
#' inputs <- c(luad_vcf_file, melanoma_vcfs, lusc_maf_file)
#' variants <- extract_variants(inputs = inputs)
#' table(variants$sample)
#'
#' # Run again but renaming samples in first four vcfs
#' new_name <- c(paste0("Sample", 1:4), NA)
#' variants <- extract_variants(inputs = inputs, rename = new_name)
#' table(variants$sample)
#'
#' @export
extract_variants <- function(inputs, id = NULL, rename = NULL,
                             sample_field = NULL,
                             filename_as_id = FALSE,
                             strip_extension = c(".vcf", ".vcf.gz", ".gz"),
                             filter = TRUE,
                             multiallele = c("expand", "exclude"),
                             fix_vcf_errors = TRUE,
                             extra_fields = NULL,
                             chromosome_col = "chr",
                             start_col = "start",
                             end_col = "end",
                             ref_col = "ref",
                             alt_col = "alt",
                             sample_col = "sample",
                             verbose = TRUE) {
  if (!is(inputs, "list")) {
    inputs <- as.list(inputs)
  }
  input_list <- vector("list", length(inputs))

  # Check arguments
  multiallele <- match.arg(multiallele)
  if (!is.null(rename)) {
    if (length(rename) != length(input_list)) {
      stop(
        "The length of 'rename' must be the same as the length of 'input'.",
        " Only vcf object or file.vcf inputs will be renamed."
      )
    }
  } else {
    rename <- rep(NULL, length(input_list))
  }
  if (!is.null(id)) {
    if (length(id) != length(input_list)) {
      stop(
        "The lenght of 'rename' must be the same as the length of 'input'.",
        " This only applies to vcf object or file.vcf inputs."
      )
    }
  } else {
    id <- rep(NULL, length(input_list))
  }

  pb <- utils::txtProgressBar(
    min = 0, max = length(input_list), initial = 0,
    style = 3
  )
  for (i in seq_along(inputs)) {
    input <- inputs[[i]]

    if (inherits(input, c("CollapsedVCF", "ExpandedVCF"))) {
      dt <- extract_variants_from_vcf(
        vcf = input,
        id = id[i],
        rename = rename[i],
        sample_field = sample_field,
        filter = filter,
        multiallele = multiallele,
        extra_fields = extra_fields
      )
    } else if (is(input, "MAF")) {
      dt <- extract_variants_from_maf(maf = input, extra_fields = extra_fields)
    } else if (inherits(input, c("matrix", "data.frame"))) {
      dt <- extract_variants_from_matrix(
        mat = input,
        extra_fields = extra_fields,
        chromosome_col = chromosome_col,
        start_col = start_col,
        end_col = end_col,
        ref_col = ref_col,
        alt_col = alt_col,
        sample_col = sample_col
      )
    } else if (is(input, "character")) {
      # Look past a compression suffix, so "x.vcf.gz" is read as a VCF
      file_type <- tools::file_ext(sub("\\.b?gz$", "", input))
      if (file_type == "vcf") {
        dt <- extract_variants_from_vcf_file(
          vcf_file = input,
          id = id[i],
          rename = rename[i],
          sample_field = sample_field,
          filename_as_id = filename_as_id,
          strip_extension = strip_extension,
          filter = filter,
          multiallele = multiallele,
          extra_fields = extra_fields,
          fix_vcf_errors = fix_vcf_errors
        )
      } else if (file_type == "maf") {
        dt <- extract_variants_from_maf_file(
          maf_file = input,
          extra_fields = extra_fields
        )
      } else {
        stop("Input file could not be automatically parsed: ", input, " ")
      }
    } else {
      stop(
        "Each input must be a collapedVCF/expandedVCF object, MAF object ",
        "file.vcf, file.maf, or a matrix/data.frame object.",
        "Item ", input, " is of class '", class(input), "'"
      )
    }
    input_list[[i]] <- dt

    utils::setTxtProgressBar(pb, i)
    if (isTRUE(verbose)) {
      message(
        "Extracted ", i, " out of ", length(inputs), " inputs: ",
        input
      )
    }
  }

  dt <- do.call("rbind", input_list)
  return(dt)
}

#' Extracts variants from a VariantAnnotation VCF object
#'
#' Variants are assigned to the sample selected by \code{id} using its GT
#' field. If the VCF has no GT field, all variants passing the filter are
#' assigned to the sample. A VCF without sample columns needs \code{rename}.
#'
#' Aaron - Need to describe differnce between ID, and name in the header, and
#' rename in terms of naming the sample. Need to describe differences in
#' multiallelic choices. Also need to describe the automatic error fixing
#'
#' @param vcf Location of vcf file
#' @param id ID of the sample to select from VCF. If \code{NULL}, then the
#' first sample will be selected. Default \code{NULL}.
#' @param rename Rename the sample to this value when extracting variants.
#' If \code{NULL}, then the sample will be named according to \code{ID}.
#' @param sample_field Some algoriths will save the name of the
#' sample in the ##SAMPLE portion of header in the VCF (e.g.
#' ##SAMPLE=<ID=TUMOR,SampleName=TCGA-01-0001>). If the ID is specified via the
#' \code{id} parameter ("TUMOR" in this example), then \code{sample_field} can
#' be used to specify the name of the tag ("SampleName" in this example).
#' Default \code{NULL}.
#' @param filter Exclude variants that do not have a \code{PASS} in the
#' \code{FILTER} column of the VCF. Default \code{TRUE}.
#' @param multiallele Multialleles are when multiple alternative variants
#' are listed in the same row in the vcf. One of \code{"expand"} or
#' \code{"exclude"}. If \code{"expand"} is selected, then each
#' alternate allele carried by the sample gets its own row: with ALT
#' \code{A,T}, a GT of \code{0/2} gives only the \code{T} allele. If
#' \code{"exclude"} is selected, then these rows will be removed; in an
#' \link[VariantAnnotation:VCF-class]{ExpandedVCF}, rows that share a
#' chromosome, position and reference allele are removed. An
#' \link[VariantAnnotation:VCF-class]{ExpandedVCF} no longer records which
#' allele each row is, so every row with a non-reference GT is kept.
#' Default \code{"expand"}.
#' @param extra_fields Optionally extract additional fields from the \code{INFO}
#' section of the VCF. Fields with one value per alternate allele
#' (\code{Number=A}) get the value for each allele; other fields with
#' several values are joined with commas. Default \code{NULL}.
#' @return Returns a data.table of variants from a vcf
#' @examples
#' vcf_file <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "musicatk"
#' )
#'
#' library(VariantAnnotation)
#' vcf <- readVcf(vcf_file)
#' variants <- extract_variants_from_vcf(vcf = vcf)
#' @export
extract_variants_from_vcf <- function(vcf, id = NULL, rename = NULL,
                                      sample_field = NULL, filter = TRUE,
                                      multiallele = c("expand", "exclude"),
                                      extra_fields = NULL) {
  multiallele <- match.arg(multiallele)

  # Configure ID
  vcf_samples <- colnames(vcf)
  id_given <- !is.null(id)
  if (length(vcf_samples) == 0) {
    if (!is.null(id)) {
      stop(
        "The VCF has no sample columns, so the id '", id, "' can't be ",
        "selected. Use 'rename' to name the sample instead."
      )
    }
    if (is.null(rename)) {
      stop(
        "The VCF has no sample columns. Use 'rename' to give a name to ",
        "the sample these variants belong to."
      )
    }
  } else if (is.null(id)) {
    id <- vcf_samples[1]
  } else {
    if (!(id %in% vcf_samples)) {
      stop(
        "The value for id, '", id, "', was not found among the list of,",
        " samples in the VCF: ",
        paste(vcf_samples, collapse = ", ")
      )
    }
  }

  # Configure sample name
  vcf_name <- id
  if (!is.null(rename)) {
    vcf_name <- rename
  } else if (!is.null(sample_field)) {
    df <- VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE
    if (!(id %in% rownames(df))) {
      stop(
        "The value for id, '", id, "' was not found in the 'SAMPLE' ",
        "metadata within the VCF."
      )
    }
    if (!sample_field %in% colnames(df)) {
      stop(
        "The value for header_sample_name_field, '",
        sample_field, "' was not found ",
        "in the 'SAMPLE' metadata within the VCF."
      )
    }
    vcf_name <- as.character(df[id, sample_field])
  }

  # Perform filtering based on FILTER column
  if (isTRUE(filter)) {
    # Remove filtered rows
    pass <- which(SummarizedExperiment::rowRanges(vcf)$FILTER == "PASS")
    if (length(pass) == 0) {
      warning(
        "No variants passed the filter for VCF with the id/name: ",
        id, "/", vcf_name
      )
      return(NULL)
    }
    vcf <- vcf[pass, ]
  }

  # Genotypes of the selected sample. Without them, every variant is kept.
  gt <- NULL
  if (!is.null(id) && "GT" %in% names(VariantAnnotation::geno(vcf))) {
    gt <- VariantAnnotation::geno(vcf)$GT[, id]
  } else {
    # Without GT the sample can't be chosen from the data
    if (length(vcf_samples) > 1 && !id_given && is.null(rename)) {
      stop(
        "The VCF has no GT field and several sample columns (",
        paste(vcf_samples, collapse = ", "), "). Use 'id' or 'rename' ",
        "to say which sample the variants belong to."
      )
    }
    message(
      "The VCF has no GT field for sample '", vcf_name, "'. All variants ",
      "passing the filter are assigned to it."
    )
  }

  if (is(vcf, "CollapsedVCF")) {
    # One row per ALT allele. Keep only the alleles listed in the sample's
    # GT (0/2 keeps the second ALT allele), since expand() copies the GT
    # unchanged to every allele of the site.
    num_alleles <- lengths(VariantAnnotation::alt(vcf))
    if (multiallele == "exclude") {
      single <- which(num_alleles == 1)
      vcf <- vcf[single, ]
      gt <- gt[single]
      num_alleles <- num_alleles[single]
    }
    site <- rep(seq_along(num_alleles), num_alleles)
    allele <- sequence(num_alleles)
    vcf <- VariantAnnotation::expand(vcf)
    if (is.null(gt)) {
      keep <- rep(TRUE, length(site))
    } else {
      gt_alleles <- strsplit(gt, "[/|]")
      carried <- paste(
        rep(seq_along(gt_alleles), lengths(gt_alleles)),
        unlist(gt_alleles)
      )
      keep <- paste(site, allele) %in% carried
    }
  } else {
    # An ExpandedVCF no longer records which allele each row is, so any
    # non-reference GT keeps the row
    if (multiallele == "exclude") {
      rows <- SummarizedExperiment::rowRanges(vcf)
      site <- paste(
        GenomicRanges::seqnames(rows), IRanges::start(rows), rows$REF
      )
      single <- which(!site %in% site[duplicated(site)])
      vcf <- vcf[single, ]
      gt <- gt[single]
    }
    if (is.null(gt)) {
      keep <- rep(TRUE, nrow(vcf))
    } else {
      keep <- stringr::str_detect(gt, "[1-9]") %in% TRUE
    }
  }

  if (sum(keep) == 0) {
    warning(
      "No variants were found for sample: ", vcf_name,
      if (!is.null(id) && !identical(id, vcf_name)) paste0(" (", id, ")"),
      ". The VCF has no variant rows left after filtering, or all ",
      "genotypes matched the reference allele."
    )
    return(NULL)
  }
  vcf <- vcf[which(keep), ]

  rows <- SummarizedExperiment::rowRanges(vcf)
  dt <- data.table::data.table(
    chr = as.character(GenomicRanges::seqnames(rows)),
    start = IRanges::start(rows),
    end = IRanges::end(rows),
    ref = as.character(rows$REF),
    alt = as.character(rows$ALT),
    sample = vcf_name
  )

  # Add extra columns if requested
  if (!is.null(extra_fields)) {
    info_field <- VariantAnnotation::info(vcf)
    temp <- setdiff(extra_fields, colnames(info_field))
    if (length(temp) > 0) {
      warning(
        "Values in extra_fields were not found in the INFO field within",
        " the VCF: ", paste(temp, collapse = ", ")
      )
    }
    temp <- intersect(extra_fields, colnames(info_field))
    for (field in temp) {
      # Fields with several values per row are joined with commas
      value <- info_field[[field]]
      if (is(value, "List")) {
        value <- S4Vectors::unstrsplit(
          methods::as(value, "CharacterList"),
          sep = ","
        )
      }
      data.table::set(dt, j = field, value = unname(value))
    }
  }
  return(dt)
}

#' Extracts variants from a vcf file
#'
#' Add Description
#'
#' @param vcf_file Path to the vcf file
#' @param id ID of the sample to select from VCF. If \code{NULL}, then the
#' first sample will be selected. Default \code{NULL}.
#' @param rename Rename the sample to this value when extracting variants.
#' If \code{NULL}, then the sample will be named according to \code{ID}.
#' @param sample_field Some algoriths will save the name of the
#' sample in the ##SAMPLE portion of header in the VCF (e.g.
#' ##SAMPLE=<ID=TUMOR,SampleName=TCGA-01-0001>). If the ID is specified via the
#' \code{id} parameter ("TUMOR" in this example), then \code{sample_field} can
#' be used to specify the name of the tag ("SampleName" in this example).
#' Default \code{NULL}.
#' @param filename_as_id If set to \code{TRUE}, the file name will be used
#' as the sample name.
#' @param strip_extension Only used if \code{filename_as_id} is set to
#' \code{TRUE}. If set to \code{TRUE}, the file extention will be stripped
#' from the filename before setting the sample name.
#' If a character vector is given, then all the strings
#' in the vector will removed from the end of the filename before setting the
#' sample name. Default \code{c(".vcf",".vcf.gz",".gz")}
#' @param filter Exclude variants that do not have a \code{PASS} in the
#' \code{FILTER} column of the VCF. Default \code{TRUE}.
#' @param multiallele Multialleles are when multiple alternative variants
#' are listed in the same row in the vcf. One of \code{"expand"} or
#' \code{"exclude"}. If \code{"expand"} is selected, then each
#' alternate allele carried by the sample gets its own row: with ALT
#' \code{A,T}, a GT of \code{0/2} gives only the \code{T} allele. If
#' \code{"exclude"} is selected, then these rows will be removed; in an
#' \link[VariantAnnotation:VCF-class]{ExpandedVCF}, rows that share a
#' chromosome, position and reference allele are removed. An
#' \link[VariantAnnotation:VCF-class]{ExpandedVCF} no longer records which
#' allele each row is, so every row with a non-reference GT is kept.
#' Default \code{"expand"}.
#' @param extra_fields Optionally extract additional fields from the \code{INFO}
#' section of the VCF. Fields with one value per alternate allele
#' (\code{Number=A}) get the value for each allele; other fields with
#' several values are joined with commas. Default \code{NULL}.
#' @param fix_vcf_errors Attempt to automatically fix VCF file
#' formatting errors. If the file can't be read, lines with a REF or ALT of
#' \code{0}, stray header lines, and lines with the wrong number of fields
#' are removed, with a warning, and the rest is read again.
#' Default \code{TRUE}.
#' @return Returns a data.table of variants extracted from a vcf
#' @examples
#' vcf <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "musicatk"
#' )
#' variants <- extract_variants_from_vcf_file(vcf_file = vcf)
#' @export
extract_variants_from_vcf_file <- function(vcf_file, id = NULL, rename = NULL,
                                           sample_field = NULL,
                                           filename_as_id = FALSE,
                                           strip_extension = c(
                                             ".vcf",
                                             ".vcf.gz",
                                             ".gz"
                                           ),
                                           filter = TRUE,
                                           multiallele = c("expand", "exclude"),
                                           extra_fields = NULL,
                                           fix_vcf_errors = TRUE) {
  vcf <- try(VariantAnnotation::readVcf(vcf_file), silent = TRUE)

  # Extract the filebase name and use it as the sample name
  if (isTRUE(filename_as_id)) {
    rename <- basename(vcf_file)
    if (!is.null(strip_extension)) {
      if (isTRUE(strip_extension)) {
        rename <- tools::file_path_sans_ext(rename)
      } else {
        rename <- gsub(
          paste(paste0(strip_extension, "$"), collapse = "|"),
          "", rename
        )
      }
    }
  }

  # Automatically try to fix some types of errors in VCF files: remove body
  # lines with a REF or ALT of 0, stray header lines, or the wrong number of
  # fields, then read the remaining lines again
  if (is(vcf, "try-error")) {
    if (!isTRUE(fix_vcf_errors)) {
      stop(
        "VCF File: ", vcf_file,
        " is malformed and automatic error fixing is disabled.",
        " \nAdditional information: \n", vcf[1]
      )
    }
    lines <- readLines(vcf_file)
    column_line <- which(startsWith(lines, "#CHROM"))[1]
    if (is.na(column_line)) {
      stop(
        "VCF File: ", vcf_file, " could not be recovered, please review.",
        " It has no '#CHROM' header line.",
        " \nAdditional information: \n", vcf[1]
      )
    }
    body <- lines[-seq_len(column_line)]
    fields <- strsplit(body, "\t", fixed = TRUE)
    num_fields <- lengths(strsplit(lines[column_line], "\t", fixed = TRUE))
    malformed <- startsWith(body, "#") |
      lengths(fields) != num_fields |
      vapply(fields, "[", character(1), 4) %in% "0" |
      vapply(fields, "[", character(1), 5) %in% "0"
    if (!any(malformed)) {
      stop(
        "VCF File: ", vcf_file, " could not be recovered, please review.",
        " \nAdditional information: \n", vcf[1]
      )
    }
    repaired_file <- tempfile(fileext = ".vcf")
    on.exit(unlink(repaired_file), add = TRUE)
    writeLines(c(lines[seq_len(column_line)], body[!malformed]), repaired_file)
    vcf <- try(VariantAnnotation::readVcf(repaired_file), silent = TRUE)
    if (is(vcf, "try-error")) {
      stop(
        "VCF File: ", vcf_file, " could not be recovered, please review.",
        " \nAdditional information: \n", vcf[1]
      )
    }
    warning(
      "VCF File: ", vcf_file, " is malformed. ", sum(malformed),
      " malformed lines were removed before extracting variants."
    )
  }

  dt <- extract_variants_from_vcf(
    vcf = vcf, id = id, rename = rename,
    sample_field = sample_field,
    multiallele = multiallele,
    filter = filter,
    extra_fields = extra_fields
  )
  return(dt)
}

#' Extract variants from a maf object
#'
#' Add description
#'
#' @param maf MAF object loaded by read.maf() from the 'maftools' package
#' @param extra_fields Optionally extract additional columns from the
#' maf object. Default \code{NULL}.
#' @return Returns a data.table of variants from a maf which can be used to
#' create a \code{musica} object.
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#'   package = "musicatk"
#' )
#' library(maftools)
#' maf <- read.maf(maf_file)
#' variants <- extract_variants_from_maf(maf = maf)
#' @export
extract_variants_from_maf <- function(maf, extra_fields = NULL) {
  if (!inherits(maf, "MAF")) {
    stop("'maf' needs to be a 'MAF' object created by the 'read.maf' function
         from package 'maftools'")
  }
  dt <- rbind(maf@data, maf@maf.silent)
  dt <- extract_variants_from_matrix(
    mat = dt, chromosome_col = "Chromosome",
    start_col = "Start_Position",
    end_col = "End_Position",
    ref_col = "Tumor_Seq_Allele1",
    alt_col = "Tumor_Seq_Allele2",
    sample_col = "Tumor_Sample_Barcode",
    extra_fields = extra_fields
  )
  return(dt)
}


#' Extract variants from matrix or data.frame like objects
#'
#' Add Description
#'
#' @param mat An object that inherits from classes "matrix" or "data.frame"
#' Examples include a matrix, data.frame, or data.table.
#' @param chromosome_col The name of the column that contains the chromosome
#' reference for each variant. Default \code{"Chromosome"}.
#' @param start_col The name of the column that contains the start
#' position for each variant. Default \code{"Start_Position"}.
#' @param end_col The name of the column that contains the end
#' position for each variant. Default \code{"End_Position"}.
#' @param ref_col The name of the column that contains the reference
#' base(s) for each variant. Default \code{"Tumor_Seq_Allele1"}.
#' @param alt_col The name of the column that contains the alternative
#' base(s) for each variant. Default \code{"Tumor_Seq_Allele2"}.
#' @param sample_col The name of the column that contains the sample
#' id for each variant. Default \code{"Tumor_Sample_Barcode"}.
#' @param extra_fields Optionally extract additional columns from the
#' object. Default \code{NULL}.
#' @return Returns a data.table of variants from a maf which can be used to
#' create a \code{musica} object.
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#'   package = "musicatk"
#' )
#' library(maftools)
#' maf <- read.maf(maf_file)
#' variants <- extract_variants_from_maf(maf = maf)
#' variants <- extract_variants_from_matrix(
#'   mat = variants,
#'   chromosome_col = "chr", start_col = "start", end_col = "end",
#'   ref_col = "ref", alt_col = "alt", sample_col = "sample"
#' )
#' @export
extract_variants_from_matrix <- function(mat, chromosome_col = "Chromosome",
                                         start_col = "Start_Position",
                                         end_col = "End_Position",
                                         ref_col = "Tumor_Seq_Allele1",
                                         alt_col = "Tumor_Seq_Allele2",
                                         sample_col = "Tumor_Sample_Barcode",
                                         extra_fields = NULL) {
  if (!inherits(mat, c("matrix", "data.frame"))) {
    stop("'mat' needs to inherit classes 'matrix' or 'data.frame'")
  }
  dt <- .check_headers(data.table::as.data.table(mat),
    chromosome = chromosome_col,
    start = start_col,
    end = end_col,
    ref = ref_col,
    alt = alt_col,
    sample = sample_col
  )

  # Add extra columns if requested
  if (!is.null(extra_fields)) {
    temp <- setdiff(extra_fields, colnames(dt))
    if (length(temp) > 0) {
      warning(
        "Values in 'extra_fields' were not found in column names in",
        " the MAF object: ", paste(temp, collapse = ", ")
      )
    }
    extra_fields <- intersect(extra_fields, colnames(dt))
  }

  all_fields <- c(.required_musica_headers(), extra_fields)
  dt <- dt[, all_fields, with = FALSE]
  return(dt)
}


#' Extracts variants from a maf file
#'
#' Add Description - Aaron
#'
#' @param maf_file Location of maf file
#' @param extra_fields Optionally extract additional columns from the
#' object. Default \code{NULL}.
#' @return Returns a data.table of variants from a maf, or \code{NULL}
#' with a warning if the file has no variant records. Files whose
#' records are all silent or non-coding are read in full.
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#'   package = "musicatk"
#' )
#' maf <- extract_variants_from_maf_file(maf_file = maf_file)
#' @export
extract_variants_from_maf_file <- function(maf_file, extra_fields = NULL) {
  if (!.maf_has_records(maf_file)) {
    warning("MAF file ", maf_file, " has no variant records.")
    return(NULL)
  }
  maf <- tryCatch(
    maftools::read.maf(maf_file, verbose = FALSE),
    error = function(e) {
      # read.maf() stops when no record is non-synonymous. Reading again
      # with every class as non-synonymous keeps silent and non-coding
      # records; any other error is raised again by this call.
      maftools::read.maf(maf_file,
        vc_nonSyn = .maf_variant_classes(), verbose = FALSE
      )
    }
  )
  return(extract_variants_from_maf(maf = maf, extra_fields = extra_fields))
}

# Whether a MAF file has at least one record after its comment and column
# header lines. Reads only as far as the first record.
.maf_has_records <- function(maf_file) {
  con <- file(maf_file, open = "r")
  on.exit(close(con))
  lines_seen <- 0L
  repeat {
    lines <- readLines(con, n = 1000L)
    if (length(lines) == 0) {
      return(FALSE)
    }
    lines_seen <- lines_seen +
      sum(!startsWith(lines, "#") & nzchar(trimws(lines)))
    if (lines_seen >= 2L) {
      return(TRUE)
    }
  }
}

# The Variant_Classification values of the MAF specification
.maf_variant_classes <- function() {
  c(
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Silent", "Splice_Site", "Splice_Region", "Translation_Start_Site",
    "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA",
    "Targeted_Region", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
    "Start_Codon_SNP", "Start_Codon_Del", "Start_Codon_Ins",
    "Stop_Codon_Del", "Stop_Codon_Ins", "lincRNA", "Unknown"
  )
}

#' Creates a musica object from a variant table
#'
#' This function creates a \linkS4class{musica} object from a variant
#' table or matrix. The \linkS4class{musica} class stores variants information,
#' variant-level annotations, sample-level annotations, and count tables and
#' is used as input to the mutational signature discovery and prediction
#' algorithms. The input variant table or matrix must have columns for
#' chromosome, start position, end position, reference allele,
#' alternate allele, and sample names. The column names in the variant table
#' can be mapped using the \code{chromosome_col}, \code{start_col},
#' \code{end_col}, \code{ref_col}, \code{alt_col}, and
#' \code{sample_col parameters}.
#'
#' @param x A data.table, matrix, or data.frame that contains columns with
#' the variant information.
#' @param genome A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from.
#' @param check_ref_chromosomes Whether to peform a check to ensure that
#' the chromosomes in the \code{variant} object match the reference
#' chromosomes in the \code{genome} object. If there are mismatches, this
#' may cause errors in downstream generation of count tables. If mismatches
#' occur, an attept to be automatically fix these with the
#' \code{\link[GenomeInfoDb]{seqlevelsStyle}} function will be made.
#' Default \code{TRUE}.
#' @param check_ref_bases Whether to check if the reference bases in the
#' \code{variant} object match the reference bases in the \code{genome}
#' object. Default \code{TRUE}.
#' @param chromosome_col The name of the column that contains the chromosome
#' reference for each variant. Default \code{"chr"}.
#' @param start_col The name of the column that contains the start
#' position for each variant. Default \code{"start"}.
#' @param end_col The name of the column that contains the end
#' position for each variant. Default \code{"end"}.
#' @param ref_col The name of the column that contains the reference
#' base(s) for each variant. Default \code{"ref"}.
#' @param alt_col The name of the column that contains the alternative
#' base(s) for each variant. Default \code{"alt"}.
#' @param sample_col The name of the column that contains the sample
#' id for each variant. Default \code{"sample"}.
#' @param extra_fields Which additional fields to extract and include in
#' the musica object. Default \code{NULL}.
#' @param standardize_indels Flag to convert indel style (e.g. `C > CAT`
#' becomes `- > AT` and `GCACA > G` becomes `CACA > -`). Deletions are moved
#' to start at the first deleted base, as in MAF files; insertions keep
#' \code{start} at the base before the insertion. Indels whose first base
#' differs between \code{ref} and \code{alt} (e.g. `A > GGT`) are complex
#' and are removed. IND83 tables need standardized indels.
#' Default \code{TRUE}.
#' @param convert_dbs Flag to convert adjacent SBS in the same sample into
#' DBS (original SBS are removed). A run of three or more adjacent SBS
#' becomes one variant of type \code{"MNV"}, which is not included in the
#' standard count tables. Default \code{TRUE}.
#' @param verbose Whether to print status messages during error checking.
#' Default \code{TRUE}.
#' @return Returns a musica object
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#'   package = "musicatk"
#' )
#' variants <- extract_variants_from_maf_file(maf_file)
#' g <- select_genome("38")
#' musica <- create_musica_from_variants(x = variants, genome = g)
#' @export
create_musica_from_variants <- function(x, genome,
                                        check_ref_chromosomes = TRUE,
                                        check_ref_bases = TRUE,
                                        chromosome_col = "chr",
                                        start_col = "start",
                                        end_col = "end",
                                        ref_col = "ref",
                                        alt_col = "alt",
                                        sample_col = "sample",
                                        extra_fields = NULL,
                                        standardize_indels = TRUE,
                                        convert_dbs = TRUE,
                                        verbose = TRUE) {
  used_fields <- c(.required_musica_headers(), extra_fields)
  if (inherits(x, c("matrix", "data.frame")) || canCoerce(x, "data.table")) {
    dt <- data.table::as.data.table(x)
  } else {
    stop(
      "'x' needs to be an object which can be coerced to a data.table. ",
      "Valid classes include but are not limited to 'matrix', 'data.frame'",
      " and 'data.table'."
    )
  }
  if (!inherits(genome, "BSgenome")) {
    stop(
      "'genome' needs to be a 'BSgenome' object containing the genome ",
      "reference that was used when calling the variants."
    )
  }

  # Check for necessary columns and change column names to stardard object
  dt <- .check_headers(dt,
    chromosome = chromosome_col,
    start = start_col,
    end = end_col,
    ref = ref_col,
    alt = alt_col,
    sample = sample_col,
    update_fields = TRUE
  )

  # Subset to necessary columns and add variant type
  all_fields <- c(.required_musica_headers(), extra_fields)
  dt <- dt[, all_fields, with = FALSE]

  # Columns may be factors (data.frame) or text (character matrix), and
  # bases may be lowercase
  for (col in c("chr", "sample")) {
    data.table::set(dt, j = col, value = as.character(dt[[col]]))
  }
  for (col in c("ref", "alt")) {
    data.table::set(dt, j = col, value = toupper(as.character(dt[[col]])))
  }
  for (col in c("start", "end")) {
    data.table::set(dt, j = col, value = as.integer(as.character(dt[[col]])))
  }
  dt <- add_variant_type(dt)

  # Some non-variants are included (e.g. T>T). These will be removed
  non_variant <- which(dt$ref == dt$alt)
  if (length(non_variant) > 0) {
    warning(
      length(non_variant), " variants has the same reference and ",
      "alternate allele. These variants were excluded."
    )
    dt <- dt[-non_variant, ]
  }

  if (isTRUE(check_ref_chromosomes)) {
    # Check for genome style and attempt to convert variants to reference
    # genome if they don't match
    if (isTRUE(verbose)) {
      message(
        "Checking that chromosomes in the 'variant' object match ",
        "chromosomes in the 'genome' object."
      )
    }
    dt <- .check_variant_genome(dt = dt, genome = genome)
  }

  # Variants on contigs that aren't in the genome can't be checked or
  # counted
  missing_chr <- !dt$chr %in% GenomeInfoDb::seqnames(genome)
  if (any(missing_chr)) {
    warning(
      sum(missing_chr), " variants are on chromosomes not found in the ",
      "'genome' object and were excluded: ",
      paste(unique(dt$chr[missing_chr]), collapse = ", ")
    )
    dt <- dt[!missing_chr, ]
  }
  if (nrow(dt) == 0) {
    stop(
      "No variants are left after removing those on chromosomes not found ",
      "in the 'genome' object. Check that the chromosome names match ",
      "'seqnames(genome)' or set 'check_ref_chromosomes = TRUE'."
    )
  }

  if (isTRUE(check_ref_bases)) {
    if (isTRUE(verbose)) {
      message(
        "Checking that the reference bases in the 'variant' object ",
        "match the reference bases in the 'genome' object."
      )
    }
    .check_variant_ref_in_genome(dt = dt, genome = genome)
  }

  if (isTRUE(standardize_indels)) {
    if (isTRUE(verbose)) {
      message("Standardizing INS/DEL style")
    }
    bases <- c("A", "T", "G", "C")
    comp_ins <- which(dt$Variant_Type == "INS" & !dt$ref %in% c(bases, "-"))
    if (length(comp_ins) > 0) {
      message("Removing ", length(comp_ins), " compound insertions")
      dt <- dt[-comp_ins, ]
    }

    comp_del <- which(dt$Variant_Type == "DEL" & !dt$alt %in% c(bases, "-"))
    if (length(comp_del) > 0) {
      message("Removing ", length(comp_del), " compound deletions")
      dt <- dt[-comp_del, ]
    }

    # VCF-style indels start with an anchor base shared by ref and alt (the
    # rule VariantAnnotation's isInsertion() and isDeletion() use). When the
    # first bases differ, the variant is complex, not an indel.
    ins <- dt$Variant_Type == "INS" & dt$ref %in% bases
    del <- dt$Variant_Type == "DEL" & dt$alt %in% bases
    complex <- which(
      (ins & substr(dt$alt, 1, 1) != dt$ref) |
        (del & substr(dt$ref, 1, 1) != dt$alt)
    )
    if (length(complex) > 0) {
      message(
        "Removing ", length(complex), " complex variants whose first ",
        "base differs between the reference and alternate alleles"
      )
      dt <- dt[-complex, ]
    }

    # Remove the anchor base. Insertions keep 'start' as the base before
    # the insertion; deletions move 'start' to the first deleted base, so
    # 'start' to 'end' covers the deleted bases as in MAF files.
    ins <- which(dt$Variant_Type == "INS" & dt$ref %in% bases)
    if (length(ins) > 0) {
      message("Converting ", length(ins), " insertions")
      dt$ref[ins] <- "-"
      ins_alt <- dt$alt[ins]
      dt$alt[ins] <- substr(ins_alt, 2, nchar(ins_alt))
    }

    del <- which(dt$Variant_Type == "DEL" & dt$alt %in% bases)
    if (length(del) > 0) {
      message("Converting ", length(del), " deletions")
      dt$alt[del] <- "-"
      del_ref <- dt$ref[del]
      dt$ref[del] <- substr(del_ref, 2, nchar(del_ref))
      dt$start[del] <- dt$start[del] + 1L
      dt$end[del] <- dt$start[del] + nchar(dt$ref[del]) - 1L
    }
  }

  if (isTRUE(convert_dbs)) {
    if (isTRUE(verbose)) {
      message("Converting adjacent SBS into DBS")
    }
    # Find runs of SBS at adjacent positions in the same sample and
    # chromosome. Sorting by position first keeps each run's rows in
    # position order.
    sbs <- which(dt$Variant_Type == "SBS")
    sbs <- sbs[order(dt$start[sbs])]
    sbs_ranges <- GenomicRanges::GRanges(
      seqnames = paste(dt$sample[sbs], dt$chr[sbs], sep = "\t"),
      ranges = IRanges::IRanges(start = dt$start[sbs], width = 1)
    )
    runs <- GenomicRanges::reduce(sbs_ranges, with.revmap = TRUE)
    revmap <- S4Vectors::mcols(runs)$revmap
    # Skip runs with two SBS at one position
    is_run <- lengths(revmap) > 1 & IRanges::width(runs) == lengths(revmap)
    revmap <- revmap[is_run]
    if (length(revmap) > 0) {
      run_id <- rep(seq_along(revmap), lengths(revmap))
      run_rows <- sbs[unlist(revmap)]
      # Each run becomes one variant in the row of its first SBS in the
      # table: a DBS for two bases, an MNV for three or more
      keep <- min(S4Vectors::splitAsList(run_rows, run_id))
      run_size <- lengths(revmap)
      dt$start[keep] <- min(
        S4Vectors::splitAsList(dt$start[run_rows], run_id)
      )
      dt$end[keep] <- dt$start[keep] + run_size - 1L
      dt$ref[keep] <- S4Vectors::unstrsplit(
        S4Vectors::splitAsList(dt$ref[run_rows], run_id)
      )
      dt$alt[keep] <- S4Vectors::unstrsplit(
        S4Vectors::splitAsList(dt$alt[run_rows], run_id)
      )
      dt$Variant_Type[keep] <- ifelse(run_size == 2, "DBS", "MNV")
      num_dbs <- sum(run_size == 2)
      num_mnv <- sum(run_size > 2)
      if (num_dbs > 0) {
        message(num_dbs, " pairs of adjacent SBS converted to DBS")
      }
      if (num_mnv > 0) {
        message(
          num_mnv, " runs of three or more adjacent SBS converted to MNV"
        )
      }
      dt <- dt[-setdiff(run_rows, keep), ]
    }
  }

  # Create and return a musica object
  s <- gtools::mixedsort(unique(dt$sample))
  annot <- data.frame(Samples = factor(s, levels = s))
  dt$sample <- factor(dt$sample, levels = s)

  musica <- new("musica", variants = dt, sample_annotations = annot)

  return(musica)
}

#' Creates a musica object from a mutation count table
#'
#' This function creates a \linkS4class{musica} object from a mutation count
#' table or matrix. The \linkS4class{musica} class stores variants information,
#' variant-level annotations, sample-level annotations, and count tables and
#' is used as input to the mutational signature discovery and prediction
#' algorithms.
#'
#' @param x A data.table, matrix, or data.frame that contains counts of mutation
#' types for each sample, with samples as columns. The column names are used
#' as the sample names and must be unique. If \code{x} has row names, they
#' must be the standard motif names (e.g. \code{"C>A_ACA"} for SBS96 or
#' \code{"AC>NN_CA"} for DBS78) or the COSMIC/SigProfiler names
#' (e.g. \code{"A[C>A]A"} or \code{"AC>CA"}), and the rows are matched to the
#' motifs by name. Without row names, the rows must be in the standard motif
#' order.
#' @param variant_class Mutations are SBS, DBS, or Indel.
#' @return Returns a musica object
#' @examples
#' data(musica)
#' count_table <- get_count_table(extract_count_tables(musica)$SBS96)
#' musica <- create_musica_from_counts(count_table, "SBS96")
#' @export
create_musica_from_counts <- function(x, variant_class) {
  if (canCoerce(x, "matrix")) {
    x <- as.matrix(x)
  } else {
    stop("'count_table' needs to be an object which can be coerced to a 
         matrix. ")
  }

  if (is.null(colnames(x))) {
    stop("'x' must have column names giving the sample names.")
  }
  if (anyDuplicated(colnames(x)) > 0 || anyNA(colnames(x))) {
    stop("The column names (sample names) of 'x' must be unique.")
  }

  # create musica object with one annotation row per sample
  musica <- new("musica",
    sample_annotations = data.frame(Samples = colnames(x))
  )

  if (variant_class %in% c("snv", "SNV", "SNV96", "SBS", "SBS96")) {
    if (nrow(x) != 96) {
      stop("SBS96 'count_table' must have 96 rows.")
    }

    # create SBS mutation type list
    forward_change <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    b1 <- rep(rep(c("A", "C", "G", "T"), each = 4), 6)
    b2 <- rep(c("C", "T"), each = 48)
    b3 <- rep(c("A", "C", "G", "T"), 24)
    mut_trinuc <- apply(cbind(b1, b2, b3), 1, paste, collapse = "")
    mut <- rep(forward_change, each = 16)
    annotation <- data.frame(
      "motif" = paste0(mut, "_", mut_trinuc),
      "mutation" = mut,
      "context" = mut_trinuc
    )
    rownames(annotation) <- annotation$motif

    # color mapping for mutation types
    color_mapping <- c(
      "C>A" = "#5ABCEBFF",
      "C>G" = "#050708FF",
      "C>T" = "#D33C32FF",
      "T>A" = "#CBCACBFF",
      "T>C" = "#ABCD72FF",
      "T>G" = "#E7C9C6FF"
    )

    # match count table rows to the SBS96 motifs by name
    x <- .match_count_rows(x, annotation$motif, "SBS96")

    # create count table object
    tab <- new("count_table",
      name = "SBS96", count_table = x,
      annotation = annotation, features = as.data.frame(annotation$motif[1]),
      type = S4Vectors::Rle("SBS"), color_variable = "mutation",
      color_mapping = color_mapping, description = paste0(
        "Single Base Substitution table with",
        " one base upstream and downstream"
      )
    )

    # add count table to musica object
    tables(musica)[["SBS96"]] <- tab
  } else if (variant_class %in% c("DBS", "dbs", "doublet")) {
    if (nrow(x) != 78) {
      stop("DBS78 'count_table' must have 78 rows.")
    }

    full_motif <- c(
      paste0("AC>NN", "_", c(
        "CA", "CG", "CT", "GA", "GG", "GT",
        "TA", "TG", "TT"
      )),
      paste0("AT>NN", "_", c("CA", "CC", "CG", "GA", "GC", "TA")),
      paste0("CC>NN", "_", c(
        "AA", "AG", "AT", "GA", "GG", "GT", "TA",
        "TG", "TT"
      )),
      paste0("CG>NN", "_", c("AT", "GC", "GT", "TA", "TC", "TT")),
      paste0("CT>NN", "_", c(
        "AA", "AC", "AG", "GA", "GC", "GG", "TA",
        "TC", "TG"
      )),
      paste0("GC>NN", "_", c("AA", "AG", "AT", "CA", "CG", "TA")),
      paste0("TA>NN", "_", c("AT", "CG", "CT", "GC", "GG", "GT")),
      paste0("TC>NN", "_", c(
        "AA", "AG", "AT", "CA", "CG", "CT", "GA",
        "GG", "GT"
      )),
      paste0("TG>NN", "_", c(
        "AA", "AC", "AT", "CA", "CC", "CT", "GA",
        "GC", "GT"
      )),
      paste0("TT>NN", "_", c(
        "AA", "AC", "AG", "CA", "CC", "CG", "GA",
        "GC", "GG"
      ))
    )

    annotation <- data.frame(
      motif = full_motif, mutation =
        unlist(lapply(strsplit(full_motif, "_"), "[[", 1)),
      context = unlist(lapply(
        strsplit(full_motif, "_"),
        "[[", 2
      )),
      row.names = full_motif
    )

    color_mapping <- .gg_color_hue(length(unique(annotation$mutation)))
    names(color_mapping) <- unique(annotation$mutation)

    # match count table rows to the DBS78 motifs by name
    x <- .match_count_rows(x, annotation$motif, "DBS78")

    # create count table object
    tab <- new("count_table",
      name = "DBS78", count_table = x,
      annotation = annotation, features = as.data.frame(annotation$motif[1]),
      type = S4Vectors::Rle("DBS"), color_variable = "mutation",
      color_mapping = color_mapping, description = paste0(
        "Standard count table for ",
        "double base substitutions",
        "using COSMIC v3 schema"
      )
    )

    # add count table to musica object
    tables(musica)[["DBS78"]] <- tab
  } else if (variant_class %in% c(
    "INDEL", "Indel", "indel", "ind", "IND",
    "ID"
  )) {
    stop("Not yet supported.")
  } else {
    stop("Only SBS and DBS classes are supported")
  }

  return(musica)
}



# Puts the rows of a count matrix in the order of 'motifs', matching them
# by row name. COSMIC/SigProfiler names are converted to musicatk names
# first: A[C>A]G becomes C>A_ACG and AC>CA becomes AC>NN_CA. Without row
# names, the rows are assumed to be in the order of 'motifs'.
.match_count_rows <- function(x, motifs, table_name) {
  if (!is.null(rownames(x))) {
    rownames(x) <- sub(
      "^([ACGT])\\[(([ACGT])>[ACGT])\\]([ACGT])$", "\\2_\\1\\3\\4",
      rownames(x)
    )
    rownames(x) <- sub("^([ACGT]{2})>([ACGT]{2})$", "\\1>NN_\\2", rownames(x))
  }
  if (is.null(rownames(x))) {
    message(
      "'x' has no row names, so its rows are assumed to be in the ",
      "standard ", table_name, " order (", motifs[1], ", ", motifs[2],
      ", ...)."
    )
    rownames(x) <- motifs
    return(x)
  }
  unmatched <- setdiff(rownames(x), motifs)
  duplicated_names <- unique(rownames(x)[duplicated(rownames(x))])
  if (length(unmatched) > 0 || length(duplicated_names) > 0) {
    stop(
      "The row names of 'x' must be the ", length(motifs), " standard ",
      table_name, " motif names, such as '", motifs[1], "'. ",
      "Not matched: ", paste(utils::head(unmatched, 5), collapse = ", "),
      ". Duplicated: ",
      paste(utils::head(duplicated_names, 5), collapse = ", "),
      ". Missing: ",
      paste(utils::head(setdiff(motifs, rownames(x)), 5), collapse = ", ")
    )
  }
  return(x[motifs, , drop = FALSE])
}


#' Load an external model into a musica object
#'
#' This function creates a \linkS4class{result_model} object from signatures,
#' exposures, and a mutation count table, and stores it in a provided musica
#' object.
#'
#' @param signatures A matrix or data.frame of signatures by mutational motifs
#' @param exposures A matrix or data.frame of signature weights, with
#' signatures as rows and samples as columns. The column names must be sample
#' names from the \code{musica} object (all or some of them); columns are
#' matched by name and put in \code{\link{sample_names}} order.
#' @param musica Existing musica object to add the new model to
#' @param result_name Name of result list entry to store the model in
#' @param modality Modality of the model
#' @param model_id Name of model
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @return A \linkS4class{musica} object
#' @examples
#' signatures <- signatures(res, "result", "SBS96", "res")
#' exposures <- exposures(res, "result", "SBS96", "res")
#' add_result(signatures, exposures, musica = musica_annot,
#' result_name = "result", modality = "SBS96", model_id = "example_model")
#' @export
add_result <- function(signatures, exposures, musica, result_name, 
                       modality, model_id, make_copy = FALSE) {
  
  # update global variable
  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }
  
  # Put exposure columns in sample order, matched by sample name
  exposures <- as.matrix(exposures)
  samples <- sample_names(musica)
  if (!is.null(colnames(exposures)) && all(colnames(exposures) %in% samples)) {
    exposures <- exposures[, samples[samples %in% colnames(exposures)],
      drop = FALSE
    ]
  }

  # create musica result object with given exposures and signatures
  model_result <- new("result_model",
    signatures = as.matrix(signatures),
    exposures = exposures,
    num_signatures = dim(signatures)[2],
    model_id = model_id, modality = modality
  )
  
  # if this result name does not exist in result_list
  if (is.null(result_list(musica)[[result_name]])) {
    # make new list entry with the desired name and assign a new
    # result_collection object
    musica@result_list[[result_name]] <- new("result_collection",
                                             modality = SimpleList(),
                                             parameter = list(), hyperparameter = list()
    )
  }
  
  # if desired modality does not exist, create entry
  if (is.null(get_modality(musica, result_name, modality))) {
    musica@result_list[[result_name]]@modality[[modality]] <- list()
  }
  
  # add result_model object in the list for the proper modality
  musica@result_list[[result_name]]@modality[[modality]][[model_id]] <- model_result
  # exposures must have the sample names as columns, in sample order
  methods::validObject(musica)
  
  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }
  
  if (make_copy == TRUE) {
    return(musica)
  }
}



.check_variant_genome <- function(dt, genome) {
  chr_header <- .required_musica_headers()["chromosome"]
  if (!chr_header %in% colnames(dt)) {
    stop("The column '", chr_header, "' was not found in the data.table.")
  }
  chr <- as.data.frame(dt[, chr_header, with = FALSE])[, 1]
  chr_u <- unique(chr)
  genome_u <- unique(GenomeInfoDb::seqnames(genome))
  diff <- setdiff(chr_u, genome_u)

  if (length(diff) > 0) {
    # Try to use GenomeInfoDb to determine style of variants and genome
    g_error <- try(
      genome_style <-
        GenomeInfoDb::seqlevelsStyle(genome_u)[1],
      silent = TRUE
    )
    v_error <- try(
      variant_style <-
        GenomeInfoDb::seqlevelsStyle(as.character(chr_u))[1],
      silent = TRUE
    )

    inter <- intersect(chr_u, genome_u)
    if (length(inter) == 0) {
      # Error 1: No matching of genome and variants
      if (is(g_error, "try-error") && is(v_error, "try-error")) {
        stop(
          "The style of the genome references in the 'variant' and ",
          "'genome' objects did not match each other or any style ",
          "from the 'GenomeInfoDb' package. Please ensure that the entries ",
          "in the'", chr_header, "' column in the variant table match ",
          "entries in 'seqnames(genome)'. First five chromosomes:\n",
          "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
          "genome: ", paste(head(genome_u, 5), collapse = ", ")
        )
      } else if (is(g_error, "try-error")) {
        stop(
          "The style of the genome references in the 'variant' and ",
          "'genome' objects did not match each other. The style for the ",
          "variant table was determined to be ", variant_style, " by the ",
          "'GenomeInfoDb' package. However, the 'genome' object did not ",
          "match any known styles and therefore could not be automatically ",
          "mapped. Please ensure that the entries ",
          "in the'", chr_header, "' column in the variant table match ",
          "entries in 'seqnames(genome)'. First five chromosomes:\n",
          "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
          "genome: ", paste(head(genome_u, 5), collapse = ", ")
        )
      } else if (is(v_error, "try-error")) {
        stop(
          "The style of the genome references in the 'variant' and ",
          "'genome' objects did not match each other. The style for the ",
          "genome object was determined to be ", genome_style, " by the ",
          "'GenomeInfoDb' package. However, the variant table did not ",
          "match any known styles and therefore could not be automatically ",
          "mapped. Please ensure that the entries ",
          "in the'", chr_header, "' column in the variant table match ",
          "entries in 'seqnames(genome)'. First five chromosomes:\n",
          "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
          "genome: ", paste(head(genome_u, 5), collapse = ", ")
        )
      } else {
        # Attempt to map variants to genome object
        new_chr <- as.character(chr)
        map_error <- try(GenomeInfoDb::seqlevelsStyle(new_chr) <-
                           genome_style, silent = TRUE)
        if (is(map_error, "try-error")) {
          stop(
            "The style of the genome references in the 'variant' and ",
            "'genome' objects did not match each other. The style for the ",
            "genome object was determined to be ", genome_style, " and the ",
            "style of the variants was determined to be ", variant_style,
            " by the 'GenomeInfoDb' package. However, they were not able ",
            "to be automatically mapped. Please ensure that the entries ",
            "in the'", chr_header, "' column in the variant table match ",
            "entries in 'seqnames(genome)'. First five chromosomes:\n",
            "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
            "genome: ", paste(head(genome_u, 5), collapse = ", ")
          )
        } else {
          # Set new chromosomes
          dt[, eval(quote(chr_header))] <- new_chr

          # Determine if mapping was complete or partial
          new_chr_u <- unique(new_chr)
          new_diff <- setdiff(new_chr_u, genome_u)
          if (length(new_diff) > 0) {
            # If conversion was partial, need to subset to overlapping variants
            dt <- subset(dt, new_chr %in% genome_u)

            new.inter <- intersect(new_chr_u, genome_u)
            inter_sum <- sum(new_chr %in% genome_u)
            diff.sum <- sum(!new_chr %in% genome_u)
            warning(
              "The style of the genome references in the 'variant' and ",
              "'genome' objects did not match each other. The style for the ",
              "genome object was determined to be ", genome_style, " and the ",
              "style of the variants was determined to be ", variant_style,
              " by the 'GenomeInfoDb' package. ", inter_sum, " variants were ",
              "automatically converted to those in the genome object. ",
              diff.sum, " variants were not able to be converted and ",
              "were excluded from the final variant table. Variant ",
              "chromosomes that were not able to be converted: \n",
              paste(new_diff, collapse = ", ")
            )
          } else {
            warning(
              "The style of the genome references in the 'variant' and ",
              "'genome' objects did not match each other. The style for the ",
              "genome object was determined to be ", genome_style, " and the ",
              "style of the variants was determined to be ", variant_style,
              " by the 'GenomeInfoDb' package. All variant chromosomes were ",
              "automatically converted to those in the genome object."
            )
          }
        }
      }
    }
  }
  return(dt)
}

.check_variant_ref_in_genome <- function(dt, genome) {
  headers <- .required_musica_headers()
  chr <- as.character(as.data.frame(dt)[, headers["chromosome"]])
  start <- as.numeric(as.data.frame(dt)[, headers["start"]])
  end <- as.numeric(as.data.frame(dt)[, headers["end"]])
  ref <- as.character(as.data.frame(dt)[, headers["ref"]])

  genome_ref <- BSgenome::getSeq(
    x = genome,
    names = chr, start = start, end = end,
    as.character = TRUE
  )

  # Only check references that have bases. For example maf files will have "-"
  # for all insertions which cannot be checked.
  ix <- grepl("[^ACGT]", ref)
  no_match_sum <- sum(ref[!ix] != genome_ref[!ix])
  if (no_match_sum > 0) {
    warning(
      "Reference bases for ", no_match_sum, " out of ",
      length(ref), " variants did not match the reference base in the ",
      "'genome'. Make sure the genome reference is correct."
    )
  }
}
.check_headers <- function(dt, chromosome = NULL,
                           start = NULL, end = NULL,
                           ref = NULL, alt = NULL, sample = NULL,
                           update_fields = TRUE) {
  # If headers are not given in arguments, then use headers from maf as default
  rfields <- .required_maf_headers()
  if (!is.null(chromosome)) rfields["chromosome"] <- chromosome
  if (!is.null(start)) rfields["start"] <- start
  if (!is.null(end)) rfields["end"] <- end
  if (!is.null(ref)) rfields["ref"] <- ref
  if (!is.null(alt)) rfields["alt"] <- alt
  if (!is.null(sample)) rfields["sample"] <- sample

  # Check for missing columns
  col.ix <- match(tolower(rfields), tolower(colnames(dt)))
  if (sum(is.na(col.ix) > 0)) {
    missing.cols <- rfields[is.na(col.ix)]
    stop(
      "Some required columns are missing in the maf: ",
      paste(missing.cols, collapse = ", ")
    )
  }

  # Adjust column case if needed
  prev_col <- colnames(dt)[col.ix]
  colnames(dt)[col.ix] <- rfields
  mismatch_ix <- colnames(dt)[col.ix] != prev_col
  if (sum(mismatch_ix) > 0) {
    warning(
      "Some columns in maf had the wrong case and were automatically ",
      "adjusted: ", paste(prev_col[mismatch_ix], collapse = ", ")
    )
  }

  # Change columns of data.table to match the required musica format
  if (isTRUE(update_fields)) {
    data.table::setnames(dt, rfields, .required_musica_headers())
  }
  return(dt)
}

.required_musica_headers <- function() {
  return(c(
    chromosome = "chr", start = "start", end = "end", ref = "ref",
    alt = "alt", sample = "sample"
  ))
}

.required_maf_headers <- function() {
  return(c(
    chromosome = "Chromosome", start = "Start_Position",
    end = "End_Position", ref = "Tumor_Seq_Allele1",
    alt = "Tumor_Seq_Allele2", sample = "Tumor_Sample_Barcode"
  ))
}
