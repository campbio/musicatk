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
  #Choose genome build version
  #Keep for now as a helper function
  if (tolower(x) %in% c("hg19", "19", "grch37")) {
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (tolower(x) %in% c("hg38", "38", "grch38")) {
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (tolower(x) %in% c("mm10", "grcm38")) {
    g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else if (tolower(x) %in% c("mm9", "mgscv37")) {
    g <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
  } else {
    stop("That genome build is not currently selectable with this function. ",
         "Run BSgenome::available.genomes() to get a full list of availble ",
         "genomes or use BSgenome::forgeBSgenomeDataPkg to create a ",
         "custome genome package.")
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
#' package = "musicatk")
#' lusc_maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#' package = "musicatk")
#' melanoma_vcfs <- list.files(system.file("extdata", package = "musicatk"),
#'   pattern = glob2rx("*SKCM*vcf"), full.names = TRUE)
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
      stop("The length of 'rename' must be the same as the length of 'input'.",
           " Only vcf object or file.vcf inputs will be renamed.")
    }
  } else {
    rename <- rep(NULL, length(input_list))
  }
  if (!is.null(id)) {
    if (length(id) != length(input_list)) {
      stop("The lenght of 'rename' must be the same as the length of 'input'.",
           " This only applies to vcf object or file.vcf inputs.")
    }
  } else {
    id <- rep(NULL, length(input_list))
  }


  pb <- utils::txtProgressBar(min = 0, max = length(input_list), initial = 0,
                                style = 3)
  for (i in seq_along(inputs)) {
    input <- inputs[[i]]

    if (inherits(input, c("CollapsedVCF", "ExpandedVCF"))) {
      dt <- extract_variants_from_vcf(vcf = input,
                                    id = id[i],
                                    rename = rename[i],
                                    sample_field = sample_field,
                                    filter = filter,
                                    multiallele = multiallele,
                                    extra_fields = extra_fields)
    } else if (is(input, "MAF")) {
      dt <- extract_variants_from_maf(maf = input, extra_fields = extra_fields)
    } else if (inherits(input, c("matrix", "data.frame"))) {
      dt <- extract_variants_from_matrix(mat = input,
                                         extra_fields = extra_fields,
                                         chromosome_col = chromosome_col,
                                         start_col = start_col,
                                         end_col = end_col,
                                         ref_col = ref_col,
                                         alt_col = alt_col,
                                         sample_col = sample_col)
    } else if (is(input, "character")) {
      if (tools::file_ext(input) %in% c("vcf", "vcf.gz")) {
        dt <- extract_variants_from_vcf_file(vcf_file = input,
                                             id = id[i],
                                             rename = rename[i],
                                             sample_field = sample_field,
                                             filename_as_id = filename_as_id,
                                             strip_extension = strip_extension,
                                             filter = filter,
                                             multiallele = multiallele,
                                             extra_fields = extra_fields,
                                             fix_vcf_errors = fix_vcf_errors)
      } else if (tools::file_ext(input) %in% c("maf", "maf.gz")) {
        dt <- extract_variants_from_maf_file(maf_file = input,
                           extra_fields = extra_fields)
      } else {
        stop("Input file could not be automatically parsed: ", input, " ")
      }
    } else {
      stop("Each input must be a collapedVCF/expandedVCF object, MAF object ",
          "file.vcf, file.maf, or a matrix/data.frame object.",
         "Item ", input, " is of class '", class(input), "'")
    }
    input_list[[i]] <- dt

    utils::setTxtProgressBar(pb, i)
    if (isTRUE(verbose)) {
      message("Extracted ", i, " out of ", length(inputs), " inputs: ",
              input)
    }
  }

  dt <- do.call("rbind", input_list)
  return(dt)
}

#' Extracts variants from a VariantAnnotation VCF object
#'
#' Aaron - Need to describe differnce between ID, and name in the header, and rename in
#' terms of naming the sample. Need to describe differences in multiallelic
#' choices. Also need to describe the automatic error fixing
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
#' alternate allele will be given their own rows. If \code{"exclude"} is
#' selected, then these rows will be removed. Default \code{"expand"}.
#' @param extra_fields Optionally extract additional fields from the \code{INFO}
#' section of the VCF. Default \code{NULL}.
#' @return Returns a data.table of variants from a vcf
#' @examples
#' vcf_file <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "musicatk")
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

  # Process MultiAllelic Sites
  num_alleles <- lengths(VariantAnnotation::fixed(vcf)[, "ALT"])
  if (multiallele == "expand") {
    vcf <- VariantAnnotation::expand(vcf)
  } else {
    multi_allelic <- which(num_alleles != 1)
    if (length(multi_allelic) > 0) {
      vcf <- vcf[-multi_allelic, ]
    }
  }

  # Configure ID
  vcf_samples <- colnames(vcf)
  if (is.null(id)) {
    id <- vcf_samples[1]
  } else {
    if (!(id %in% vcf_samples)) {
      stop("The value for id, '", id, "', was not found among the list of,",
           " samples in the VCF: ",
           paste(vcf_samples, collapse = ", "))
    }
  }

  # Configure sample name
  vcf_name <- id
  if (!is.null(rename)) {
    vcf_name <- rename
  } else if (!is.null(sample_field)) {
    df <- VariantAnnotation::meta(VariantAnnotation::header(vcf))$SAMPLE
    if (!(id %in% rownames(df))) {
      stop("The value for id, '", id, "' was not found in the 'SAMPLE' ",
           "metadata within the VCF.")
    }
    if (!sample_field %in% colnames(df)) {
      stop("The value for header_sample_name_field, '",
           sample_field, "' was not found ",
           "in the 'SAMPLE' metadata within the VCF.")
    }
    vcf_name <- as.character(df[id, sample_field])
  }

  # Perform filtering based on FILTER column
  if (isTRUE(filter)) {
    # Remove filtered rows
    pass <- which(SummarizedExperiment::rowRanges(vcf)$FILTER == "PASS")
    if (length(pass) == 0) {
      warning("No variants passed the filter for VCF with the id/name: ",
              id, "/", vcf_name)
      return(NULL)
    }
    vcf <- vcf[pass, ]
  }

  # Find variants in that sample (alt allele having number other than 0)
  pass <- stringr::str_detect(VariantAnnotation::geno(vcf)$GT[, id], "[1-9]")
  if (sum(pass) == 0) {
    if (id != vcf_name) {
      warning("All variants matched the reference allele ",
              "for: ", vcf_name, " (", id, ")")

    } else {
      warning("All variants matched the reference allele ",
              "for id: ", id)
    }
    return(NULL)
  }
  vcf <- vcf[pass, ]

  # What does this do?
  rows <- SummarizedExperiment::rowRanges(vcf)
  dt <- cbind(data.table::as.data.table(rows)[, c("seqnames", "REF", "ALT",
                                                 "start", "end")],
              "sample" = vcf_name)

  data.table::setnames(dt, c("seqnames", "start", "end",
                             "REF", "ALT", "sample"),
                       .required_musica_headers())

  # Add extra columns if requested
  if (!is.null(extra_fields)) {
    info_field <- VariantAnnotation::info(vcf)
    temp <- setdiff(extra_fields, colnames(info_field))
    if (length(temp) > 0) {
      warning("Values in extra_fields were not found in the INFO field within",
              " the VCF: ", paste(temp, collapse = ", "))
    }
    temp <- intersect(extra_fields, colnames(info_field))
    dt <- cbind(dt, info_field[, temp])
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
#' alternate allele will be given their own rows. If \code{"exclude"} is
#' selected, then these rows will be removed. Default \code{"expand"}.
#' @param extra_fields Optionally extract additional fields from the \code{INFO}
#' section of the VCF. Default \code{NULL}.
#' @param fix_vcf_errors Attempt to automatically fix VCF file
#' formatting errors.
#' @return Returns a data.table of variants extracted from a vcf
#' @examples
#' vcf <- system.file("extdata", "public_LUAD_TCGA-97-7938.vcf",
#'   package = "musicatk")
#' variants <- extract_variants_from_vcf_file(vcf_file = vcf)
#' @export
extract_variants_from_vcf_file <- function(vcf_file, id = NULL, rename = NULL,
                                           sample_field = NULL,
                                           filename_as_id = FALSE,
                                           strip_extension = c(".vcf",
                                                               ".vcf.gz",
                                                               ".gz"),
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
        rename <- gsub(paste(paste0(strip_extension, "$"), collapse = "|"),
                       "", rename)
      }
    }
  }

  # Automatically try to fix some types of errors in VCF files
  if (is(vcf, "try-error") && fix_vcf_errors) {
    alt_input <- utils::read.table(vcf_file, stringsAsFactors = FALSE,
                            check.names = FALSE, comment.char = "", skip = 7,
                            header = TRUE)
    sample_header_name <- names(alt_input[10])
    vcf_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT")
    if (!all(vcf_columns %in% names(alt_input))) {
      stop(paste("VCF File: ", vcf_file,
                 " could not be recovered, please review.",
                 " \nAdditional information: \n", vcf[1], sep = ""))
    }
    non_letters <- rawToChar(as.raw(c(32:47, 58:64, 91, 93:96, 123:126)),
                             multiple = TRUE)
    malformed_rows <- c(which(alt_input$REF == 0), which(alt_input$ALT == 0),
                        grep("#", alt_input$`#CHROM`))
    if (length(malformed_rows) > 0) {
      alt_input <- alt_input[-malformed_rows, ]
    } else {
      stop(paste("VCF File: ", vcf_file,
                 " could not be recovered, please review.",
                 " \nAdditional information: \n", vcf[1], sep = ""))
    }
    alt_input[, "End_Position"] <- alt_input[, "POS"]
    alt_input <- alt_input[, colnames(alt_input)[c(seq_len(2), ncol(alt_input),
                                                   seq(4, (ncol(alt_input) - 1))
                                                   )]]

    # Needs to be changed so that it creates a VCF object which can be
    # passed to extract_variants_from_vcf
    #dt <- cbind(data.table::as.data.table(alt_input[, c("#CHROM", "POS",
    #                                                    "End_Position", "REF",
    #                                                    "ALT", "QUAL",
    #                                                    "FILTER")]),
    #            Tumor_Sample_Barcode = vcf_name, data.table::as.data.table(
    #              alt_input[, c("INFO", "FORMAT", sample_header_name), drop =
    #                          FALSE])) %>% dplyr::rename(
    #                            "Tumor_Seq_Allele1" = "REF",
    #                            "Tumor_Seq_Allele2" = "ALT", "Chromosome" =
    #                              "#CHROM", "Start_Position" = "POS",
    #                            "Alleles" = sample_header_name)
    #warning(paste("\nVCF File: ", vcf_file,
    #           " is malformed but could be recovered, review optional.",
    #           " \nAdditional information: \n", vcf[1],
    #           sep = ""))
  } else if (is(vcf, "try-error")) {
    stop("VCF File: ", vcf_file,
                  " is malformed and automatic error fixing is disabled.",
                  " \nAdditional information: \n", vcf[1])
  } else {
    dt <- extract_variants_from_vcf(vcf = vcf, id = id, rename = rename,
                                    sample_field = sample_field,
                                    multiallele = multiallele,
                                    filter = filter,
                                    extra_fields = extra_fields)
  }
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
#' package = "musicatk")
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
  dt <- extract_variants_from_matrix(mat = dt, chromosome_col = "Chromosome",
                                     start_col = "Start_Position",
                                     end_col = "End_Position",
                                     ref_col = "Tumor_Seq_Allele1",
                                     alt_col = "Tumor_Seq_Allele2",
                                     sample_col = "Tumor_Sample_Barcode",
                                     extra_fields = extra_fields)
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
#' package = "musicatk")
#' library(maftools)
#' maf <- read.maf(maf_file)
#' variants <- extract_variants_from_maf(maf = maf)
#' @export
extract_variants_from_matrix <- function(mat, chromosome_col = "chr",
                                         start_col = "start",
                                         end_col = "end",
                                         ref_col = "ref",
                                         alt_col = "alt",
                                         sample_col = "sample",
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
                       sample = sample_col)

  # Add extra columns if requested
  if (!is.null(extra_fields)) {
    temp <- setdiff(extra_fields, colnames(dt))
    if (length(temp) > 0) {
      warning("Values in 'extra_fields' were not found in column names in",
              " the MAF object: ", paste(temp, collapse = ", "))
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
#' @return Returns a data.table of variants from a maf
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#' package = "musicatk")
#' maf <- extract_variants_from_maf_file(maf_file = maf_file)
#' @export
extract_variants_from_maf_file <- function(maf_file, extra_fields = NULL) {
  maf <- maftools::read.maf(maf_file, verbose = FALSE)
  return(extract_variants_from_maf(maf = maf, extra_fields = extra_fields))
}

#' Creates a musica object from a variant table
#'
#' Add description
#'
#' @param x Any object that can be coerced to a data.table including a matrix
#' or data.frame.
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
#' @param extra_fields Which additional fields to extract and include in
#' the musica object. Default \code{NULL}.
#' @param verbose Whether to print status messages during error checking.
#' Default \code{TRUE}.
#' @return Returns a musica object
#' @examples
#' maf_file <- system.file("extdata", "public_TCGA.LUSC.maf",
#' package = "musicatk")
#' variants <- extract_variants_from_maf_file(maf_file)
#' g <- select_genome("38")
#' musica <- create_musica(x = variants, genome = g)
#' @export
create_musica <- function(x, genome,
                         check_ref_chromosomes = TRUE,
                         check_ref_bases = TRUE,
                         chromosome_col = "chr",
                         start_col = "start",
                         end_col = "end",
                         ref_col = "ref",
                         alt_col = "alt",
                         sample_col = "sample",
                         extra_fields = NULL,
                         verbose = TRUE) {

  used_fields <- c(.required_musica_headers(), extra_fields)
  if (canCoerce(x, "data.table")) {
    dt <- data.table::as.data.table(x)
  } else {
    stop("'x' needs to be an object which can be coerced to a data.table. ",
         "Valid classes include but are not limited to 'matrix', 'data.frame'",
         " and 'data.table'.")
  }
  if (!inherits(genome, "BSgenome")) {
    stop("'genome' needs to be a 'BSgenome' object containing the genome ",
         "reference that was used when calling the variants.")
  }

  # Check for necessary columns and change column names to stardard object
  dt <- .check_headers(dt,
                 chromosome = chromosome_col,
                 start = start_col,
                 end = end_col,
                 ref = ref_col,
                 alt = alt_col,
                 sample = sample_col,
                 update_fields = TRUE)

  # Subset to necessary columns and add variant type
  all_fields <- c(.required_musica_headers(), extra_fields)
  dt <- dt[, all_fields, with = FALSE]
  dt <- add_variant_type(dt)

  # Some non-variants are included (e.g. T>T). These will be removed
  non_variant <- which(dt$ref == dt$alt)
  if (length(non_variant) > 0) {
    warning(length(non_variant), " variants has the same reference and ",
            "alternate allele. These variants were excluded.")
    dt <- dt[-non_variant, ]
  }

  if (isTRUE(check_ref_chromosomes)) {
    # Check for genome style and attempt to convert variants to reference
    # genome if they don't match
    if (isTRUE(verbose)) {
      message("Checking that chromosomes in the 'variant' object match ",
              "chromosomes in the 'genome' object.")
    }
    dt <- .check_variant_genome(dt = dt, genome = genome)
  }

  if (isTRUE(check_ref_bases)) {
    if (isTRUE(verbose)) {
      message("Checking that the reference bases in the 'variant' object ",
              "match the reference bases in the 'genome' object.")
    }
    .check_variant_ref_in_genome(dt = dt, genome = genome)
  }

  # Create and return a musica object
  s <- gtools::mixedsort(unique(dt$sample))
  annot <- data.frame(Samples = factor(s, levels = s))
  dt$sample <- factor(dt$sample, levels = s)
  
  musica <- new("musica", variants = dt, sample_annotations = annot)
  return(musica)
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
    g_error <- try(genome_style <-
                     GenomeInfoDb::seqlevelsStyle(genome_u)[1],
                   silent = TRUE)
    v_error <- try(variant_style <-
                     GenomeInfoDb::seqlevelsStyle(as.character(chr_u))[1],
                   silent = TRUE)

    inter <- intersect(chr_u, genome_u)
    if (length(inter) == 0) {
      # Error 1: No matching of genome and variants
      if (is(g_error, "try-error") & is(v_error, "try-error")) {
        stop("The style of the genome references in the 'variant' and ",
             "'genome' objects did not match each other or any style ",
             "from the 'GenomeInfoDb' package. Please ensure that the entries ",
             "in the'", chr_header, "' column in the variant table match ",
             "entries in 'seqnames(genome)'. First five chromosomes:\n",
             "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
             "genome: ", paste(head(genome_u, 5), collapse = ", "))
      } else if (is(g_error, "try-error")) {
        stop("The style of the genome references in the 'variant' and ",
             "'genome' objects did not match each other. The style for the ",
             "variant table was determined to be ", variant_style, " by the ",
             "'GenomeInfoDb' package. However, the 'genome' object did not ",
             "match any known styles and therefore could not be automatically ",
             "mapped. Please ensure that the entries ",
             "in the'", chr_header, "' column in the variant table match ",
             "entries in 'seqnames(genome)'. First five chromosomes:\n",
             "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
             "genome: ", paste(head(genome_u, 5), collapse = ", "))
      } else if (is(v_error, "try-error")) {
        stop("The style of the genome references in the 'variant' and ",
             "'genome' objects did not match each other. The style for the ",
             "genome object was determined to be ", genome_style, " by the ",
             "'GenomeInfoDb' package. However, the variant table did not ",
             "match any known styles and therefore could not be automatically ",
             "mapped. Please ensure that the entries ",
             "in the'", chr_header, "' column in the variant table match ",
             "entries in 'seqnames(genome)'. First five chromosomes:\n",
             "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
             "genome: ", paste(head(genome_u, 5), collapse = ", "))
      } else {
        # Attempt to map variants to genome object
        new_chr <- chr
        map_error <- try(GenomeInfoDb::seqlevelsStyle(new_chr) <-
                           genome_style, silent = TRUE)
        if (is(map_error, "try-error")) {
          stop("The style of the genome references in the 'variant' and ",
               "'genome' objects did not match each other. The style for the ",
               "genome object was determined to be ", genome_style, " and the ",
               "style of the variants was determined to be ", variant_style,
               " by the 'GenomeInfoDb' package. However, they were not able ",
               "to be automatically mapped. Please ensure that the entries ",
               "in the'", chr_header, "' column in the variant table match ",
               "entries in 'seqnames(genome)'. First five chromosomes:\n",
               "variant: ", paste(head(chr_u, 5), collapse = ", "), "\n",
               "genome: ", paste(head(genome_u, 5), collapse = ", "))
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
            warning("The style of the genome references in the 'variant' and ",
              "'genome' objects did not match each other. The style for the ",
              "genome object was determined to be ", genome_style, " and the ",
              "style of the variants was determined to be ", variant_style,
              " by the 'GenomeInfoDb' package. ", inter_sum, " variants were ",
              "automatically converted to those in the genome object. ",
              diff.sum, " variants were not able to be converted and ",
              "were excluded from the final variant table. Variant ",
              "chromosomes that were not able to be converted: \n",
              paste(new_diff, collapse = ", "))
          } else {
            warning("The style of the genome references in the 'variant' and ",
              "'genome' objects did not match each other. The style for the ",
              "genome object was determined to be ", genome_style, " and the ",
              "style of the variants was determined to be ", variant_style,
              " by the 'GenomeInfoDb' package. All variant chromosomes were ",
              "automatically converted to those in the genome object.")
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

  genome_ref <- BSgenome::getSeq(x = genome,
                   names = chr, start = start, end = end,
                   as.character = TRUE)

  # Only check references that have bases. For example maf files will have "-"
  # for all insertions which cannot be checked.
  ix <- grepl("[^ACGT]", ref)
  no_match_sum <- sum(ref[!ix] != genome_ref[!ix])
  if (no_match_sum > 0) {
    warning("Reference bases for ", no_match_sum, " out of ",
         length(ref), " variants did not match the reference base in the ",
         "'genome'. Make sure the genome reference is correct.")
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
  if (!is.null(alt))  rfields["alt"] <- alt
  if (!is.null(sample))  rfields["sample"] <- sample

  # Check for missing columns
  col.ix <- match(tolower(rfields), tolower(colnames(dt)))
  if (sum(is.na(col.ix) > 0)) {
    missing.cols <- rfields[is.na(col.ix)]
    stop("Some required columns are missing in the maf: ",
         paste(missing.cols, collapse = ", "))
  }

  # Adjust column case if needed
  prev_col <- colnames(dt)[col.ix]
  colnames(dt)[col.ix] <- rfields
  mismatch_ix <- colnames(dt)[col.ix] != prev_col
  if (sum(mismatch_ix) > 0) {
    warning("Some columns in maf had the wrong case and were automatically ",
            "adjusted: ", paste(prev_col[mismatch_ix], collapse = ", "))
  }

  # Change columns of data.table to match the required musica format
  if (isTRUE(update_fields)) {
    data.table::setnames(dt, rfields, .required_musica_headers())
  }
  return(dt)
}

.required_musica_headers <- function() {
  return(c(chromosome = "chr", start = "start", end = "end", ref = "ref",
           alt = "alt", sample = "sample"))
}

.required_maf_headers <- function() {
  return(c(chromosome = "Chromosome", start = "Start_Position",
           end = "End_Position", ref = "Tumor_Seq_Allele1",
           alt = "Tumor_Seq_Allele2", sample = "Tumor_Sample_Barcode"))
}
