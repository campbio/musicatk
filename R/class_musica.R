# Primary variant object/methods -------------------------------

#' The primary object that contains variants, count_tables,
#' and samples annotations
#'
#' @slot variants \code{data.table} of variants
#' @slot count_tables Summary table with per-sample unnormalized motif counts
#' @slot sample_annotations Sample-level annotations (e.g. age, sex, primary)
#' @slot result_list Results from various algorithms, modalities, and models
#' @export
#' @exportClass musica
setClass("musica",
  slots = c(
    variants = "data.table",
    count_tables = "list",
    sample_annotations = "data.frame",
    result_list = "SimpleList"
  ),
  prototype = list(
    variants = data.table::data.table(),
    count_tables = list(),
    sample_annotations = data.frame(),
    result_list = S4Vectors::SimpleList()
  )
)

# Sample-Level object/methods -------------------------------

#' Creates a new musica subsetted to only samples with enough variants
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of table used for subsetting
#' @param num_counts Minimum sum count value to drop samples
#' @return Returns a new musica object with sample annotations, count tables,
#' and variants subsetted to only contains samples with the specified minimum
#' number of counts (column sums) in the specified table
#' @examples
#' data(musica_sbs96)
#' subset_musica_by_counts(musica_sbs96, "SBS96", 20)
#' @export
subset_musica_by_counts <- function(musica, table_name, num_counts) {
  tab <- .extract_count_table(musica, table_name)
  min_samples <- colnames(tab)[which(colSums(tab) >= num_counts)]

  tables(musica) <- .subset_count_tables(musica, min_samples)

  # Subset variants
  variants(musica) <- variants(musica)[
    which(variants(musica)$sample %in% min_samples),
  ]

  # Subset sample annotations
  if (nrow(samp_annot(musica)) != 0) {
    .overwrite_samp_annot(
      musica = musica,
      new_annot =
        samp_annot(musica)[
          which(samp_annot(musica)$Samples
                %in% min_samples), ,
          drop = FALSE
        ]
    )
    # samp_annot(musica) <- samp_annot(musica)[which(
    #  samp_annot(musica)$Samples %in% min_samples), ]
  }
  return(musica)
}

#' Creates a new musica object subsetted to only one value of a sample
#' annotation
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param annot_col Annotation class to use for subsetting
#' @param annot_names Annotational value to subset to
#' @return Returns a new musica object with sample annotations, count tables,
#' and variants subsetted to only contains samples of the specified annotation
#' type
#' @examples
#' data(musica_sbs96)
#' annot <- read.table(system.file("extdata", "sample_annotations.txt",
#'   package = "musicatk"
#' ), sep = "\t", header = TRUE)
#'
#' samp_annot(musica_sbs96, "Tumor_Subtypes") <- annot$Tumor_Subtypes
#'
#' musica_sbs96 <- subset_musica_by_annotation(
#'   musica_sbs96, "Tumor_Subtypes",
#'   "Lung"
#' )
#' @export
subset_musica_by_annotation <- function(musica, annot_col, annot_names) {
  if (!all(annot_col %in% colnames(samp_annot(musica)))) {
    stop(paste(annot_col, " not found in annotation columns, please review.",
      sep = ""
    ))
  }
  annotation_index <- which(samp_annot(musica)[[which(colnames(
    samp_annot(musica)
  ) %in% annot_col)]] %in% annot_names)
  if (length(annotation_index) == 0) {
    stop(paste(annot_names, " not present in ", annot_col,
      " column, please review.",
      sep = "", collapse = TRUE
    ))
  }
  .overwrite_samp_annot(musica, samp_annot(musica)[annotation_index, ])
  annotation_samples <- samp_annot(musica)$"Samples"
  tables(musica) <- .subset_count_tables(musica, samples = annotation_samples)
  variants(musica) <- variants(musica)[
    which(variants(musica)$sample %in% annotation_samples),
  ]
  return(musica)
}

.overwrite_samp_annot <- function(musica, new_annot) {
  eval.parent(substitute(musica@sample_annotations <- new_annot))
}

drop_na_variants <- function(variants, annot_col) {
  if (!annot_col %in% colnames(variants)) {
    stop(paste(annot_col, " not found in annotation columns, please review.",
      sep = ""
    ))
  }
  if (length(which(is.na(variants[[annot_col]]))) == 0) {
    return(variants)
  } else {
    return(variants[-which(is.na(variants[[annot_col]])), ])
  }
}


# Variant-Level object/methods -------------------------------

#' Return sample from musica_variant object
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param sample_name Sample name to subset by
#' @return Returns sample data.frame subset to a single sample
#' @examples
#' data(musica)
#' subset_variants_by_samples(musica, "TCGA-94-7557-01A-11D-2122-08")
#' @export
subset_variants_by_samples <- function(musica, sample_name) {
  return(variants(musica)[
    which(variants(musica)$sample == sample_name),
  ])
}

################################################################################

#' Rename signatures for a model
#'
#' @param musica A \code{\linkS4class{musica}} object containing a mutational
#' signature discovery or prediction.
#' @param model_id The name of the model to rename signatures for.
#' @param name_vector Vector of user-defined signature names
#' @param modality The modality of the model. Must be "SBS96", "DBS78", or
#' "IND83". Default \code{"SBS96"}.
#' @param result_name Name of the result list entry containing the model.
#' Default \code{"result"}.
#' @return Musica object with user-defined signatures names
#' @examples
#' data(res)
#' name_signatures(res,
#'   model_id = "res",
#'   name_vector = c("smoking", "apobec", "unknown")
#' )
#' @export
name_signatures <- function(musica, model_id, name_vector, modality = "SBS96",
                            result_name = "result") {
  # Get result object from musica object
  result <- get_model(musica,
    result = result_name,
    modality = modality,
    model = model_id
  )

  num_sigs <- length(colnames(signatures(result)))
  if (length(name_vector) != num_sigs) {
    stop(
      "Please provide a full list of signatures names (length = ",
      num_sigs, ")."
    )
  }
  eval.parent(substitute(colnames(signatures(musica, result_name, modality,
                                             model_id)) <- name_vector))
  eval.parent(substitute(rownames(exposures(musica, result_name, modality,
                                            model_id)) <- name_vector))
}

################################################################################

#' @title Retrieve variants from a musica object
#' @description The \code{variants} \code{data.table} contains the variants
#' and variant-level annotations
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @rdname variants
#' @return A data.table of variants
#' @export
#' @examples
#' data(res)
#' variants(res)
setGeneric(
  name = "variants",
  def = function(object) {
    standardGeneric("variants")
  }
)

#' @rdname variants
setMethod(
  f = "variants",
  signature = "musica",
  definition = function(object) {
    return(object@variants)
  }
)

#' @rdname variants
#' @param musica A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function
#' @param value A \code{\linkS4class{data.table}} of mutational variants and
#' variant-level annotations
#' @export
setGeneric(
  name = "variants<-",
  def = function(musica, value) {
    standardGeneric("variants<-")
  }
)

#' @rdname variants
setReplaceMethod(
  f = "variants",
  signature = c("musica", "data.table"),
  definition = function(musica, value) {
    musica@variants <- value
    return(musica)
  }
)

#' @title Retrieve the list of count_tables from a musica object
#' @description The \code{count_tables} contains standard and/or custom
#' count tables created from variants
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @rdname tables
#' @return A list of count_tables
#' @export
#' @examples
#' data(res)
#' tables(res)
setGeneric(
  name = "tables",
  def = function(object) {
    standardGeneric("tables")
  }
)

#' @rdname tables
setMethod(
  f = "tables",
  signature = "musica",
  definition = function(object) {
    return(object@count_tables)
  }
)

#' @rdname tables
#' @param musica A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @param value A list of \code{\linkS4class{count_table}} objects representing
#' counts of motifs in samples
#' @export
setGeneric(
  name = "tables<-",
  def = function(musica, value) {
    standardGeneric("tables<-")
  }
)

#' @rdname tables
setReplaceMethod(
  f = "tables",
  signature = c("musica", "list"),
  definition = function(musica, value) {
    if (length(value) > 0) {
      if (!all(unlist(lapply(value, is, "count_table")))) {
        stop("All objects in the list must be count tables.")
      }
    }
    musica@count_tables <- value
    return(musica)
  }
)

#' @title Get or set sample annotations from a musica object
#' @description  Sample annotations can be used to store information about
#' each sample such as tumor type or treatment status. These are used in
#' downstream plotting functions such as \code{\link{plot_exposures}} or
#' \code{\link{plot_umap}} to group or color samples by a particular annotation.
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @param name The name of the new annotation to add.
#' @param value A vector containing the new sample annotations. Needs to be
#' the same length as the number of samples in the object.
#' @rdname samp_annot
#' @return A new object with the sample annotations added to the table in the
#' \code{sample_annotations} slot.
#' @seealso See \code{\link{sample_names}} to get a vector of sample names in
#' the \code{\linkS4class{musica}} object.
#' @export
#' @examples
#' data(res_annot)
#' samp_annot(res_annot)
#'
#' # Add new annotation
#' samp_annot(res_annot, "New_Annotation") <- rep(c("A", "B"), c(3, 4))
#' samp_annot(res_annot)
setGeneric(
  name = "samp_annot",
  def = function(object) {
    standardGeneric("samp_annot")
  }
)

#' @rdname samp_annot
setMethod(
  f = "samp_annot",
  signature = "musica",
  definition = function(object) {
    return(object@sample_annotations)
  }
)

#' @rdname samp_annot
#' @examples
#' data(musica)
#' samp_annot(musica, "example") <- rep("ex", 7)
#' @export
setGeneric(
  name = "samp_annot<-",
  def = function(object, name, value) {
    standardGeneric("samp_annot<-")
  }
)

#' @rdname samp_annot
setReplaceMethod(
  f = "samp_annot",
  signature = c("musica", "character", "vector"),
  definition = function(object, name, value) {
    # Initialize sample annotations in musica object
    if (nrow(samp_annot(object)) == 0) {
      samples <- unique(object@variants$sample)
      sample_dt <- data.table::data.table(Samples = samples)
      object@sample_annotations <- sample_dt
    }

    if (length(name) != 1 || !is.character(name)) {
      stop("The 'name' parameter must be a character of length 1.")
    }
    if (length(value) != nrow(object@sample_annotations)) {
      stop(
        "The new sample annotation vector in 'value' must have the ",
        "length as the number of samples in the 'musica' object. Number ",
        "of samples in 'musica' object: ", nrow(object@sample_annotations),
        ". Length of new sample annotation: ",
        length(value), "."
      )
    }
    object@sample_annotations[[name]] <- value
    return(object)
  }
)

#' @title Retrieve sample names from a musica object
#' @description Sample names were included in the \code{sample} column
#' in the variant object passed to \code{\link{create_musica_from_variants}}, or
#' in the colnames of the count table object passed to
#' \code{\link{create_musica_from_counts}}. This returns a unique list of
#' samples names in the order they are inside the \code{\linkS4class{musica}}
#' object.
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @rdname sample_names
#' @return A character vector of sample names
#' @export
#' @examples
#' data(res)
#' sample_names(res)
setGeneric(
  name = "sample_names",
  def = function(object) {
    standardGeneric("sample_names")
  }
)

#' @rdname sample_names
setMethod(
  f = "sample_names",
  signature = "musica",
  definition = function(object) {
    if (is.null(object@sample_annotations) ||
        nrow(object@sample_annotations) == 0 ||
        is.null(object@sample_annotations$Sample)) {
      s <- gtools::mixedsort(unique(object@variants$sample))
    } else {
      s <- object@sample_annotations$Sample
    }
    return(s)
  }
)

#' @title Retrieve result_list from a musica object
#' @description The \code{result_list} contains results from various runs
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @rdname result_list
#' @return A list of results
#' @export
#' @examples
#' data(res)
#' result_list(res)
setGeneric(
  name = "result_list",
  def = function(object) {
    standardGeneric("result_list")
  }
)

#' @rdname result_list
setMethod(
  f = "result_list",
  signature = "musica",
  definition = function(object) {
    return(object@result_list)
  }
)

#' @rdname result_list
#' @param musica A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function
#' @param value A list of results
#' @export
setGeneric(
  name = "result_list<-",
  def = function(musica, value) {
    standardGeneric("result_list<-")
  }
)

#' @rdname result_list
setReplaceMethod(
  f = "result_list",
  signature = c("musica", "SimpleList"),
  definition = function(musica, value) {
    musica@result_list <- value
    return(musica)
  }
)

#' @title Retrieve result_list entry from a musica object
#' @description The \code{result_list} contains results from various runs
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @param result_name The name of the result_list entry.
#' @rdname get_result_list_entry
#' @return A list of results
#' @export
#' @examples
#' data(res)
#' get_result_list_entry(res, "result")
setGeneric(
  name = "get_result_list_entry",
  def = function(object, result_name) {
    standardGeneric("get_result_list_entry")
  }
)

#' @rdname get_result_list_entry
setMethod(
  f = "get_result_list_entry",
  signature = c("musica", "character"),
  definition = function(object, result_name) {
    return(object@result_list[[result_name]])
  }
)


################################################ 33

#' @title Retrieve the names of count_tables from a musica object
#' @description The \code{count_tables} contains standard and/or custom
#' count tables created from variants
#' @param object A \code{\linkS4class{musica}} object generated by
#' the \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function.
#' @rdname built_tables
#' @return The names of created count_tables
#' @export
#' @examples
#' data(res)
#' built_tables(res)
setGeneric(
  name = "built_tables",
  def = function(object) {
    standardGeneric("built_tables")
  }
)

#' @rdname built_tables
setMethod(
  f = "built_tables",
  signature = "musica",
  definition = function(object) {
    return(names(object@count_tables))
  }
)
