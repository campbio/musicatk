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
setClass("musica", slots = c(variants = "data.table",
                             count_tables = "list",
                             sample_annotations = "data.frame",
                             result_list = "SimpleList"),
         prototype = list(variants = data.table::data.table(),
                          count_tables = list(),
                          sample_annotations = data.frame(),
                          result_list = S4Vectors::SimpleList()))

# setMethod("show", "musica_variants",
#           function(object)cat(cat("musica object containing \n**Variants: \n"),
#                               if (!all(is.na(object@variants))) {
#                                 cat(methods::show(object@variants))
#                                 }else{
#                                   cat("Empty")
#                                     },
#                               cat("\n**Count_Tables Object containing: \n"),
#                               if (length(object@count_tables@table_name) > 0) {
#                                 cat("\n**Count Tables: \n",
#                                     apply(cbind(do.call("rbind", lapply(
#                                       object@count_tables@table_list, dim)),
#                                       "\n"), 1, paste),
#                                     "\n**Names: \n", paste(
#                                         unlist(object@count_tables@table_name),
#                                         "\n", sep = ""), "\n**Descriptions: \n",
#                                     paste(unlist(
#                                       object@count_tables@description), "\n",
#                                       sep = ""))
#                                 }else{
#                                   cat("Empty")
#                                   },
#                               cat("\n**Sample Level Annotations: \n"),
#                               if (!all(is.na(object@sample_annotations))) {
#                                 cat(methods::show(object@sample_annotations))
#                               }else{
#                                 cat("Empty")
#                               })
# )

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
  
  #Subset variants
  variants(musica) <- variants(musica)[
    which(variants(musica)$sample %in% min_samples), ]
  
  #Subset sample annotations
  if (nrow(samp_annot(musica)) != 0) {
    .overwrite_samp_annot(musica = musica, 
                          new_annot = 
                            samp_annot(musica)[which(samp_annot(musica)$Samples 
                                                     %in% min_samples), , 
                                               drop = FALSE])
    #samp_annot(musica) <- samp_annot(musica)[which(
    #  samp_annot(musica)$Samples %in% min_samples), ]
  }
  return(musica)
}

#' Creates a new musica object subsetted to only one value of a sample annotation
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
#' package = "musicatk"), sep = "\t", header=TRUE)
#'
#' samp_annot(musica_sbs96, "Tumor_Subtypes") <- annot$Tumor_Subtypes
#'
#' musica_sbs96 <- subset_musica_by_annotation(musica_sbs96, "Tumor_Subtypes", 
#' "Lung")
#' @export
subset_musica_by_annotation <- function(musica, annot_col, annot_names) {
  if (!all(annot_col %in% colnames(samp_annot(musica)))) {
    stop(paste(annot_col, " not found in annotation columns, please review.",
               sep = ""))
  }
  annotation_index <- which(samp_annot(musica)[[which(colnames(
    samp_annot(musica)) %in% annot_col)]] %in% annot_names)
  if (length(annotation_index) == 0) {
    stop(paste(annot_names, " not present in ", annot_col,
               " column, please review.", sep = "", collapse = TRUE))
  }
  .overwrite_samp_annot(musica, samp_annot(musica)[annotation_index, ])
  annotation_samples <- samp_annot(musica)$"Samples"
  tables(musica) <- .subset_count_tables(musica, samples = annotation_samples)
  variants(musica) <- variants(musica)[
    which(variants(musica)$sample %in% annotation_samples), ]
  return(musica)
}

.overwrite_samp_annot <- function(musica, new_annot) {
  eval.parent(substitute(musica@sample_annotations <- new_annot))
}

drop_na_variants <- function(variants, annot_col) {
  if (!annot_col %in% colnames(variants)) {
    stop(paste(annot_col, " not found in annotation columns, please review.",
               sep = ""))
  }
  if (length(which(is.na(variants[[annot_col]]))) == 0) {
    return(variants)
  } else {
    return(variants[-which(is.na(variants[[annot_col]])), ])
  }
}
