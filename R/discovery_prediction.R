#' @importFrom NMF nmf
NULL

# Return 'model_id', or '<model_id>.<n>' with the smallest n not already
# used by a model in 'result_name' and 'modality' of 'musica'
.unique_model_id <- function(musica, result_name, modality, model_id) {
  existing <- NULL
  if (!is.null(result_list(musica)[[result_name]])) {
    existing <- names(get_modality(musica, result_name, modality))
  }
  if (!model_id %in% existing) {
    return(model_id)
  }
  tag <- 1
  while (paste0(model_id, ".", tag) %in% existing) {
    tag <- tag + 1
  }
  new_id <- paste0(model_id, ".", tag)
  message(
    "model_id ", model_id, " already exists. model_id updated to ", new_id
  )
  return(new_id)
}

# Names of the samples (columns) of 'counts_table' with at least one
# mutation; unless 'quiet', a message names the others, which are left out
# of the fit and its exposures
.samples_with_mutations <- function(counts_table, quiet = FALSE) {
  has_mutations <- colSums(counts_table) > 0
  if (!all(has_mutations) && !quiet) {
    message(
      "These samples have no mutations in the count table and are left ",
      "out: ", paste(colnames(counts_table)[!has_mutations], collapse = ", ")
    )
  }
  return(colnames(counts_table)[has_mutations])
}

# Distinct values of the sample annotation column 'annotation', as
# character. Stops when the column doesn't exist or has no values; a
# message names samples with a missing value, which are left out.
.annotation_groups <- function(musica, annotation) {
  available <- setdiff(colnames(samp_annot(musica)), "Samples")
  if (length(annotation) != 1 || !annotation %in% available) {
    stop(
      "Sample annotation '", paste(annotation, collapse = ", "),
      "' not found. Available annotations: ",
      paste(available, collapse = ", ")
    )
  }
  groups <- as.character(samp_annot(musica)[[annotation]])
  if (anyNA(groups)) {
    message(
      "These samples have no '", annotation, "' value and are left out: ",
      paste(sample_names(musica)[is.na(groups)], collapse = ", ")
    )
  }
  group_names <- unique(groups[!is.na(groups)])
  if (length(group_names) == 0) {
    stop("Sample annotation '", annotation, "' has no values.")
  }
  return(group_names)
}

# 'musica' without its stored results, limited to the samples whose
# 'annotation' value is 'group', or all samples when 'group' is NA. Models
# fit on it are stored under ids chosen against the full object.
.musica_for_group <- function(musica, annotation, group) {
  result_list(musica) <- S4Vectors::SimpleList()
  if (is.na(group)) {
    return(musica)
  }
  return(subset_musica_by_annotation(musica, annotation, group))
}

# Reorder the rows of 'signature' to the motifs (rows) of 'counts_table';
# stop when the two don't hold the same motifs
.match_motifs <- function(signature, counts_table, modality) {
  signature_motifs <- rownames(signature)
  table_motifs <- rownames(counts_table)
  if (is.null(signature_motifs) || is.null(table_motifs)) {
    stop(
      "The signatures and the '", modality, "' count table both need motif ",
      "names as row names."
    )
  }
  if (anyDuplicated(signature_motifs)) {
    stop(
      "The signatures have duplicated motif names, for example ",
      signature_motifs[anyDuplicated(signature_motifs)], "."
    )
  }
  if (!setequal(signature_motifs, table_motifs)) {
    differ <- union(
      setdiff(signature_motifs, table_motifs),
      setdiff(table_motifs, signature_motifs)
    )
    stop(
      "The motifs of 'signature_res' don't match the '", modality,
      "' count table (", length(signature_motifs), " and ",
      length(table_motifs), " motifs; differing motifs include ",
      paste(utils::head(differ, 3), collapse = ", "),
      "). Use signatures of the same type as the count table."
    )
  }
  return(signature[table_motifs, , drop = FALSE])
}

# Exposures of the signatures in 'signature' in each sample of
# 'counts_table', scaled to the sample's mutation count, with the LDA start
# seeded by 'seed' (NULL uses the caller's RNG). Samples with no
# mutations are left out. Returns the signatures, with rows in the table's
# motif order, and the exposures.
.predict_exposures <- function(counts_table, signature, algorithm, modality,
                               seed = 1, verbose = FALSE, quiet = FALSE) {
  algorithm <- match.arg(algorithm, c("lda", "decompTumor2Sig"))
  signature <- .match_motifs(signature, counts_table, modality)
  fit_samples <- .samples_with_mutations(counts_table, quiet = quiet)
  if (length(fit_samples) == 0) {
    stop("No sample in the '", modality, "' count table has mutations.")
  }
  counts_table <- counts_table[, fit_samples, drop = FALSE]

  if (algorithm == "lda") {
    run_lda <- function() {
      lda_posterior(
        counts_table = counts_table, signature = signature,
        max.iter = 100, verbose = verbose
      )
    }
    # Seed the random start, leaving the caller's RNG state unchanged
    if (is.null(seed)) {
      lda_res <- run_lda()
    } else {
      lda_res <- withr::with_seed(seed, run_lda())
    }
    exposures <- t(lda_res$samp_sig_prob_mat)
  } else {
    if (nrow(counts_table) != 96 ||
          !all(grepl("^[ACGT]>[ACGT]_[ACGT]{3}$", rownames(counts_table)))) {
      stop(
        "algorithm = \"decompTumor2Sig\" only supports SBS96 count tables; ",
        "use algorithm = \"lda\" for '", modality, "'."
      )
    }
    decomp_res <- predict_decompTumor2Sig(counts_table, signature)
    exposures <- t(do.call(rbind, decomp_res))
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
  }

  # Multiply proportions by sample counts
  exposures <- sweep(exposures, 2, colSums(counts_table), FUN = "*")
  return(list(signatures = signature, exposures = exposures))
}

#' @title Discover mutational signatures
#' @description Mutational signatures and exposures will be discovered using
#' methods such as Latent Dirichlet Allocation (lda) or Non-Negative
#' Matrix Factorization (nmf). These algorithms will deconvolute a matrix of
#' counts for mutation types in each sample to two matrices: 1) a "signature"
#' matrix containing the probability of each mutation type in each sample and
#' 2) an "exposure" matrix containing the estimated counts for each signature
#' in each sample. Before mutational discovery can be performed, samples first
#' need to be stored in a \code{\linkS4class{musica}} object using the
#' \link{create_musica_from_variants} or \link{create_musica_from_counts}
#' function and mutation count tables need to be created using functions such as
#' \link{build_standard_table} if \link{create_musica_from_counts} was not used.
#'
#' Samples with no mutations in the count table are left out of the fit and
#' of the exposures, with a message naming them. At least two samples need
#' mutations. If \code{model_id} is already used, \code{.1}, \code{.2}, and
#' so on is added to it.
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality to use for signature discovery. Needs
#' to be the same name supplied to the table building functions such as
#' \link{build_standard_table}.
#' @param num_signatures Number of signatures to discover.
#' @param algorithm Method to use for mutational signature discovery. One of
#' \code{"lda"} or \code{"nmf"}. Default \code{"lda"}.
#' @param result_name Name for result_list entry to save the results to. Default
#' \code{"result"}.
#' @param model_id Identifier for the result. If \code{NULL}, will be
#' automatically set to the algorithm and number of signatures. Default
#' \code{NULL}.
#' @param seed Seed to be used for the random number generators in the
#' signature discovery algorithms. Default \code{1}.
#' @param nstart Number of independent random starts used in the mutational
#' signature algorithms. Default \code{10}.
#' @param par_cores Number of parallel cores to use. Only used if
#' \code{method = "nmf"}. Default \code{1}.
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica)
#' g <- select_genome("hg38")
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#' discover_signatures(
#'   musica = musica, modality = "SBS96",
#'   num_signatures = 3, algorithm = "lda", seed = 12345, nstart = 1
#' )
#' @export
discover_signatures <- function(musica, modality, num_signatures,
                                algorithm = "lda", result_name = "result",
                                model_id = NULL, seed = 1, nstart = 10,
                                par_cores = 1, make_copy = FALSE,
                                table_name = NULL) {
  if (!is.null(table_name)) {
    warning("table_name argument deprecated. Use modality instead.")
    modality <- table_name
  }

  # update global variable
  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }

  if (!methods::is(musica, "musica")) {
    stop("Input to discover_signatures must be a 'musica' object.")
  }

  algorithm <- match.arg(algorithm, c("lda", "nmf"))
  counts_table <- .extract_count_table(musica, modality)
  fit_samples <- .samples_with_mutations(counts_table)
  if (length(fit_samples) < 2) {
    stop(
      "The 'musica' object inputted must contain at least two samples ",
      "with mutations."
    )
  }
  counts_table <- counts_table[, fit_samples, drop = FALSE]

  # if no model_id provided, generate one
  if (is.null(model_id)) {
    model_id <- paste(algorithm, num_signatures, sep = "")
  }
  model_id <- .unique_model_id(musica, result_name, modality, model_id)

  if (algorithm == "lda") {
    lda_counts_table <- t(counts_table)
    if (is.null(seed)) {
      control <- list(nstart = nstart)
    } else {
      control <- list(seed = (seq_len(nstart) - 1) + seed, nstart = nstart)
    }

    lda_out <- topicmodels::LDA(lda_counts_table, num_signatures,
      control = control
    )

    lda_sigs <- exp(t(lda_out@beta))
    rownames(lda_sigs) <- colnames(lda_counts_table)
    colnames(lda_sigs) <- paste0("Signature", seq_len(num_signatures))

    weights <- t(lda_out@gamma)
    rownames(weights) <- paste0("Signature", seq_len(num_signatures))
    colnames(weights) <- rownames(lda_counts_table)

    result <- methods::new("result_model",
      signatures = lda_sigs,
      exposures = weights, num_signatures = num_signatures,
      other_parameters = S4Vectors::SimpleList(),
      credible_intervals = S4Vectors::SimpleList(),
      metrics = S4Vectors::SimpleList(), umap = matrix(),
      model_id = model_id, modality = modality
    )
    exposures(result) <- sweep(exposures(result), 2, colSums(exposures(result)),
      FUN = "/"
    )
  } else if (algorithm == "nmf") {
    # Needed to prevent error with entirely zero rows
    epsilon <- 0.00000001

    decomp <- NMF::nmf(counts_table + epsilon, num_signatures,
      seed = seed,
      nrun = nstart, .options = paste("p", par_cores,
        sep = ""
      )
    )

    rownames(decomp@fit@H) <- paste("Signature", seq_len(num_signatures),
      sep = ""
    )
    colnames(decomp@fit@W) <- paste("Signature", seq_len(num_signatures),
      sep = ""
    )
    result <- methods::new("result_model",
      signatures = decomp@fit@W,
      exposures = decomp@fit@H, num_signatures = num_signatures,
      model_id = model_id, modality = modality
    )
    signatures(result) <- sweep(signatures(result), 2,
      colSums(signatures(result)),
      FUN = "/"
    )
    exposures(result) <- sweep(exposures(result), 2, colSums(exposures(result)),
      FUN = "/"
    )
  } else {
    stop(
      "That method is not supported. Please select 'lda' or 'nmf' ",
      "to generate signatures."
    )
  }
  # Multiply Weights by sample counts
  exposures(result) <- sweep(exposures(result), 2, colSums(counts_table),
    FUN = "*"
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
  musica@result_list[[result_name]]@modality[[modality]][[model_id]] <- result

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}

#' @title Prediction of exposures in new samples using pre-existing signatures
#' @description Exposures for samples will be predicted using an existing set
#' of signatures stored in a \code{\linkS4class{result_model}} object.
#' Algorithms available for prediction include a modify version of \code{"lda"},
#' and \code{"decompTumor2Sig"}.
#'
#' Signatures are matched to the count table by motif name (row name), so
#' their rows may be in any order; prediction stops when the signatures and
#' the count table don't have the same motifs. \code{"decompTumor2Sig"} only
#' supports SBS96 tables. Samples with no mutations are left out of the
#' exposures, with a message naming them.
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality for posterior prediction. Must match the table type
#' used to generate the prediction signatures
#' @param signature_res Signatures used to predict exposures for the samples
#' \code{musica} object. Existing signatures need to stored in a
#' \code{\linkS4class{result_model}} object.
#' @param algorithm Algorithm to use for prediction of exposures. One of
#' \code{"lda"} or \code{"decompTumor2Sig"}.
#' @param result_name Name for result_list entry to save the results to. Default
#' \code{"result"}.
#' @param model_id Identifier for the result. If \code{NULL}, will be
#' automatically set to the algorithm and number of signatures. Default
#' \code{NULL}.
#' @param signatures_to_use Which signatures in the \code{signature_res} result
#' object to use. Default is to use all signatures.
#' @param verbose If \code{TRUE}, progress will be printing. Only used if
#' \code{algorithm = "lda"}. Default \code{FALSE}.
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @param seed Seed for the random start of \code{algorithm = "lda"}, so
#' results are reproducible; the caller's random number state is left
#' unchanged. If \code{NULL}, the current random number state is used.
#' Default \code{1}.
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica)
#' data(cosmic_v2_sigs)
#' g <- select_genome("hg38")
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#' result <- predict_exposure(
#'   musica = musica, modality = "SBS96",
#'   signature_res = cosmic_v2_sigs, algorithm = "lda"
#' )
#'
#' # Predict using LDA-like algorithm with seed set to 2
#' predict_exposure(
#'   musica = musica, modality = "SBS96",
#'   signature_res = cosmic_v2_sigs, algorithm = "lda", seed = 2
#' )
#' @export
predict_exposure <- function(musica, modality, signature_res,
                             algorithm = c("lda", "decompTumor2Sig"),
                             result_name = "result", model_id = NULL,
                             signatures_to_use = seq_len(ncol(
                               signatures(signature_res)
                             )), verbose = FALSE,
                             make_copy = FALSE, seed = 1, table_name = NULL) {
  if (!is.null(table_name)) {
    warning("table_name argument deprecated. Use modality instead.")
    modality <- table_name
  }

  # update global variable
  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }

  algorithm <- match.arg(algorithm)
  signature <- signatures(signature_res)[, signatures_to_use, drop = FALSE]
  counts_table <- .extract_count_table(musica, modality)
  predicted <- .predict_exposures(counts_table, signature, algorithm,
    modality,
    seed = seed, verbose = verbose
  )

  # if not model_id provided, fill one in
  if (is.null(model_id)) {
    model_id <-
      paste(algorithm, length(signatures_to_use), "_exp_pred", sep = "")
  }
  model_id <- .unique_model_id(musica, result_name, modality, model_id)

  result <- methods::new("result_model",
    signatures = predicted$signatures,
    exposures = predicted$exposures,
    num_signatures = length(signatures_to_use),
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

  # if modality does not exist, create entry
  if (is.null(get_modality(musica, result_name, modality))) {
    musica@result_list[[result_name]]@modality[[modality]] <- list()
  }

  # add result_model object in the list for the proper modality
  musica@result_list[[result_name]]@modality[[modality]][[model_id]] <- result

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}

lda_posterior <- function(counts_table, signature, max.iter = 100,
                          theta = 0.1, verbose) {
  k <- ncol(signature) # number of signatures/topics
  num_samples <- ncol(counts_table) # number of samples

  if (length(theta) == 1) {
    theta <- rep(theta, k) # symmetric singular value converted to vector
  }
  sample_count_sums <- colSums(counts_table)

  # Initialize signature proportion matrix
  samp_sig_prob_mat <- matrix(NA, nrow = num_samples, ncol = k)
  sig_mut_counts <- matrix(NA, nrow = num_samples, ncol = k)
  rownames(samp_sig_prob_mat) <-
    rownames(sig_mut_counts) <- colnames(counts_table)
  colnames(samp_sig_prob_mat) <-
    colnames(sig_mut_counts) <- colnames(signature)

  for (s in seq_len(num_samples)) {
    sig_mut_counts[s, ] <- base::tabulate(sample(
      x = seq_len(k), size =
        sample_count_sums[s],
      replace = TRUE
    ), k)
  }

  # Update signature proportion matrix
  if (verbose) {
    message("Calculating Signature Proportions")
  }
  for (i in seq_len(max.iter)) {
    for (s in seq_len(num_samples)) {
      # updating each mutation probability to reassign to a signature
      log_prob_mut_reassignment <-
        digamma(sig_mut_counts[s, ] + theta) -
        digamma(sample_count_sums[s] + sum(theta))
      # updating present sample topic probability
      sig_sample_weights <- t(signature + 1e-20) *
        exp(log_prob_mut_reassignment) # avoid 0 in norm
      sig_sample_weights <- sweep(sig_sample_weights,
        MARGIN = 2, STATS =
          colSums(sig_sample_weights), FUN = "/"
      )
      # assigned counts for a topic for a sample
      updated_topic_motifs <- counts_table[, s] * t(sig_sample_weights)

      # Update nN.SbyT[s, ] sample counts assigned to signature
      sig_mut_counts[s, ] <- colSums(updated_topic_motifs)

      # Update p.SbyT[s, ]
      samp_sig_prob_mat[s, ] <- (sig_mut_counts[s, ]) / (sample_count_sums[s])
    }
    # Update theta; a Dirichlet can't be fit to one sample, so keep the prior
    if (num_samples > 1) {
      theta <- MCMCprecision::fit_dirichlet(x = samp_sig_prob_mat +
                                              .Machine$double.eps)$alpha
    }

    if (verbose) {
      message(theta)
    }
  }
  return(list(samp_sig_prob_mat = samp_sig_prob_mat, theta.poster = theta))
}

predict_decompTumor2Sig <- function(sample_mat, signature_mat) {
  # Alexandrov-type prediction
  input_signatures_normalized <- apply(
    signature_mat, 2,
    function(x) {
      x / sum(x)
    }
  )
  signatures <- split(
    input_signatures_normalized,
    col(input_signatures_normalized)
  )
  # signatures_ref <- decompTumor2Sig::readAlexandrovSignatures()
  ns <- as.matrix(row.names(signature_mat))
  ns <- apply(ns, 1, function(x) {
    stringr::str_c(
      substr(x, 5, 5), "[", substr(x, 1, 3), "]",
      substr(x, 7, 7)
    )
  })
  signatures <- lapply(signatures, setNames, ns)

  input_samples_normalized <- apply(sample_mat, 2, function(x) {
    x / sum(x)
  })
  input_samples1 <- split(
    input_samples_normalized,
    col(input_samples_normalized)
  )
  genomes <- lapply(input_samples1, setNames, ns)

  sample_weight_mat <- decompTumor2Sig::decomposeTumorGenomes(genomes,
    signatures,
    verbose = FALSE
  )
  return(sample_weight_mat)
}

# placeholder
.multi_modal_discovery <- function(musica, num_signatures, motif96_name,
                                   rflank_name, lflank_name, max.iter = 125) {
  motif96 <- .extract_count_table(musica, motif96_name)
  rflank <- .extract_count_table(musica, rflank_name)
  lflank <- .extract_count_table(musica, lflank_name)
  message(dim(motif96))
  message(dim(rflank))
  message(dim(lflank))
}

#' Generate result_grid from musica based on annotation and range of k
#'
#' Discovers signatures for each number of signatures from \code{k_start} to
#' \code{k_end}, either in all samples or separately in each group of samples
#' defined by an annotation column.
#'
#' @details
#' Models are stored in \code{result_name} as \code{<algorithm><k>}, for
#' example \code{lda3}, or with an annotation as
#' \code{<algorithm><k>_<group>}, for example \code{lda3_Lung}. If a name is
#' already used, for example by an earlier run, \code{.1}, \code{.2}, and so
#' on is added to it, and earlier models are kept.
#'
#' Each model's \code{metrics()} holds a \code{reconstruction_error}: for
#' each sample, the mean over motifs of the squared difference between the
#' observed and reconstructed counts, then the mean of the squares of those
#' per-sample values. It is on the scale of squared counts, so it is not
#' comparable with the reconstruction error of \link{compare_k_vals}.
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality used for signature discovery
#' @param algorithm Algorithm for signature discovery
#' @param annotation Name of a sample annotation column. Signatures are
#' discovered separately in the samples of each of its values; samples with
#' a missing value are left out. If \code{NA}, all samples are used together.
#' Default \code{NA}.
#' @param k_start Lower range of number of signatures for discovery
#' @param k_end Upper range of number of signatures for discovery
#' @param result_name Name for result_list entry to save the results to. Default
#' \code{"result_grid"}.
#' @param n_start Number of times to discover signatures and compare based on
#' posterior loglikihood
#' @param seed Seed to use for reproducible results, set to null to disable
#' @param par_cores Number of parallel cores to use (NMF only)
#' @param verbose Whether to output loop iterations
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica_sbs96)
#' grid <- generate_result_grid(musica_sbs96, "SBS96", "lda",
#'   k_start = 2,
#'   k_end = 5
#' )
#' @export
generate_result_grid <- function(musica, modality, algorithm = "lda",
                                 annotation = NA, k_start, k_end,
                                 result_name = "result_grid",
                                 n_start = 1,
                                 seed = NULL, par_cores = FALSE,
                                 verbose = FALSE, make_copy = FALSE,
                                 table_name = NULL) {
  if (!is.null(table_name)) {
    warning("table_name argument deprecated. Use modality instead.")
    modality <- table_name
  }

  # update global variable
  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }

  # if this result name does not exist in result_list
  if (is.null(result_list(musica)[[result_name]])) {
    # make new list entry with the desired name and assign a new
    # result_collection object
    musica@result_list[[result_name]] <- new("result_collection",
      modality = SimpleList(),
      parameter = list(), hyperparameter = list()
    )
  }

  # if modality does not exist, create entry
  if (is.null(get_modality(musica, result_name, modality))) {
    musica@result_list[[result_name]]@modality[[modality]] <- list()
  }

  # One group of samples per value of the annotation, or all samples
  if (!is.na(annotation)) {
    group_names <- .annotation_groups(musica, annotation)
  } else {
    group_names <- NA
  }

  for (group in group_names) {
    if (is.na(group)) {
      id_suffix <- ""
    } else {
      if (verbose) {
        message("Current Annotation: ", group)
      }
      id_suffix <- paste0("_", group)
    }
    cur_musica <- .musica_for_group(musica, annotation, group)
    # Used for reconstruction error
    cur_counts <- .extract_count_table(cur_musica, modality)

    # Define new results
    for (cur_k in k_start:k_end) {
      cur_model_name <- .unique_model_id(
        musica, result_name, modality,
        paste0(algorithm, cur_k, id_suffix)
      )
      cur_musica <- tryCatch(
        discover_signatures(
          musica = cur_musica, modality = modality,
          num_signatures = cur_k, algorithm = algorithm,
          result_name = result_name, model_id = cur_model_name,
          nstart = n_start, seed = seed, par_cores = par_cores,
          make_copy = TRUE
        ),
        error = function(e) {
          if (is.na(group)) {
            stop(e)
          }
          stop("In group '", group, "': ", conditionMessage(e), call. = FALSE)
        }
      )

      cur_result <- get_model(cur_musica, result_name, modality, cur_model_name)

      # Samples with no mutations are not in the model
      fit_samples <- colnames(exposures(cur_result))
      recon_error <- mean(vapply(fit_samples, function(x) {
        mean((cur_counts[, x, drop = FALSE] -
                reconstruct_sample(cur_result, x))^2)
      }, FUN.VALUE = 0)^2)

      # add result_model object in the list for the proper modality
      musica@result_list[[result_name]]@
        modality[[modality]][[cur_model_name]] <- cur_result

      metrics(musica, result_name, modality, cur_model_name) <-
        SimpleList("reconstruction_error" = recon_error)
    }
  }

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}

reconstruct_sample <- function(result, sample_number) {
  reconstruction <- matrix(apply(
    sweep(signatures(result), 2,
      exposures(result)[, sample_number,
        drop = FALSE
      ],
      FUN = "*"
    ),
    1, sum
  ), dimnames = list(
    rownames(signatures(result)),
    "Reconstructed"
  ))
  return(reconstruction)
}

#' Automatic filtering of signatures for exposure prediction gridded across
#' specific annotation
#'
#' @param musica Input samples to predict signature weights
#' @param modality Modality used for posterior prediction (e.g. SBS96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", and decompTumor2Sig
#' @param model_id Name of model
#' @param result_name Name for result_list entry to save the results to. Default
#' \code{"result"}.
#' @param sample_annotation Annotation to grid across, if none given,
#' prediction subsetting on all samples together
#' @param min_exists Threshold to consider a signature active in a sample
#' @param proportion_samples Threshold of samples to consider a signature
#' active in the cohort
#' @param rare_exposure A sample will be considered active in the cohort if at
#' least one sample has more than this threshold proportion
#' @param verbose Print current annotation value being predicted on
#' @param combine_res Automatically combines a list of annotation results
#' into a single result object with zero exposure values for signatures not
#' found in a given annotation's set of samples
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @param seed Seed for the random start of \code{algorithm = "lda"}; see
#' \link{predict_exposure}. Default \code{1}.
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica_annot)
#' data(cosmic_v2_sigs)
#' auto_predict_grid(
#'   musica = musica_annot, modality = "SBS96",
#'   signature_res = cosmic_v2_sigs, algorithm = "lda",
#'   sample_annotation = "Tumor_Subtypes"
#' )
#' auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda")
#' @export
auto_predict_grid <- function(musica, modality, signature_res, algorithm,
                              model_id = NULL, result_name = "result",
                              sample_annotation = NULL, min_exists = 0.05,
                              proportion_samples = 0.25, rare_exposure = 0.4,
                              verbose = TRUE, combine_res = TRUE,
                              make_copy = FALSE, seed = 1,
                              table_name = NULL) {
  if (!is.null(table_name)) {
    warning("table_name argument deprecated. Use modality instead.")
    modality <- table_name
  }

  # update global variable
  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }

  # if this result name does not exist in result_list
  if (is.null(result_list(musica)[[result_name]])) {
    # make new list entry with the desired name and assign a new
    # result_collection object
    musica@result_list[[result_name]] <- new("result_collection",
      modality = SimpleList(),
      parameter = list(), hyperparameter = list()
    )
  }

  # if modality does not exist, create entry
  if (is.null(get_modality(musica, result_name, modality))) {
    musica@result_list[[result_name]]@modality[[modality]] <- list()
  }

  if (is.null(sample_annotation)) {
    combine_res <- FALSE
    musica <- auto_subset_sigs(
      musica = musica, modality = modality,
      signature_res = signature_res, algorithm = algorithm,
      min_exists = min_exists,
      proportion_samples = proportion_samples,
      rare_exposure = rare_exposure, result_name = result_name,
      model_id = model_id, seed = seed
    )
  } else {
    group_names <- .annotation_groups(musica, sample_annotation)
    group_ids <- character(0)
    for (group in group_names) {
      if (verbose) {
        message(group)
      }
      cur_model_name <- .unique_model_id(musica, result_name, modality, group)
      current_musica <- .musica_for_group(musica, sample_annotation, group)
      current_musica <- tryCatch(
        auto_subset_sigs(
          musica = current_musica, modality = modality,
          signature_res = signature_res, min_exists = min_exists,
          proportion_samples = proportion_samples,
          rare_exposure = rare_exposure, algorithm = algorithm,
          result_name = result_name, model_id = cur_model_name, seed = seed
        ),
        error = function(e) {
          stop("In group '", group, "': ", conditionMessage(e), call. = FALSE)
        }
      )

      cur_result <- get_model(current_musica, result_name, modality,
                              cur_model_name)

      musica@result_list[[result_name]]@
        modality[[modality]][[cur_model_name]] <- cur_result
      group_ids <- c(group_ids, cur_model_name)
    }
  }
  if (combine_res) {
    musica <- combine_predict_grid(musica, modality, signature_res, group_ids,
      result_name = result_name, make_copy = TRUE,
      model_rename = model_id
    )
  }

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}

#' Automatic filtering of inactive signatures
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality used for posterior prediction (e.g. SBS96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior" and decompTumor2Sig
#' @param min_exists Threshold to consider a signature active in a sample
#' @param proportion_samples Threshold of samples to consider a signature
#' active in the cohort
#' @param rare_exposure A sample will be considered active in the cohort if at
#' least one sample has more than this threshold proportion
#' @param result_name Name for result_list entry to save the results to. Default
#' \code{"result"}.
#' @param model_id Identifier for the result. If \code{NULL}, will be
#' automatically set to the algorithm and number of signatures. Default
#' \code{NULL}.
#' @param seed Seed for the random start of \code{algorithm = "lda"}.
#' Default \code{1}.
#' @return Returns new \code{\linkS4class{musica}} object with the results.
#' @keywords internal
auto_subset_sigs <- function(musica, modality, signature_res, algorithm,
                             min_exists = 0.05, proportion_samples = 0.25,
                             rare_exposure = 0.4, result_name = "result",
                             model_id = NULL, seed = 1) {
  # First pass with every signature; not stored. Messages about samples
  # with no mutations come from the second pass.
  first_pass <- .predict_exposures(
    .extract_count_table(musica, modality), signatures(signature_res),
    algorithm, modality,
    seed = seed, quiet = TRUE
  )

  # Proportions in each sample
  exposures <- first_pass$exposures
  num_samples <- ncol(exposures)
  exposures <- sweep(exposures, 2, colSums(exposures), "/")
  to_use <- as.numeric(which(apply(exposures, 1, function(x) {
    sum(x > min_exists) / num_samples
  }) > proportion_samples |
    apply(exposures, 1, max) > rare_exposure))
  if (length(to_use) == 0) {
    stop(
      "No signature passed the thresholds (min_exists = ", min_exists,
      ", proportion_samples = ", proportion_samples, ", rare_exposure = ",
      rare_exposure, "). Lower them to keep more signatures."
    )
  }

  final_inferred <- predict_exposure(
    musica = musica, modality = modality,
    signature_res = signature_res,
    signatures_to_use = to_use,
    algorithm = algorithm,
    result_name = result_name,
    model_id = model_id,
    make_copy = TRUE,
    seed = seed
  )

  return(final_inferred)
}

#' Combine signatures and exposures of different models. Exposure values are
#' zero for samples in an annotation where that signature was not predicted
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality used for prediction.
#' @param signature_res Signatures to automatically subset from for prediction
#' @param model_ids Vector of ids for the models to combine. If null, all models
#' in the modality and result_list entry will be combined. Default \code{NULL}.
#' @param result_name Name of the result list entry containing the signatures
#' to plot. Default \code{"result"}.
#' @param model_rename New model identifier. If null, will be combination of
#' the ids for the models being combined. Deafult \code{NULL}.
#' @param make_copy If \code{FALSE}, the inputted \code{\linkS4class{musica}}
#' object is updated and nothing is returned. If \code{TRUE}, a new
#' \code{\linkS4class{musica}} object is created and returned. Default
#' \code{FALSE}.
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica_annot)
#' data(cosmic_v2_sigs)
#' grid <- auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda",
#'   "Tumor_Subtypes",
#'   combine_res = FALSE, make_copy = TRUE
#' )
#' combined <- combine_predict_grid(grid, "SBS96", cosmic_v2_sigs,
#' make_copy = TRUE)
#' @export
combine_predict_grid <- function(musica, modality, signature_res,
                                 model_ids = NULL, result_name = "result",
                                 model_rename = NULL, make_copy = FALSE,
                                 table_name = NULL) {
  if (!is.null(table_name)) {
    warning("table_name argument deprecated. Use modality instead.")
    modality <- table_name
  }

  if (make_copy == FALSE) {
    var_name <- deparse(substitute(musica))
  }

  # check if valid result_name
  if (!(result_name %in% names(result_list(musica)))) {
    stop(
      result_name, " does not exist in the result_list. Current names are: ",
      paste(names(result_list(musica)), collapse = ", ")
    )
  }

  # check if valid modality
  if (!(modality %in%
        names(get_result_list_entry(musica, result_name)@modality))) {
    stop(
      modality, " is not a valid modality. Current modalities are: ",
      paste(names(get_result_list_entry(musica, result_name)@modality),
            collapse = ", ")
    )
  }

  if (is.null(model_ids)) {
    model_ids <- names(get_modality(musica, result_name, modality))
  }

  grid_list <- get_modality(musica, result_name, modality)[model_ids]

  sig_names <- NULL
  for (id in model_ids) {
    sig_names <-
      c(sig_names, rownames(exposures(musica, result_name, modality, id)))
  }
  sig_names <- unique(sig_names)
  sig_names <- sig_names[order(sig_names)]

  comb <- NULL
  for (i in seq_len(length(grid_list))) {
    if (!table_selected(grid_list[[i]]) %in% modality) {
      stop(
        "Result number: ", i, " was not in selected table_name: ",
        modality
      )
    }
    samp <- exposures(grid_list[[i]])
    missing <- sig_names[!sig_names %in% rownames(samp)]
    missing_mat <- matrix(0, length(missing), ncol(samp))
    rownames(missing_mat) <- missing
    samp <- rbind(samp, missing_mat)
    samp <- samp[order(rownames(samp)), , drop = FALSE]
    comb <- cbind(comb, samp)
  }
  # Each sample must come from exactly one model; then restore sample order
  shared <- unique(colnames(comb)[duplicated(colnames(comb))])
  if (length(shared) > 0) {
    stop(
      "These samples are in more than one model being combined: ",
      paste(shared, collapse = ", "), ". Use 'model_ids' to choose models ",
      "that hold different samples."
    )
  }
  samples <- sample_names(musica)
  unknown <- setdiff(colnames(comb), samples)
  if (length(unknown) > 0) {
    stop(
      "These exposure samples are not in the 'musica' object: ",
      paste(unknown, collapse = ", ")
    )
  }
  comb <- comb[, samples[samples %in% colnames(comb)], drop = FALSE]

  musica@result_list[[result_name]]@modality[[modality]] <-
    get_modality(musica, result_name, modality)[!names(
      get_modality(musica, result_name, modality)) %in% model_ids]

  if (is.null(model_rename)) {
    model_id <- paste(model_ids, collapse = ".")
  } else {
    model_id <- model_rename
  }
  model_id <- .unique_model_id(musica, result_name, modality, model_id)

  # Signatures in the same motif order as the models being combined, which
  # predict_exposure() put in the count table's order
  comb_signatures <- signatures(signature_res)[, sig_names, drop = FALSE]
  model_motifs <- rownames(signatures(grid_list[[1]]))
  if (!is.null(model_motifs) &&
        setequal(model_motifs, rownames(comb_signatures))) {
    comb_signatures <- comb_signatures[model_motifs, , drop = FALSE]
  }

  grid_res <- new("result_model",
    exposures = comb,
    signatures = comb_signatures,
    num_signatures = length(sig_names),
    model_id = model_id,
    modality = modality
  )

  musica@result_list[[result_name]]@modality[[modality]][[model_id]] <- grid_res

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}
