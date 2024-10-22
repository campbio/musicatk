#' @importFrom NMF nmf
NULL

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
#' g <- select_genome("19")
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
  if (dim(counts_table)[2] < 2) {
    stop("The 'musica' object inputted must contain at least two samples.")
  }
  present_samples <- which(colSums(counts_table) > 0)
  counts_table <- counts_table[, present_samples]

  # if no model_id provided, generate one
  if (is.null(model_id)) {
    model_id <- paste(algorithm, num_signatures, sep = "")
  }

  # check if model_id is unique and update if not
  if (!is.null(names(result_list(musica)))) {
    if (!is.null(get_result_list_entry(musica, result_name))) {
      if (!is.null(names(get_result_list_entry(musica, result_name)
                         @modality))) {
        if (model_id %in% names(get_modality(musica, result_name, modality))) {
          original_id <- model_id
          tag <- 1
          while (model_id %in%
                 names(get_modality(musica, result_name, modality))) {
            if (tag > 1) {
              model_id <- substr(model_id, 1, nchar(model_id) - 2)
            }
            model_id <- paste(model_id, ".", tag, sep = "")
            tag <- tag + 1
          }
          message("model_id ", original_id,
                  " already exists. model_id updated to ", model_id)
        }
      }
    }
  }

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
  sample_counts <- colSums(counts_table)
  matched <- match(colnames(counts_table), names(sample_counts))
  exposures(result) <- sweep(exposures(result), 2, sample_counts[matched],
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
#' @param table_name Use modality instead
#' @return Returns nothing or a new \code{\linkS4class{musica}} object,
#' depending on the \code{make_copy} parameter.
#' @examples
#' data(musica)
#' data(cosmic_v2_sigs)
#' g <- select_genome("19")
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#' result <- predict_exposure(
#'   musica = musica, modality = "SBS96",
#'   signature_res = cosmic_v2_sigs, algorithm = "lda"
#' )
#'
#' # Predict using LDA-like algorithm with seed set to 1
#' set.seed(1)
#' predict_exposure(
#'   musica = musica, modality = "SBS96",
#'   signature_res = cosmic_v2_sigs, algorithm = "lda"
#' )
#' @export
predict_exposure <- function(musica, modality, signature_res,
                             algorithm = c("lda", "decompTumor2Sig"),
                             result_name = "result", model_id = NULL,
                             signatures_to_use = seq_len(ncol(
                               signatures(signature_res)
                             )), verbose = FALSE,
                             make_copy = FALSE, table_name = NULL) {
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
  present_samples <- which(colSums(counts_table) > 0)
  counts_table <- counts_table[, present_samples, drop = FALSE]

  if (algorithm %in% c("lda_posterior", "lda", "lda_post")) {
    lda_res <- lda_posterior(
      counts_table = counts_table, signature =
        signature, max.iter = 100, verbose = verbose
    )
    exposures <- t(lda_res$samp_sig_prob_mat)
    algorithm_name <- "posterior_LDA"
  } else if (algorithm %in% c("decomp", "decompTumor2Sig")) {
    decomp_res <- predict_decompTumor2Sig(counts_table, signature)
    exposures <- t(do.call(rbind, decomp_res))
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
    algorithm_name <- "decompTumor2Sig"
  } else {
    stop("Type must be lda or decomp")
  }

  # if not model_id provided, fill one in
  if (is.null(model_id)) {
    model_id <-
      paste(algorithm, length(signatures_to_use), "_exp_pred", sep = "")
  }

  # check if model_id is unique and update if not
  if (!is.null(names(result_list(musica)))) {
    if (!is.null(get_result_list_entry(musica, result_name))) {
      if (!is.null(
        names(get_result_list_entry(musica, result_name)@modality))) {
        if (model_id %in% names(get_modality(musica, result_name, modality))) {
          original_id <- model_id
          tag <- 1
          while (model_id %in%
                 names(get_modality(musica, result_name, modality))) {
            if (tag > 1) {
              model_id <- substr(model_id, 1, nchar(model_id) - 2)
            }
            model_id <- paste(model_id, ".", tag, sep = "")
            tag <- tag + 1
          }
          message("model_id ", original_id,
                  " already exists. model_id updated to ", model_id)
        }
      }
    }
  }

  result <- methods::new("result_model",
    signatures = signature,
    exposures = exposures, num_signatures = length(signatures_to_use),
    model_id = model_id, modality = modality
  )

  # Multiply Weights by sample counts
  sample_counts <- colSums(counts_table)
  matched <- match(colnames(counts_table), names(sample_counts))
  exposures(result) <- sweep(exposures(result), 2, sample_counts[matched],
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
    # Update theta
    theta <- MCMCprecision::fit_dirichlet(x = samp_sig_prob_mat +
                                            .Machine$double.eps)$alpha

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
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modality used for signature discovery
#' @param algorithm Algorithm for signature discovery
#' @param annotation Sample annotation to split results into
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

  # Generate and set result_list
  if (!is.na(annotation)) {
    annot_samples <- samp_annot(musica)$Samples
    annot <- samp_annot(musica)$Tumor_Type
    annot_names <- unique(annot)
    num_annotation <- length(annot_names)
  } else {
    annot_names <- NA
    num_annotation <- 1
  }

  # Define new musicas
  for (i in seq_len(num_annotation)) {
    if (!is.na(annotation)) {
      if (verbose) {
        cat(paste("Current Annotation: ", annot_names[i], "\n", sep = ""))
      }
      cur_ind <- which(annot == annot_names[i])
      cur_annot_samples <- annot_samples[cur_ind]
      cur_annot_variants <- variants(musica)[which(
        variants(musica)$sample %in% cur_annot_samples
      ), ]

      cur_musica <- methods::new("musica",
        variants = cur_annot_variants,
        sample_annotations =
          samp_annot(musica)[cur_ind, ],
        count_tables = .subset_count_tables(
          musica,
          cur_annot_samples
        )
      )
    } else {
      cur_musica <- musica
      cur_annot_samples <- unique(variants(musica)$sample)
    }
    # Used for reconstruction error
    cur_counts <- .extract_count_table(cur_musica, modality)

    # Define new results
    for (cur_k in k_start:k_end) {
      discover_signatures(
        musica = cur_musica, modality = modality,
        num_signatures = cur_k, algorithm = algorithm,
        result_name = result_name, nstart = n_start,
        seed = seed, par_cores = par_cores
      )

      cur_model_name <- paste(algorithm, cur_k, sep = "")

      cur_result <- get_model(cur_musica, result_name, modality, cur_model_name)

      recon_error <- mean(vapply(seq_len(ncol(cur_counts)), function(x) {
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
                              make_copy = FALSE, table_name = NULL) {
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
      model_id = model_id
    )
  } else {
    available_annotations <- setdiff(
      colnames(samp_annot(musica)),
      "Samples"
    )
    if (!sample_annotation %in% available_annotations) {
      stop(paste0(
        "Sample annotation ", sample_annotation, " not found, ",
        "available annotations: ", available_annotations
      ))
    }
    annot <- unique(samp_annot(musica)[[sample_annotation]])
    for (i in seq_along(annot)) {
      if (verbose) {
        message(as.character(annot[i]))
      }
      current_musica <- subset_musica_by_annotation(
        musica = musica, annot_col =
          sample_annotation,
        annot_names = annot[i]
      )
      current_musica <- auto_subset_sigs(
        musica = current_musica, modality = modality,
        signature_res = signature_res, min_exists = min_exists,
        proportion_samples = proportion_samples,
        rare_exposure = rare_exposure, algorithm = algorithm,
        result_name = result_name, model_id = annot[i]
      )

      cur_model_name <- annot[i]

      cur_result <- get_model(current_musica, result_name, modality,
                              cur_model_name)

      musica@result_list[[result_name]]@
        modality[[modality]][[cur_model_name]] <- cur_result
    }
  }
  if (combine_res) {
    musica <- combine_predict_grid(musica, modality, signature_res, annot,
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
#' @return Returns new \code{\linkS4class{musica}} object with the results.
#' @keywords internal
auto_subset_sigs <- function(musica, modality, signature_res, algorithm,
                             min_exists = 0.05, proportion_samples = 0.25,
                             rare_exposure = 0.4, result_name = "result",
                             model_id = NULL) {
  test_predicted <- predict_exposure(
    musica = musica, modality = modality,
    signature_res = signature_res,
    algorithm = algorithm,
    model_id = "temp", make_copy = TRUE
  )

  test_predicted_model <-
    get_model(test_predicted, result_name, modality, "temp")

  exposures <- exposures(test_predicted_model)
  num_samples <- ncol(exposures)
  exposures <- sweep(exposures, 2, colSums(exposures), "/")
  to_use <- as.numeric(which(apply(exposures, 1, function(x) {
    sum(x > min_exists) / num_samples
  }) > proportion_samples |
    apply(exposures, 1, max) > rare_exposure))

  final_inferred <- predict_exposure(
    musica = musica, modality = modality,
    signature_res = signature_res,
    signatures_to_use = to_use,
    algorithm = algorithm,
    result_name = result_name,
    model_id = model_id,
    make_copy = TRUE
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

  musica@result_list[[result_name]]@modality[[modality]] <-
    get_modality(musica, result_name, modality)[!names(
      get_modality(musica, result_name, modality)) %in% model_ids]

  grid_res <- new("result_model",
    exposures = comb,
    signatures = signatures(signature_res)[, sig_names],
    modality = modality
  )

  if (is.null(model_rename)) {
    model_id <- paste(model_ids, collapse = ".")
  } else {
    model_id <- model_rename
  }

  # check if model_id is unique and update if not
  if (!is.null(names(result_list(musica)))) {
    if (!is.null(get_result_list_entry(musica, result_name))) {
      if (!is.null(
        names(get_result_list_entry(musica, result_name)@modality))) {
        if (model_id %in% names(get_modality(musica, result_name, modality))) {
          original_id <- model_id
          tag <- 1
          while (model_id %in%
                 names(get_modality(musica, result_name, modality))) {
            if (tag > 1) {
              model_id <- substr(model_id, 1, nchar(model_id) - 2)
            }
            model_id <- paste(model_id, ".", tag, sep = "")
            tag <- tag + 1
          }
          message("model_id ", original_id,
                  " already exists. model_id updated to ", model_id)
        }
      }
    }
  }

  musica@result_list[[result_name]]@modality[[modality]][[model_id]] <- grid_res

  if (make_copy == FALSE) {
    assign(var_name, musica, envir = parent.frame())
  }

  if (make_copy == TRUE) {
    return(musica)
  }
}
