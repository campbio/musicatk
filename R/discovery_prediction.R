#' @importFrom NMF nmf
NULL

#' @title Discover mutational signatures
#' @description Mutational signatures and exposures will be discovered using
#' methods such as Latent Dirichlet Allocation (lda) or Non-Negative 
#' Matrix Factorization (nmf). These algorithms will deconvolute a matrix of
#' counts for mutation types in each sample to two matrices: 1) a "signature"
#' matrix containing the probability of each mutation type in each sample and
#' 2) an "exposure" matrix containing the estimated counts for each signature
#' in each sample. Before mutational discovery can be performed,
#' variants from samples first need to be stored in a 
#' \code{\linkS4class{musica}} object using the \link{create_musica} function 
#' and mutation count tables need to be created using functions such as
#' \link{build_standard_table}.
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of the table to use for signature discovery. Needs
#' to be the same name supplied to the table building functions such as
#' \link{build_standard_table}.
#' @param num_signatures Number of signatures to discover.
#' @param algorithm Method to use for mutational signature discovery. One of 
#' \code{"lda"} or \code{"nmf"}. Default \code{"lda"}.
#' @param seed Seed to be used for the random number generators in the
#' signature discovery algorithms. Default \code{1}.
#' @param nstart Number of independent random starts used in the mutational
#' signature algorithms. Default \code{10}.
#' @param par_cores Number of parallel cores to use. Only used if
#' \code{method = "nmf"}. Default \code{1}.
#' @return Returns a A \code{\linkS4class{musica_result}} object containing
#' signatures and exposures.
#' @examples
#' data(musica)
#' g <- select_genome("19")
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#' discover_signatures(musica = musica, table_name = "SBS96",
#' num_signatures = 3, algorithm = "lda", seed = 12345, nstart = 1)
#' @export
discover_signatures <- function(musica, table_name, num_signatures,
                                algorithm = "lda", seed = 1, nstart = 10,
                                par_cores = 1) {
  if (!methods::is(musica, "musica")) {
    stop("Input to discover_signatures must be a 'musica' object.")
  }
  algorithm <- match.arg(algorithm, c("lda", "nmf"))
  counts_table <- .extract_count_table(musica, table_name)
  present_samples <- which(colSums(counts_table) > 0)
  counts_table <- counts_table[, present_samples]

  if (algorithm  == "lda") {
    lda_counts_table <- t(counts_table)
    if (is.null(seed)) {
      control <- list(nstart = nstart)
    } else {
      control <- list(seed = (seq_len(nstart) - 1) + seed, nstart = nstart)
    }
    
    lda_out <- topicmodels::LDA(lda_counts_table, num_signatures,
                                control = control)
    
    lda_sigs <- exp(t(lda_out@beta))
    rownames(lda_sigs) <- colnames(lda_counts_table)
    colnames(lda_sigs) <- paste0("Signature", seq_len(num_signatures))

    weights <- t(lda_out@gamma)
    rownames(weights) <- paste0("Signature", seq_len(num_signatures))
    colnames(weights) <- rownames(lda_counts_table)

    result <- methods::new("musica_result", signatures = lda_sigs,
                               table_name = table_name,
                               exposures = weights, algorithm = "LDA", 
                           musica = musica)
    exposures(result) <- sweep(exposures(result), 2, colSums(exposures(result)),
                              FUN = "/")
  } else if (algorithm  == "nmf") {
    #Needed to prevent error with entirely zero rows
    epsilon <- 0.00000001

    decomp <- NMF::nmf(counts_table + epsilon, num_signatures, seed = seed,
                       nrun = nstart, .options = paste("p", par_cores,
                                                       sep = ""))
    
    rownames(decomp@fit@H) <- paste("Signature", seq_len(num_signatures),
                                 sep = "")
    colnames(decomp@fit@W) <- paste("Signature", seq_len(num_signatures),
                                    sep = "")
    result <- methods::new("musica_result", signatures = decomp@fit@W,
                               table_name = table_name,
                               exposures = decomp@fit@H, algorithm = "NMF",
                               musica = musica)
    signatures(result) <- sweep(signatures(result), 2, 
                                colSums(signatures(result)), FUN = "/")
    exposures(result) <- sweep(exposures(result), 2, colSums(exposures(result)),
                              FUN = "/")
  } else{
    stop("That method is not supported. Please select 'lda' or 'nmf' ",
         "to generate signatures.")
  }
  # Multiply Weights by sample counts
  sample_counts <- colSums(counts_table)
  matched <- match(colnames(counts_table), names(sample_counts))
  exposures(result) <- sweep(exposures(result), 2, sample_counts[matched],
                            FUN = "*")
  return(result)
}

#' @title Prediction of exposures in new samples using pre-existing signatures
#' @description Exposures for samples will be predicted using an existing set
#' of signatures stored in a \code{\linkS4class{musica_result}} object. 
#' Algorithms available for prediction include a modify version of \code{"lda"},
#' \code{"decompTumor2Sig"}, and \code{"deconstructSigs"}.
#' @param musica A \code{\linkS4class{musica}} object.
#' @param g A \linkS4class{BSgenome} object indicating which genome
#' reference the variants and their coordinates were derived from. Only used
#' if \code{algorithm = "deconstructSigs"}
#' @param table_name Name of table used for posterior prediction.
#' Must match the table type used to generate the prediction signatures
#' @param signature_res Signatures used to predict exposures for the samples
#' \code{musica} object. Existing signatures need to stored in a
#' \code{\linkS4class{musica_result}} object.
#' @param algorithm Algorithm to use for prediction of exposures. One of
#' \code{"lda"}, \code{"decompTumor2Sig"}, or
#' \code{"deconstructSigs"}.
#' @param signatures_to_use Which signatures in the \code{signature_res} result
#' object to use. Default is to use all signatures.
#' @param verbose If \code{TRUE}, progress will be printing. Only used if
#' \code{algorithm = "lda"}. Default \code{FALSE}.
#' @return Returns a A \code{\linkS4class{musica_result}} object containing
#' signatures given by the \code{signature_res} parameter and exposures
#' predicted from these signatures.
#' @examples
#' data(musica)
#' data(cosmic_v2_sigs)
#' g <- select_genome("19")
#' build_standard_table(musica, g, "SBS96", overwrite = TRUE)
#' result <- predict_exposure(musica = musica, table_name = "SBS96",
#' signature_res = cosmic_v2_sigs, algorithm = "lda")
#' 
#' # Predict using LDA-like algorithm with seed set to 1
#' set.seed(1)
#' predict_exposure(musica = musica, table_name = "SBS96",
#' signature_res = cosmic_v2_sigs, algorithm = "lda")
#' @export
predict_exposure <- function(musica, g, table_name, signature_res, 
                             algorithm = c("lda", "decompTumor2Sig", 
                                           "deconstructSigs"),
                             signatures_to_use = seq_len(ncol(
                               signatures(signature_res))), verbose = FALSE) {
  algorithm <- match.arg(algorithm)
  signature <- signatures(signature_res)[, signatures_to_use, drop = FALSE]
  counts_table <- .extract_count_table(musica, table_name)
  present_samples <- which(colSums(counts_table) > 0)
  counts_table <- counts_table[, present_samples, drop = FALSE]

  if (algorithm %in% c("lda_posterior", "lda", "lda_post")) {
     lda_res <- lda_posterior(counts_table = counts_table, signature = 
                                signature, max.iter = 100, verbose = verbose)
    exposures <- t(lda_res$samp_sig_prob_mat)
    algorithm_name <- "posterior_LDA"
  }else if (algorithm %in% c("decomp", "decompTumor2Sig")) {
    decomp_res <- predict_decompTumor2Sig(counts_table, signature)
    exposures <- t(do.call(rbind, decomp_res))
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
    algorithm_name <- "decompTumor2Sig"
  }else if (algorithm %in% c("ds", "deconstruct", "deconstructSigs")) {
    sigs.input <- deconstructSigs::mut.to.sigs.input(mut.ref = variants(musica),
                              sample.id = "sample", chr = "chr", pos = "start", 
                              ref = "ref", alt = "alt", bsg = g)
    sig_all <- t(signature)
    middle <- unlist(lapply(strsplit(colnames(sig_all), "_"), "[", 1))
    context <- lapply(strsplit(colnames(sig_all), "_"), "[", 2)
    first <- unlist(lapply(context, substr, 1, 1))
    last <- unlist(lapply(context, substr, 3, 3))
    new_cols <- paste(first, "[", middle, "]", last, sep = "")
    colnames(sig_all) <- new_cols

    ds_res <- vapply(rownames(sigs.input), function(x) {
      ds_result <- whichSignatures(tumor_ref = sigs.input,
                                   contexts_needed = TRUE,
                                   signatures_limit = ncol(signature),
                                   tri_counts_method = "default",
                                   sample_id = x, signatures_ref = sig_all)
      return(as.matrix(ds_result$weights))
    }, FUN.VALUE = rep(0, ncol(signature)))
    exposures <- ds_res
    colnames(exposures) <- colnames(counts_table)
    rownames(exposures) <- colnames(signature)
    algorithm_name <- "deconstructSigs"
  } else {
    stop("Type must be lda or decomp")
  }
  result <- methods::new("musica_result", signatures = signature, 
                         exposures = exposures, algorithm = algorithm_name, 
                         musica = musica, table_name = table_name)

  # Multiply Weights by sample counts
  sample_counts <- colSums(counts_table)
  matched <- match(colnames(counts_table), names(sample_counts))
  exposures(result) <- sweep(exposures(result), 2, sample_counts[matched], 
                             FUN = "*")
  return(result)
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
    sig_mut_counts[s, ] <- base::tabulate(sample(x = seq_len(k), size =
                                                   sample_count_sums[s],
                                                 replace = TRUE), k)
  }

  # Update signature proportion matrix
  if (verbose) {
    message("Calculating Signature Proportions")
  }
  for (i in seq_len(max.iter)) {
    for (s in seq_len(num_samples)) {
      #updating each mutation probability to reassign to a signature
      log_prob_mut_reassignment <-
        digamma(sig_mut_counts[s, ] + theta) -
        digamma(sample_count_sums[s] + sum(theta))
      #updating present sample topic probability
      sig_sample_weights <- t(signature + 1e-20) *
        exp(log_prob_mut_reassignment) # avoid 0 in norm
      sig_sample_weights <- sweep(sig_sample_weights, MARGIN = 2, STATS =
                                    colSums(sig_sample_weights), FUN = "/")
      #assigned counts for a topic for a sample
      updated_topic_motifs <- counts_table[, s] * t(sig_sample_weights)

      # Update nN.SbyT[s, ] sample counts assigned to signature
      sig_mut_counts[s, ] <- colSums(updated_topic_motifs)

      # Update p.SbyT[s, ]
      samp_sig_prob_mat[s, ] <- (sig_mut_counts[s, ]) / (sample_count_sums[s])
    }
    # Update theta
    theta <- MCMCprecision::fit_dirichlet(x = samp_sig_prob_mat + 
                                            .Machine$double.eps)$alpha
    #if(any(is.nan(theta))){
    #  browser()
    #}
    if (verbose) {
      message(theta)
    }
  }
  return(list(samp_sig_prob_mat = samp_sig_prob_mat, theta.poster = theta))
}

predict_decompTumor2Sig <- function(sample_mat, signature_mat) {
  #Alexandrov-type prediction
  input_signatures_normalized <- apply(signature_mat, 2,
                                       function(x) {x / sum(x)})
  signatures <- split(input_signatures_normalized,
                      col(input_signatures_normalized))
  #signatures_ref <- decompTumor2Sig::readAlexandrovSignatures()
  ns <- as.matrix(row.names(signature_mat))
  ns <- apply(ns, 1, function(x) {
    stringr::str_c(substr(x, 5, 5), "[", substr(x, 1, 3), "]",
                   substr(x, 7, 7))})
  signatures <- lapply(signatures, setNames, ns)

  input_samples_normalized <- apply(sample_mat, 2, function(x) {x / sum(x)})
  input_samples1 <- split(input_samples_normalized,
                          col(input_samples_normalized))
  genomes <- lapply(input_samples1, setNames, ns)

  sample_weight_mat <- decompTumor2Sig::decomposeTumorGenomes(genomes,
                                                              signatures,
                                                              verbose = FALSE)
  return(sample_weight_mat)
}

#placeholder
.multi_modal_discovery <- function(musica, num_signatures, motif96_name,
                                  rflank_name, lflank_name, max.iter = 125) {
  motif96 <- .extract_count_table(musica, motif96_name)
  rflank <- .extract_count_table(musica, rflank_name)
  lflank <- .extract_count_table(musica, lflank_name)
  message(dim(motif96))
  message(dim(rflank))
  message(dim(lflank))
}

whichSignatures <- function(tumor_ref = NA,
                           sample_id,
                           signatures_ref,
                           associated = c(),
                           signatures_limit = NA,
                           signature_cutoff = 0.06,
                           contexts_needed = FALSE,
                           tri_counts_method = "default") {
  if (is(tumor_ref, 'matrix')) {
    stop(paste("Input tumor.ref needs to be a data frame or location of ", 
               "input text file", sep = ""))
  }

  if (exists("tumor.ref", mode = "list") | is(tumor_ref, "data.frame")) {
    tumor     <- tumor_ref
    if(contexts_needed == TRUE) {
      tumor   <- deconstructSigs::getTriContextFraction(mut.counts.ref = tumor,
                                                        trimer.counts.method =
                                                          tri_counts_method)
    }
  } else {
    if (file.exists(tumor_ref)) {
      tumor   <- utils::read.table(tumor_ref, sep = "\t", header = TRUE,
                                   as.is = TRUE, check.names = FALSE)
      if (contexts_needed == TRUE) {
        tumor <- deconstructSigs::getTriContextFraction(tumor,
                                                        trimer.counts.method =
                                                          tri_counts_method)
      }
    } else {
      message("tumor.ref is neither a file nor a loaded data frame")
    }
  }

  if (missing(sample_id) && nrow(tumor) == 1) {
    sample_id <- rownames(tumor)[1]
  }
  # Take patient id given
  tumor <- as.matrix(tumor)
  if (!sample_id %in% rownames(tumor)) {
    stop(paste(sample_id, " not found in rownames of tumor.ref", sep = ""))
  }
  tumor <- subset(tumor, rownames(tumor) == sample_id)
  if (round(rowSums(tumor), digits = 1) != 1) {
    stop(paste0("Sample: ", sample_id, " is not normalized. Consider using ", 
    "contexts.needed = TRUE", sep = " "))
  }
  signatures <- signatures_ref

  signatures    <- as.matrix(signatures)
  original_sigs <- signatures

  # Check column names are formatted correctly
  if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < 
      length(colnames(signatures))) {
    colnames(tumor) <- deconstructSigs::changeColumnNames(colnames(tumor))
    if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < 
        length(colnames(signatures))) {
      stop("Check column names on input file")
    }
  }

  # Ensure that columns in tumor match the order of those in signatures
  tumor <- tumor[, colnames(signatures), drop = FALSE]

  #Take a subset of the signatures
  if (!is.null(associated)) {
    signatures <- signatures[rownames(signatures) %in% associated, ]
  }

  if (is.na(signatures_limit)) {
    signatures_limit <- nrow(signatures)
  }

  #Set the weights matrix to 0
  weights         <- matrix(0, nrow = nrow(tumor), ncol = nrow(signatures),
                            dimnames = list(rownames(tumor),
                                           rownames(signatures)))

  seed            <- deconstructSigs::findSeed(tumor, signatures)
  weights[seed]   <- 1
  w               <- weights * 10

  error_diff      <- Inf
  error_threshold <- 1e-3

  num <- 0
  while (error_diff > error_threshold) {
    num        <- num + 1
    #message(num)
    error_pre  <- deconstructSigs::getError(tumor, signatures, w)
    if (error_pre == 0) {
      break
    }
    w          <- deconstructSigs::updateW_GR(tumor, signatures, w,
                                              signatures.limit =
                                                signatures_limit)
    error_post <- deconstructSigs::getError(tumor, signatures, w)
    error_diff <- (error_pre - error_post) / error_pre
  }

  weights <- w / sum(w)
  unknown <- 0

  ## filtering on a given threshold value (0.06 default)
  weights[weights < signature_cutoff] <- 0
  unknown <- 1 - sum(weights)

  product <- weights %*% signatures
  diff    <- tumor - product

  x       <- matrix(data = 0, nrow = 1, ncol = nrow(original_sigs),
                    dimnames = list(rownames(weights), rownames(original_sigs)))
  x       <- data.frame(x)
  x[colnames(weights)] <- weights
  weights <- x

  out        <- list(weights, tumor, product, diff, unknown)
  names(out) <- c("weights", "tumor", "product", "diff", "unknown")
  return(out)
}

#' Generate result_grid from musica based on annotation and range of k
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of table used for signature discovery
#' @param algorithm Algorithm for signature discovery
#' @param annotation Sample annotation to split results into
#' @param k_start Lower range of number of signatures for discovery
#' @param k_end Upper range of number of signatures for discovery
#' @param n_start Number of times to discover signatures and compare based on
#' posterior loglikihood
#' @param seed Seed to use for reproducible results, set to null to disable
#' @param par_cores Number of parallel cores to use (NMF only)
#' @param verbose Whether to output loop iterations
#' @return A result object containing signatures and sample weights
#' @examples
#' data(musica_sbs96)
#' grid <- generate_result_grid(musica_sbs96, "SBS96", "lda", k_start = 2, 
#' k_end = 5)
#' @export
generate_result_grid <- function(musica, table_name, algorithm = "lda",
                                 annotation = NA, k_start, k_end, n_start = 1,
                                 seed = NULL, par_cores = FALSE,
                                 verbose = FALSE) {
  result_grid <- methods::new("musica_result_grid")

  #Set Parameters
  params <- data.table::data.table("discovery_type" = algorithm,
                                   "annotation_used" = annotation, "k_start" =
                                     k_start, "k_end" = k_end,
                                   "total_num_samples" =
                                     nrow(samp_annot(musica)),
                                   "nstart" = n_start, seed = seed)
  set_grid_params(result_grid, params)

  #Initialize grid_table and result_list
  grid_table <- data.table::data.table(annotation = character(), k =
                                         numeric(), num_samples = numeric(),
                                       reconstruction_error = numeric())
  result_list <- list()
  list_elem <- 1

  #Generate and set result_list
  if (!is.na(annotation)) {
    annot_samples <- samp_annot(musica)$Samples
    annot <- samp_annot(musica)$Tumor_Type
    annot_names <- unique(annot)
    num_annotation <- length(annot_names)
  } else {
    annot_names <- NA
    num_annotation <- 1
  }

  #Define new musicas
  for (i in seq_len(num_annotation)) {
    if (!is.na(annotation)) {
      if (verbose) {
        cat(paste("Current Annotation: ", annot_names[i], "\n", sep = ""))
      }
      cur_ind <- which(annot == annot_names[i])
      cur_annot_samples <- annot_samples[cur_ind]
      cur_annot_variants <- variants(musica)[which(
        variants(musica)$sample %in% cur_annot_samples), ]

      cur_musica <- methods::new("musica", variants = cur_annot_variants,
                       sample_annotations =
                         samp_annot(musica)[cur_ind, ],
                       count_tables = .subset_count_tables(musica, 
                                                           cur_annot_samples))
    } else {
      cur_musica <- musica
      cur_annot_samples <- unique(variants(musica)$sample)
    }
    #Used for reconstruction error
    cur_counts <- .extract_count_table(cur_musica, table_name)

    #Define new results
    for (cur_k in k_start:k_end) {
      cur_result <- discover_signatures(musica = cur_musica, table_name =
                                          table_name, num_signatures = cur_k,
                                        algorithm = algorithm, nstart =
                                          n_start, seed = seed, par_cores =
                                          par_cores)
      result_list[[list_elem]] <- cur_result
      list_elem <- list_elem + 1

      recon_error <- mean(vapply(seq_len(ncol(cur_counts)), function(x)
        mean((cur_counts[, x, drop = FALSE] -
                reconstruct_sample(cur_result, x))^2), FUN.VALUE = 0)^2)

      grid_table <- rbind(grid_table, data.table::data.table(
        annotation = annot_names[i], k = cur_k, num_samples =
          length(cur_annot_samples), reconstruction_error = recon_error))
    }
  }
  set_grid_list(result_grid, result_list)
  set_grid_table(result_grid, grid_table)
  return(result_grid)
}

reconstruct_sample <- function(result, sample_number) {
  reconstruction <- matrix(apply(sweep(signatures(result), 2,
                                       exposures(result)[, sample_number,
                                                      drop = FALSE], FUN = "*"),
                                 1, sum), dimnames = list(
                                   rownames(signatures(result)), 
                                   "Reconstructed"))
  return(reconstruction)
}

#' Automatic filtering of signatures for exposure prediction gridded across
#' specific annotation
#'
#' @param musica Input samples to predit signature weights
#' @param table_name Name of table used for posterior prediction (e.g. SBS96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", decompTumor2Sig, and deconstructSigs
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
#' @return A list of results, one per unique annotation value, if no
#' annotation value is given, returns a single result for all samples, or
#' combines into a single result if combines_res = TRUE
#' @examples
#' data(musica_annot)
#' data(cosmic_v2_sigs)
#' auto_predict_grid(musica = musica_annot, table_name = "SBS96",
#' signature_res = cosmic_v2_sigs, algorithm = "lda",
#' sample_annotation = "Tumor_Subtypes")
#' auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda")
#' @export
auto_predict_grid <- function(musica, table_name, signature_res, algorithm,
                              sample_annotation = NULL, min_exists = 0.05,
                              proportion_samples = 0.25, rare_exposure = 0.4,
                              verbose = TRUE, combine_res = TRUE) {
  if (is.null(sample_annotation)) {
    combine_res <- FALSE
    result <- auto_subset_sigs(musica = musica, table_name =
                       table_name, signature_res =
                       signature_res, algorithm = algorithm,
                       min_exists = min_exists, proportion_samples =
                       proportion_samples, rare_exposure =
                       rare_exposure)
  } else {
    available_annotations <- setdiff(colnames(samp_annot(musica)),
                                     "Samples")
    if (!sample_annotation %in% available_annotations) {
      stop(paste0("Sample annotation ", sample_annotation, " not found, ",
                 "available annotations: ", available_annotations))
    }
    annot <- unique(samp_annot(musica)[[sample_annotation]])
    result <- list()
    for (i in seq_along(annot)) {
      if (verbose) {
        message(as.character(annot[i]))
      }
      current_musica <- subset_musica_by_annotation(musica = musica, annot_col =
                                                    sample_annotation,
                                                  annot_names = annot[i])
      current_predicted <- auto_subset_sigs(musica = current_musica, table_name =
                                              table_name, signature_res =
                                              signature_res, min_exists =
                                              min_exists, proportion_samples =
                                              proportion_samples,
                                            rare_exposure = rare_exposure,
                                            algorithm = algorithm)
      result[[as.character(annot[i])]] <- current_predicted
    }
  }
  if (combine_res) {
    result <- combine_predict_grid(result, musica, table_name, signature_res)
  }
  return(result)
}

#' Automatic filtering of inactive signatures
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of table used for posterior prediction (e.g. SBS96)
#' @param signature_res Signatures to automatically subset from for prediction
#' @param algorithm Algorithm to use for prediction. Choose from
#' "lda_posterior", decompTumor2Sig, and deconstructSigs
#' @param min_exists Threshold to consider a signature active in a sample
#' @param proportion_samples Threshold of samples to consider a signature
#' active in the cohort
#' @param rare_exposure A sample will be considered active in the cohort if at
#' least one sample has more than this threshold proportion
#' @return A result object containing automatically subset signatures
#' and corresponding sample weights
#' @keywords internal
auto_subset_sigs <- function(musica, table_name, signature_res, algorithm,
                             min_exists = 0.05, proportion_samples = 0.25,
                             rare_exposure = 0.4) {
  test_predicted <- predict_exposure(musica = musica, table_name = table_name,
                                    signature_res = signature_res,
                                    algorithm = algorithm)
  exposures <- exposures(test_predicted)
  num_samples <- ncol(exposures)
  exposures <- sweep(exposures, 2, colSums(exposures), "/")
  to_use <- as.numeric(which(apply(exposures, 1, function(x)
    sum(x > min_exists) / num_samples) > proportion_samples |
      apply(exposures, 1, max) > rare_exposure))
  final_inferred <- predict_exposure(musica = musica, table_name = table_name,
                                     signature_res = signature_res,
                                     signatures_to_use = to_use,
                                     algorithm = algorithm)
  return(final_inferred)
}

#' Combine prediction grid list into a result object. Exposure values are zero
#' for samples in an annotation where that signature was not predicted
#'
#' @param grid_list A list of result objects from the prediction grid to
#' combine into a single result
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Table name used for prediction
#' @param signature_res Signatures to automatically subset from for prediction
#' @return A result object combining all samples and signatures from a
#' prediction grid. Samples have zero exposure value for signatures not found
#' in that annotation type.
#' @examples
#' data(musica_annot)
#' data(cosmic_v2_sigs)
#' grid <- auto_predict_grid(musica_annot, "SBS96", cosmic_v2_sigs, "lda",
#' "Tumor_Subtypes", combine_res = FALSE)
#' combined <- combine_predict_grid(grid, musica_annot, "SBS96", cosmic_v2_sigs)
#' plot_exposures(combined, group_by = "annotation", 
#' annotation="Tumor_Subtypes")
#' @export
combine_predict_grid <- function(grid_list, musica, table_name, signature_res) {
  sig_names <- NULL
  for (i in seq_len(length(grid_list))) {
    sig_names <- c(sig_names, rownames(grid_list[[i]]@exposures))
  }
  sig_names <- unique(sig_names)
  sig_names <- sig_names[order(sig_names)]

  comb <- NULL
  for (i in seq_len(length(grid_list))) {
    if (!table_selected(grid_list[[i]]) %in% table_name) {
      stop("Result number: ", i, " was not in selected table_name: ", 
           table_name)
    }
    samp <- exposures(grid_list[[i]])
    missing <- sig_names[!sig_names %in% rownames(samp)]
    missing_mat <- matrix(0, length(missing), ncol(samp))
    rownames(missing_mat) <- missing
    samp <- rbind(samp, missing_mat)
    samp <- samp[order(rownames(samp)), , drop = FALSE]
    comb <- cbind(comb, samp)
  }
  grid_res <- new("musica_result", musica = musica, exposures = comb,
                  signatures = signatures(signature_res)[, sig_names], 
                  table_name = table_name)
}
