#' @importFrom lsei nnls

#' @title Create a full_benchmark object 
#' @description Initialize a \code{\linkS4class{full_benchmark}} object for benchmarking
#'
#' @param true_signatures A matrix of true signatures by mutational motifs
#' @param true_exposures A matrix or data.frame of true signature weights,
#' with samples as rows and signatures as columns (signatures as rows and
#' samples as columns is also accepted). Sample names must match the column
#' names of \code{count_table}, and signature names the column names of
#' \code{true_signatures}; they are matched by name, not position.
#' @param count_table Summary table with per-sample un-normalized motif counts
#' @param variant_class The count table: \code{"SBS96"} or \code{"DBS78"}.
#' Aliases accepted by \link{create_musica_from_counts}, such as
#' \code{"SBS"} or \code{"DBS"}, are stored under the table name, which is
#' the \code{prediction_modality} to give to \link{benchmark}.
#'
#' @return A \code{\linkS4class{full_benchmark}} object
#' @examples
#' data(synthetic_breast_counts)
#' data(synthetic_breast_true_exposures)
#' true_sigs <- c("SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26")
#' true_signatures <- signatures(cosmic_v2_sigs)[,true_sigs]
#' full_benchmark <- create_benchmark(true_signatures, 
#'   synthetic_breast_true_exposures, synthetic_breast_counts, "SBS96")
#' 
#' @export
create_benchmark <- function(true_signatures, true_exposures, count_table, variant_class){
  
  modality <- .benchmark_modality(variant_class)
  musica_truth <- create_musica_from_counts(count_table, modality)
  true_signatures <- as.matrix(true_signatures)
  # as.matrix() drops automatic data.frame row names, which are sample names
  exposure_matrix <- as.matrix(true_exposures)
  if (is.data.frame(true_exposures)) {
    rownames(exposure_matrix) <- rownames(true_exposures)
  }
  true_exposures <- .orient_true_exposures(
    exposure_matrix, colnames(true_signatures), sample_names(musica_truth)
  )
  
  add_result(true_signatures, true_exposures, musica_truth, "ground_truth",
             modality, "ground_truth")
  
  # create benchmark object
  full_benchmark <- new("full_benchmark", ground_truth = musica_truth)
  
  return(full_benchmark)
  
}

# The count-table name for a benchmark 'variant_class': "SBS96" or "DBS78"
.benchmark_modality <- function(variant_class) {
  sbs <- c("snv", "SNV", "SNV96", "SBS", "SBS96")
  dbs <- c("DBS", "dbs", "doublet", "DBS78")
  if (is.character(variant_class) && length(variant_class) == 1) {
    if (variant_class %in% sbs) {
      return("SBS96")
    }
    if (variant_class %in% dbs) {
      return("DBS78")
    }
  }
  stop(
    "'variant_class' must be \"SBS96\" or \"DBS78\" (or an alias such as ",
    "\"SBS\" or \"DBS\")."
  )
}

# True exposures as signatures x samples, rows and columns matched by name
# to the true signatures and the count-table samples
.orient_true_exposures <- function(true_exposures, signature_names, samples) {
  if (setequal(rownames(true_exposures), samples)) {
    true_exposures <- t(true_exposures)
  }
  sample_match <- anyDuplicated(colnames(true_exposures)) == 0 &&
    setequal(colnames(true_exposures), samples)
  if (!sample_match) {
    stop(
      "The sample names of 'true_exposures' must match the column names ",
      "of 'count_table'."
    )
  }
  signature_match <- anyDuplicated(rownames(true_exposures)) == 0 &&
    setequal(rownames(true_exposures), signature_names)
  if (!signature_match) {
    stop(
      "The signature names of 'true_exposures' must match the column ",
      "names of 'true_signatures'."
    )
  }
  return(true_exposures[signature_names, samples, drop = FALSE])
}

# True exposures of one signature for 'samples', matched by sample name
.true_exposure <- function(result_true, signature, samples) {
  true_exposures <- exposures(result_true)
  if (!signature %in% rownames(true_exposures)) {
    stop(
      "The ground truth has no exposures for signature '", signature, "' ",
      "in its rows. A full_benchmark made with an older version of ",
      "musicatk stores them the other way round; recreate it with ",
      "create_benchmark()."
    )
  }
  return(as.numeric(true_exposures[signature, as.character(samples)]))
}

#' @title Run the benchmark framework on a prediction
#' @description Perform benchmarking on a signature discovery prediction compared
#' to a ground truth. Potential errors in the predicted signatures, such as
#' composite or duplicate signatures are adjusted, and a summary of the accuracy
#' of the prediction is given.
#'
#' @param full_benchmark An object of class \code{\linkS4class{full_benchmark}}
#' created with the \link{create_benchmark} function or returned from a previous
#' \code{benchmark} run.
#' @param prediction_musica An object of class \code{\linkS4class{musica}}
#' containing the predicted signatures and exposures to benchmark.
#' @param prediction_model_id The model id of the prediction to be located within
#' prediction_musica
#' @param prediction_modality The modality of the prediction to benchmark,
#' \code{"SBS96"} or \code{"DBS78"}. It must be the modality of the ground
#' truth (the table name given by \link{create_benchmark}).
#' @param prediction_result_name Name of the result list entry prediction to
#' benchmark. Default \code{"result"}.
#' @param method_id An identifier for the prediction being benchmarked. If not
#' supplied, it will be automatically set to the variable name of the prediction
#' provided. Default \code{NULL}.
#' @param threshold Cosine similarity cutoff for comparing predicted and true
#' signatures. Default \code{0.8}.
#' @param adjustment_threshold Cosine similarity value of high confidence.
#' Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold. Default \code{0.9}.
#' @param re_method Method for reconstruction error calculation, either
#' proportion or raw. Default \code{"prop"}.
#' @param description Further details about the prediction being benchmarked.
#' Default \code{NULL}.
#' @param plot If \code{FALSE}, plots will be suppressed. Default \code{TRUE}.
#' @param make_copy If \code{TRUE}, the \code{full_benchmark} object provided
#' will not be modified and a new object will be returned. If \code{FALSE}, the
#' object provided will be modified and nothing will be returned. Default
#' \code{FALSE}.
#'
#' @return If \code{make_copy == TRUE}, a new \code{full_benchmark} object is
#' returned. If \code{make_copy == FALSE}, nothing is returned.
#' @examples
#' data(synthetic_breast_counts)
#' data(synthetic_breast_true_exposures)
#' data(cosmic_v2_sigs)
#' true_sigs <- c("SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26")
#' true_signatures <- signatures(cosmic_v2_sigs)[,true_sigs]
#' full_benchmark <- create_benchmark(true_signatures, 
#'   synthetic_breast_true_exposures, synthetic_breast_counts, "SBS96")
#' prediction_musica <- create_musica_from_counts(synthetic_breast_counts, "SBS96")
#' discover_signatures(prediction_musica, "SBS96", 8, "nmf")
#' benchmark(full_benchmark, prediction_musica, "nmf8", "SBS96")
#' 
#' data(synthetic_breast_counts)
#' data(synthetic_breast_true_exposures)
#' data(example_predicted_sigs)
#' data(example_predicted_exp)
#' true_sigs <- c("SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26")
#' true_signatures <- signatures(cosmic_v2_sigs)[,true_sigs]
#' full_benchmark <- create_benchmark(true_signatures, 
#'   synthetic_breast_true_exposures, synthetic_breast_counts, "SBS96")
#' prediction_musica <- create_musica_from_counts(synthetic_breast_counts, "SBS96")
#' add_result(example_predicted_sigs, example_predicted_exp, prediction_musica, 
#'   "result", "SBS96", "example_k8")
#' benchmark(full_benchmark, prediction_musica, "example_k8", "SBS96")
#' 
#' @export
benchmark <- function(full_benchmark, prediction_musica, prediction_model_id, 
                      prediction_modality, prediction_result_name = "result", 
                      method_id = NULL, threshold = 0.8, 
                      adjustment_threshold = 0.9, re_method = "prop", 
                      description = NULL, plot = TRUE, make_copy = FALSE){
  
  # save variable name to update full benchmark object
  if (make_copy == FALSE){
    var_name <- deparse(substitute(full_benchmark))
  }
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object. ", 
         "Use function 'create_benchmark' to initialize a 'full_benchmark' object")
  }
  
  # summaries were matrices before musicatk 2.7.18 (ADR 0004)
  if (!is.data.frame(method_view_summary(full_benchmark))) {
    stop(
      "This full_benchmark was made with an older version of musicatk, ",
      "which stored its summaries as matrices. Recreate it with ",
      "create_benchmark()."
    )
  }

  # check that prediction is a result model object
  if (!methods::is(prediction_musica, "musica")) {
    stop("'prediction' must be a 'musica' object.")
  }
  
  # the prediction must be of the ground truth's modality
  truth_musica <- ground_truth(full_benchmark)
  if (is.null(get_modality(truth_musica, "ground_truth",
                           prediction_modality))) {
    stop(
      "'prediction_modality' is '", prediction_modality, "', but the ground ",
      "truth is '", paste(names(tables(truth_musica)), collapse = "', '"),
      "'."
    )
  }

  # extract the result_model object containing the prediction to benchmark
  prediction <- get_model(prediction_musica, prediction_result_name, 
                          prediction_modality, prediction_model_id)
  
  # if no ID provided, try to make one automatically
  if (is.null(method_id)){
    
    # create ID
    method_id <- model_id(prediction)
    
    if(is.null(method_id)){
      method_id <- "prediction"
    }
    
    # display message that ID was automatically generated
    message("No method_id provided, automatically generated method_id: ", method_id)
  }
  
  # check if ID is unique
  if (method_id %in% names(indv_benchmarks(full_benchmark))){
    original_id <- method_id
    
    # update id to be unique
    tag <- 1
    while (method_id %in% names(indv_benchmarks(full_benchmark))){
      if (tag > 1){
        method_id <- substr(method_id, 1, nchar(method_id)-2)
      }
      method_id <- paste(method_id, ".", tag, sep = "")
      tag <- tag + 1
    }
    
    # display message that ID was not unique and was updated
    message("method_id ", original_id, " already exists. method_id updated to ", 
            method_id)
    
  }
  
  if (threshold != 0.8){
    warning("Default threshold overriden. Interpret results with caution if 
            comparing benchmark runs with inconsistent thresholds.")
  }
  
  if (adjustment_threshold != 0.9){
    warning("Default adjustment threshold overriden. Interpret results with caution if 
            comparing benchmark runs with inconsistent thresholds.")
  }
  
  # extract ground truth musica and model
  truth <- truth_musica
  truth_model <- get_model(truth, "ground_truth", prediction_modality,
                           "ground_truth")

  # initial comparison of predicted and true signatures
  message("\nComparing to true signatures (initial)...")
  initial_comparison <- .benchmark_compare(
    prediction, truth_model, threshold, adjustment_threshold
  )

  # correct duplicate signatures
  message("Correcting duplicates...")
  duplicates_corrected <- .correct_duplicates(
    prediction, initial_comparison
  )

  # intermediate comparison of predicted and true signatures
  message("Comparing to true signatures (post duplicates corrected)...")
  duplicates_corrected_comparison <- .benchmark_compare(
    duplicates_corrected, truth_model, threshold, adjustment_threshold
  )

  # correct composite signatures
  message("Correcting composites...")
  composites_corrected <- .correct_composites(
    duplicates_corrected, duplicates_corrected_comparison,
    musica = prediction_musica, truth_model
  )

  # final comparison between predicted and true signatures
  message("Comparing to true signatures (post composites corrected)...")
  composites_corrected_comparison <- .benchmark_compare(
    composites_corrected, truth_model, threshold, adjustment_threshold
  )

  # extract count table
  count_table <- extract_count_tables(prediction_musica)
  count_table <- get_count_table(count_table[[prediction_modality]])
  
  # create the one-row summary data.frame
  message("Creating summary...")
  single_summary <- .generate_summary(
    method_id, prediction, truth_model, initial_comparison, count_table,
    composites_corrected, composites_corrected_comparison, re_method
  )
  
  # create single benchmark object
  message("Creating individual benchmark object...")
  indv_benchmark <- new("single_benchmark", 
                        initial_pred = get_model(prediction_musica, prediction_result_name,
                                                 prediction_modality, prediction_model_id),
                        intermediate_pred = duplicates_corrected, 
                        final_pred = composites_corrected, initial_comparison = initial_comparison,
                        intermediate_comparison = duplicates_corrected_comparison,
                        final_comparison = composites_corrected_comparison, 
                        single_summary = single_summary, method_id = method_id, threshold = threshold,
                        adjustment_threshold = adjustment_threshold, description = description)
  
  # update full benchmark object
  message("Updating full benchmark object...")
  full_benchmark <- .update_benchmark(full_benchmark, indv_benchmark, single_summary, prediction_modality)
  
  # update global variable
  if (make_copy == FALSE){
    assign(var_name, full_benchmark, envir = parent.frame())
  }
  
  # plot all figures
  if (plot == TRUE){
    
    # plots
    message("Generating plots...")
    
    # initial signatures plot
    initial_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Initial")
    print(initial_sig_plot)
    # initial comparison plot
    .plot_stage_comparison(full_benchmark, method_id, "Initial",
                           initial_comparison)

    # Duplicate signature exposures before/afters
    duplicate_plot <- benchmark_plot_duplicate_exposures(full_benchmark, method_id)
    print(duplicate_plot)
    
    # Intermediate signatures plot
    intermediate_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Intermediate")
    print(intermediate_sig_plot)
    # Intermediate comparison plot
    .plot_stage_comparison(full_benchmark, method_id, "Intermediate",
                           duplicates_corrected_comparison)
    
    # Composite signature exposures before/afters
    composite_plot <- benchmark_plot_composite_exposures(full_benchmark, method_id)
    print(composite_plot)
    
    # Final signatures plot
    final_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Final")
    print(final_sig_plot)
    # Final comparison and exposures plots
    if (.plot_stage_comparison(full_benchmark, method_id, "Final",
                               composites_corrected_comparison)) {
      exposure_plot <- benchmark_plot_exposures(full_benchmark, method_id,
                                                prediction = "Final")
      print(exposure_plot)
    }
    
  }
  
  # display the single summary
  print(single_summary)
  
  message("\nDone.\n")
  
  if (make_copy == TRUE){
    return(full_benchmark)
  }
  
}

#' @title Predict and benchmark
#' @description This function will discover signatures from a musica object using a given
#' algorithm and a range of k values. A new discovery is done for each k value.
#' As each discovery is completed, the prediction is benchmarked.
#' 
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality Modalitye to use for signature discovery. Needs
#' to be the same name supplied to the table building functions such as
#' \link{build_standard_table}.
#' @param k_min Minimum number of singatures to predict
#' @param k_max Maximum number of signatures to predict
#' @param full_benchmark An object of class \code{\linkS4class{full_benchmark}}
#' created with the \link{create_benchmark} function or returned from a previous
#' \code{benchmark} run.
#' @param algorithm Method to use for mutational signature discovery. One of 
#' \code{"lda"} or \code{"nmf"}. Default \code{"lda"}.
#' @param result_name Name of the result list entry to save prediction. 
#' Default \code{"result"}.
#' @param threshold Cosine similarity cutoff for comparing predicted and true
#' signatures. Default \code{0.8}.
#' @param adjustment_threshold Cosine similarity value of high confidence.
#' Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold. Default \code{0.9}.
#' @param re_method Method for reconstruction error calculation, either
#' proportion or raw. Default \code{"prop"}.
#' @param plot If \code{FALSE}, plots will be suppressed. Default \code{TRUE}.
#' @param seed Seed to be used for the random number generators in the
#' signature discovery algorithms. Default \code{1}.
#' @param nstart Number of independent random starts used in the mutational
#' signature algorithms. Default \code{10}.
#' @param par_cores Number of parallel cores to use. Only used if
#' \code{method = "nmf"}. Default \code{1}. 
#'
#' @return If \code{make_copy == TRUE}, a new \code{full_benchmark} object is
#' returned. If \code{make_copy == FALSE}, nothing is returned.
#' @examples
#' data(synthetic_breast_counts)
#' data(synthetic_breast_true_exposures)
#' true_sigs <- c("SBS1", "SBS2", "SBS3", "SBS8", "SBS13", "SBS17", "SBS18", "SBS26")
#' true_signatures <- signatures(cosmic_v2_sigs)[,true_sigs]
#' full_benchmark <- create_benchmark(true_signatures, 
#'   synthetic_breast_true_exposures, synthetic_breast_counts, "SBS96")
#' prediction_musica <- create_musica_from_counts(synthetic_breast_counts, "SBS96")
#' predict_and_benchmark(prediction_musica, "SBS96", 7, 9, full_benchmark, "nmf")
#' 
#' @export
predict_and_benchmark <- function(musica, modality, 
                                  k_min, k_max, full_benchmark, algorithm = "lda",
                                  result_name = "result", threshold = 0.8,
                                  adjustment_threshold = 0.9, re_method = "prop",
                                  plot = FALSE, seed = 1, nstart = 10, par_cores = 1){
  
  # loop through k values
  for (k in c(k_min:k_max)){
    
    message("Discovering signatures for k = ", k, "...")
    
    # discover signatures for a k value and save to musica object provided
    discover_signatures(musica = musica, modality = modality,
                        num_signatures = k, algorithm = algorithm,
                        result_name = result_name,
                        seed = seed, nstart = nstart, 
                        par_cores = par_cores)
    
    message("Benchmarking prediction for k = ", k, "...\n")
    
    # benchmark the prediction
    full_benchmark <- benchmark(full_benchmark = full_benchmark, 
                                prediction_musica = musica, 
                                prediction_model_id = paste0(algorithm, k),
                                prediction_modality = modality,
                                prediction_result_name = result_name,
                                method_id = paste0(algorithm, k),
                                threshold = threshold, 
                                adjustment_threshold = adjustment_threshold,
                                re_method = re_method,
                                description = paste0("Prediction from predict_and_benchmark function. Algorithm: ", algorithm, ". K = ", k, "."),
                                plot = plot, make_copy = TRUE)
  }
  
  
  print(method_view_summary(full_benchmark))
  
  return(full_benchmark)
  
}


#' @title Get a single_benchmark object
#' @description Access a \code{\linkS4class{single_benchmark}} object containing information from
#' an individual benchmark run from a \code{\linkS4class{full_benchmark}} object.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object that contains
#' the desired \code{\linkS4class{single_benchmark}} object
#' @param method_id The identifier for the desired \code{\linkS4class{single_benchmark}}
#' object
#'
#' @return A \code{\linkS4class{single_benchmark}} object
#' @examples
#' data(full_benchmark_example)
#' indv_benchmark <- benchmark_get_entry(full_benchmark_example, "example_k8")
#' 
#' @export
benchmark_get_entry <- function(full_benchmark, method_id){
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # extract the single_benchmark
  benchmark <- indv_benchmarks(full_benchmark)[[method_id]]
  
  return(benchmark)
}


#' @title Get a benchmark prediction
#' @description Access a \code{\linkS4class{result_model}} object containing a particular prediction
#' from a \code{\linkS4class{single_benchmark}} object.
#'
#' @param indv_benchmark A \code{\linkS4class{single_benchmark}} object containing the
#' desired prediction. This can be accessed using the \link{benchmark_get_entry}
#' function.
#' @param prediction \code{"Initial"} for the prediction before any benchmarking
#' adjustments have been made, \code{"Intermediate"} for the prediction after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{"Final"} for the prediction at the end of the benchmarking adjustments. 
#'
#' @return A \code{\linkS4class{result_model}} object
#' @examples
#' data(full_benchmark_example)
#' indv_benchmark <- benchmark_get_entry(full_benchmark_example, "example_k8")
#' initial_prediction <- benchmark_get_prediction(indv_benchmark, "Initial")
#' 
#' @export
benchmark_get_prediction <- function(indv_benchmark, prediction){
  
  # "Initial", "Intermediate" or "Final" (aliases accepted)
  prediction <- .match_benchmark_stage(prediction)
  
  # access musica object for desired prediction
  if (prediction == "Initial") {
    result <- initial_pred(indv_benchmark)
  } else if (prediction == "Intermediate") {
    result <- intermediate_pred(indv_benchmark)
  } else {
    result <- final_pred(indv_benchmark)
  }
  return(result)
  
}


#' @title Plot signatures from a benchmarking analysis
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' this function can be used to plot signatures from any step in the benchmarking
#' process. Comparable to the \code{plot_signatures} function but compatible with
#' benchmarking objects.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the signatures to be plotted
#' @param prediction \code{"Initial"} for the signatures before any benchmarking
#' adjustments have been made, \code{"Intermediate"} for the signatures after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{"Final"} for the signatures at the end of the benchmarking adjustments.
#' @param plotly If \code{TRUE}, the the plot will be made interactive
#' using \code{\link[plotly]{plotly}}. Default \code{FALSE}.
#' @param color_variable Name of the column in the variant annotation data.frame
#' to use for coloring the mutation type bars. The variant annotation data.frame 
#' can be found within the count table of the \code{\linkS4class{musica}}
#' object. If \code{NULL}, then the default column specified in the count
#' table will be used. Default \code{NULL}.
#' @param color_mapping A character vector used to map items in the
#' \code{color_variable} to a color. The items in \code{color_mapping}
#' correspond to the colors. The names of the items in \code{color_mapping}
#' should correspond to the uniqeu items in \code{color_variable}. If
#' \code{NULL}, then the default \code{color_mapping} specified in the count
#' table will be used. Default \code{NULL}.
#' @param text_size Size of axis text. Default \code{10}.
#' @param show_x_labels If \code{TRUE}, the labels for the mutation types
#' on the x-axis will be shown. Default \code{TRUE}.
#' @param show_y_labels If \code{TRUE}, the y-axis ticks and labels will be 
#' shown. Default \code{TRUE}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' signature will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each signature. Default \code{FALSE}.
#' @param y_max Vector of maximum y-axis limits for each signature. One value 
#' may also be provided to specify a constant y-axis limit for all signatures.
#' Vector length must be 1 or equivalent to the number of signatures. Default 
#' \code{NULL}.
#' @param annotation Vector of annotations to be displayed in the top right
#' corner of each signature. Vector length must be equivalent to the number of
#' signatures. Default \code{NULL}.
#' @param percent If \code{TRUE}, the y-axis will be represented in percent 
#' format instead of mutation counts. Default \code{TRUE}.
#'
#' @return Generates a ggplot or plotly object
#' @examples
#' data(full_benchmark_example)
#' benchmark_plot_signatures(full_benchmark_example, "example_k8", "Final")
#' 
#' @export
benchmark_plot_signatures <- function(full_benchmark, method_id, prediction,
                                      plotly = FALSE, color_variable = NULL, color_mapping = NULL, text_size = 10,
                                      show_x_labels = TRUE, show_y_labels = TRUE, same_scale = FALSE, y_max = NULL, 
                                      annotation = NULL, percent = TRUE){
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # "Initial", "Intermediate" or "Final" (aliases accepted)
  prediction <- .match_benchmark_stage(prediction)
  
  # access individual benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
    
  # access desired prediction
  result <- benchmark_get_prediction(indv_benchmark, prediction)
  
  # create dummy musica object for proper plotting format
  musica_temp <- ground_truth(full_benchmark)
  add_result(signatures(result), exposures(result), musica_temp, "to_plot",
             modality(result), "to_plot")
  
  # If annotation is NULL, fill in the prediction category
  if(is.null(annotation)){
    num_sigs <- dim(signatures(result))[2]
    annotation <- c(prediction, rep("", num_sigs-1))
  }
  
  # Plot
  signatures_plot <- plot_signatures(musica_temp, "to_plot", modality(result), "to_plot", plotly = plotly, color_variable = color_variable, color_mapping = color_mapping,
                                     text_size = text_size, show_x_labels = show_x_labels, show_y_labels = show_y_labels,
                                     same_scale = same_scale, y_max = y_max, annotation = annotation, percent = percent)
  
  return(signatures_plot)
  
}


#' @title Get benchmark comparison table
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' this function can be used to extract the comparison table between true and
#' predicted signatures from any step in the benchmarking process.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the comparison of interest
#' @param prediction \code{"Initial"} for the comparison before any benchmarking
#' adjustments have been made, \code{"Intermediate"} for the comparison after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{"Final"} for the comparison at the end of the benchmarking adjustments.
#'
#' @return A data.frame containing the comparison between true and predicted
#' signatures
#' @examples
#' data(full_benchmark_example)
#' final_comparison <- benchmark_get_comparison(full_benchmark_example, "example_k8", "Final")
#' 
#' @export
benchmark_get_comparison <- function(full_benchmark, method_id, prediction){
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # "Initial", "Intermediate" or "Final" (aliases accepted)
  prediction <- .match_benchmark_stage(prediction)
  
  # access individual benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # access desired prediction
  if (prediction == "Initial") {
    comparison <- initial_comparison(indv_benchmark)
  } else if (prediction == "Intermediate") {
    comparison <- intermediate_comparison(indv_benchmark)
  } else {
    comparison <- final_comparison(indv_benchmark)
  }
  
  return(comparison)
  
}

#' @title Plot a signature comparison from a benchmarking analysis
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' the comparison between the true and predicted signatures at any step of the
#' benchmarking process can be plotted.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the comparison of interest
#' @param prediction \code{"Initial"} for the comparison before any benchmarking
#' adjustments have been made, \code{"Intermediate"} for the comparison after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{"Final"} for the comparison at the end of the benchmarking adjustments.
#' @param decimals Specifies rounding for similarity metric displayed. Default
#' \code{2}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' comparison will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each comparison. Default \code{FALSE}.
#'
#' @return Returns the comparison plot
#' @examples
#' data(full_benchmark_example)
#' benchmark_plot_comparison(full_benchmark_example, "example_k8", "Final")
#' 
#' @export
benchmark_plot_comparison <- function(full_benchmark, method_id, prediction,
                                      decimals = 2, same_scale = FALSE){
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # "Initial", "Intermediate" or "Final" (aliases accepted)
  prediction <- .match_benchmark_stage(prediction)
  
  comparison <- benchmark_get_comparison(full_benchmark, method_id, prediction)
  .check_comparison_has_matches(comparison, prediction)
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  result <- benchmark_get_prediction(indv_benchmark, prediction)

  truth_musica <- ground_truth(full_benchmark)
  truth <- get_model(truth_musica, "ground_truth", modality(result), "ground_truth")
  
  result_subset <- methods::new("result_model",
                                signatures =
                                  signatures(result)[, comparison$x_sig_index, drop = FALSE],
                                exposures = matrix(),
                                num_signatures = dim(signatures(result))[2],
                                modality = modality(result), model_id = "result_subset"
  )
  other_subset <- methods::new("result_model",
                               signatures =
                                 signatures(truth)[, comparison$y_sig_index, drop = FALSE],
                               exposures = matrix(), modality = modality(truth),
                               model_id = "result_subset"
  )
  
  result_subset_maxes <- NULL
  other_subset_maxes <- NULL
  for (index in seq_len(dim(comparison)[1])) {
    result_subset_maxes <- c(result_subset_maxes,
                             max(signatures(result_subset)[, index]))
  }
  for (index in seq_len(dim(comparison)[1])) {
    other_subset_maxes <- c(other_subset_maxes,
                            max(signatures(other_subset)[, index]))
  }
  maxes <- pmax(result_subset_maxes, other_subset_maxes) * 100
  
  if (same_scale == TRUE) {
    maxes <- rep(max(maxes), length(maxes))
  }
  
  .plot_compare_result_signatures(result_subset, other_subset, comparison,
                                  truth_musica, res1_name = paste(prediction, "Signatures"),
                                  res2_name = "True Signatures",
                                  decimals = decimals, same_scale = same_scale,
                                  maxes = maxes
  )
  
  #return(comparison_plot) #COME BACK TO THIS
  
}

#' @title Plot exposure comparison from a benchmarking analysis
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' the comparison between the true and predicted exposures at any stage of the
#' benchmarking process can be plotted.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#' @param prediction \code{"Initial"} for the exposures before any benchmarking
#' adjustments have been made, \code{"Intermediate"} for the exposures after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{"Final"} for the exposures at the end of the benchmarking adjustments.
#'
#' @return Generates a  ggplot object
#' @examples
#' data(full_benchmark_example)
#' benchmark_plot_exposures(full_benchmark_example, "example_k8", "Final")
#' 
#' @export
benchmark_plot_exposures <- function(full_benchmark, method_id, prediction){
  
  Predicted <- NULL
  True <- NULL
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # "Initial", "Intermediate" or "Final" (aliases accepted)
  prediction <- .match_benchmark_stage(prediction)
  
  # access individual benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # access comparison for desired prediction
  comparison <- benchmark_get_comparison(full_benchmark, method_id, prediction)
  .check_comparison_has_matches(comparison, prediction)
  
  # access desired prediction
  prediction <- benchmark_get_prediction(indv_benchmark, prediction)
  
  # get ground truth
  truth <- ground_truth(full_benchmark)
  truth <- get_model(truth, "ground_truth", modality(prediction), "ground_truth")
  
  # get exposures for all signatures
  predicted <- c()
  true <- c()
  sig <- c()
  index <- 1
  for (true_sig in comparison$y_sig_name){
    predicted_sig <- 
      comparison[index,4]
    predicted <- c(predicted, exposures(prediction)[predicted_sig,])
    true <- c(true, .true_exposure(truth, true_sig,
                                   colnames(exposures(prediction))))
    sig <- c(sig, rep(true_sig, ncol(exposures(prediction))))
    index <- index + 1
  }
  
  # create dataframe for plotting
  plot_df <- data.frame(Predicted = predicted, True = true, Sig = sig)
  
  # plot
  compare_exposure_plot <- ggplot(plot_df, aes(x = Predicted, y = True)) + 
    geom_point(size = 3) + 
    facet_wrap(~Sig, scales = "free") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    theme_classic() + 
    labs(title = "Predicted vs True Activity Levels for Matched Signatures",
         x="Predicted Activity", y = "True Activity") + 
    theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
    geom_abline() + 
    geom_smooth(method = "lm") +
    theme(legend.title=element_blank())
  
  return(compare_exposure_plot)
  
}


#' @title Plot effect of duplicate correction
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' the true and predicted exposures can be plotted both before and after the
#' duplicate signature adjustment. The effect of the adjustment can then be
#' observed.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#'
#' @return A list of ggplot objects
#' @examples
#' data(full_benchmark_example)
#' benchmark_plot_duplicate_exposures(full_benchmark_example, "example_k8")
#' 
#' @export
benchmark_plot_duplicate_exposures <- function(full_benchmark, method_id){
  
  Predicted <- NULL
  True <- NULL
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # get single benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # get ground truth
  result_true <- ground_truth(full_benchmark)
  result_true <- get_model(result_true, "ground_truth", modality(intermediate_pred(indv_benchmark)), "ground_truth")
  
  # before correction results
  before <- initial_pred(indv_benchmark)
  # after correction results
  after <- intermediate_pred(indv_benchmark)
  
  # comparison between initial signatures and true signatures
  comparison <- initial_comparison(indv_benchmark)
  
  # number of samples
  num_samples <- ncol(exposures(before))
  
  # duplicate groups, as merged by benchmark()
  groups <- .duplicate_groups(comparison)

  final_figures <- list()

  # loop through duplicated signatures
  for (duplicated_sig in names(groups)) {

    before_exposures <- NULL

    # signatures to combine
    sigs_to_merge <- groups[[duplicated_sig]]
    
    # get pre-merge exposures
    for (signature in sigs_to_merge){
      
      #tmp_exposures <- exposures(after)[signature,]
      tmp_exposures <- as.data.frame(as.vector(exposures(before)[signature,]))
      true_exp <- .true_exposure(result_true, duplicated_sig,
                                 colnames(exposures(before)))
      tmp_exposures_all <- cbind(tmp_exposures, as.data.frame(true_exp))
      colnames(tmp_exposures_all) <- c("Predicted", "True")
      tmp_exposures_all$Source <- signature
      
      before_exposures <- rbind(before_exposures, tmp_exposures_all)
    }
    
    rownames(before_exposures) <- c(1: (num_samples * length(sigs_to_merge)))
    
    # plot exposures before correction
    before_plot <- ggplot(before_exposures, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(before_exposures$Predicted), max(before_exposures$True)))) +
      scale_y_continuous(labels = scales::comma,limits = c(0, max(max(before_exposures$Predicted), max(before_exposures$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of duplicated signature, ", duplicated_sig, sep = ""), 
           subtitle = "Before merging",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(before_plot)
    
    # name of the merged signature
    new_sig_name <- .merged_signature_name(sigs_to_merge)
    
    # get post-merge exposures
    tmp_exposures_list2 <- as.data.frame(as.vector(exposures(after)[new_sig_name,]))
    true_exp <- .true_exposure(result_true, duplicated_sig,
                               colnames(exposures(after)))
    after_exposures <- cbind(tmp_exposures_list2, as.data.frame(true_exp))
    colnames(after_exposures) <- c("Predicted", "True")
    
    rownames(after_exposures) <- c(1:num_samples)
    
    # Plot exposures after correction
    after_plot <- ggplot(after_exposures, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(after_exposures$Predicted), max(after_exposures$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(after_exposures$Predicted), max(after_exposures$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of duplicated signature, ", duplicated_sig, sep = ""), 
           subtitle = "After merging",
           x="Predicted Activity (Merged)", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(after_plot)
    
    figure <- ggpubr::ggarrange(before_plot, after_plot, ncol = 2, nrow = 1)
    
    final_figures <- append(final_figures, list(figure))
    
  }
  
  return(final_figures)
  
}


#' @title Plot effect of composite correction
#' @description After a prediction has been benchmarked with the \link{benchmark} function,
#' the true and predicted exposures can be plotted both before and after the
#' composite signature adjustment. The effect of the adjustment can then be
#' observed.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#'
#' @return A list of ggplot objects
#' @examples
#' data(full_benchmark_example)
#' benchmark_plot_composite_exposures(full_benchmark_example, "example_k8")
#' 
#' @export
benchmark_plot_composite_exposures <- function(full_benchmark, method_id){
  
  Predicted <- NULL
  True <- NULL
  
  # check that full_benchmark is a full_benchmark class object
  if (!methods::is(full_benchmark, "full_benchmark")) {
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # get single benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # get ground truth
  result_true <- ground_truth(full_benchmark)
  result_true <- get_model(result_true, "ground_truth", modality(final_pred(indv_benchmark)), "ground_truth")
  
  # pre correction results
  before <- intermediate_pred(indv_benchmark)
  # post correction results
  after <- final_pred(indv_benchmark)
  
  # pre correction comparison between predicted and true signatures
  comparison <- intermediate_comparison(indv_benchmark)
  
  # number of samples
  num_samples <- ncol(exposures(before))
  
  # find composite signatures
  freq <- table(comparison$x_sig_name)
  composite_signatures <- names(freq[freq > 1])
  
  count <- 1
  
  final_figures <- list()
  
  # loop through composite signatures
  for (composite_sig in composite_signatures){
    
    before_exposures <- NULL
    
    # list of SBS signatures that are the components of this composite signature
    sig_components <- comparison[comparison$x_sig_name == composite_sig, 5]
    
    for (component_index in 1:length(sig_components)){
      
      # pre-correction exposures
      tmp_exp <- exposures(before)[composite_sig,]
      tmp_exp_list <- as.data.frame(as.vector(tmp_exp))
      true_exp <- .true_exposure(result_true, sig_components[component_index],
                                 colnames(exposures(before)))
      tmp_exp_all <- cbind(tmp_exp_list, as.data.frame(true_exp))
      colnames(tmp_exp_all) <- c("Predicted", "True")
      tmp_exp_all$Source <- sig_components[component_index]
      
      before_exposures <- rbind(before_exposures, tmp_exp_all)
      
    }
    
    rownames(before_exposures) <- c(1: (num_samples * length(sig_components)))
    
    # plot pre-correction exposures
    before_plot <- ggplot(before_exposures, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(before_exposures$Predicted), max(before_exposures$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(before_exposures$Predicted), max(before_exposures$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of composite signature, ", composite_sig, sep = ""), 
           subtitle = "Before decomposing",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(before_plot)
    
    colnames <- NULL
    for (component in sig_components){
      #colnames <- c(colnames, paste("Signature", component, "_like", sep = ""))
      colnames <- c(colnames, paste(component, "_like", sep = ""))
      
    }
    colnames <- colnames[colnames %in% rownames(exposures(after))]
    sig_components <- gsub("_like", "", colnames)
    
    after_exposures <- NULL
    
    for (component_index in 1:length(sig_components)){
      
      # post-correction exposure
      tmp_exp <- exposures(after)[colnames[component_index],]
      tmp_exp_list <- as.data.frame(as.vector(tmp_exp))
      true_exp <- .true_exposure(result_true, sig_components[component_index],
                                 colnames(exposures(after)))
      tmp_exp_all <- cbind(tmp_exp_list, as.data.frame(true_exp))
      colnames(tmp_exp_all) <- c("Predicted", "True")
      tmp_exp_all$Source <- sig_components[component_index]
      
      after_exposures <- rbind(after_exposures, tmp_exp_all)
      
    }
    
    # plot post-correction exposures
    after_plot <- ggplot(after_exposures, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(after_exposures$Predicted), max(after_exposures$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(after_exposures$Predicted), max(after_exposures$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of composite signature, ", composite_sig, sep = ""), 
           subtitle = "After decomposing",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(after_plot)
    
    figure <- ggpubr::ggarrange(before_plot, after_plot, ncol = 2, nrow = 1)
    
    final_figures <- append(final_figures, list(figure))
    
  }
  
  return(final_figures)
  
}


# The benchmark stage that 'prediction' names: "Initial", "Intermediate" or
# "Final". Also accepts "Init", "Inter" and "Fin", in any case.
.match_benchmark_stage <- function(prediction) {
  stages <- c(
    initial = "Initial", init = "Initial",
    intermediate = "Intermediate", inter = "Intermediate",
    final = "Final", fin = "Final"
  )
  valid <- is.character(prediction) && length(prediction) == 1 &&
    tolower(prediction) %in% names(stages)
  if (!valid) {
    stop(
      "'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'."
    )
  }
  return(unname(stages[tolower(prediction)]))
}

# Duplicate groups in 'comparison': each predicted signature joins the group
# of the true signature it matches best (ties go to the earlier row), and a
# true signature with two or more predicted signatures is duplicated. A list
# named by true signature, in alphabetical order, of the predicted
# signatures to merge, best match first.
.duplicate_groups <- function(comparison) {
  comparison <- comparison[order(comparison$cosine, decreasing = TRUE), ,
                           drop = FALSE]
  best <- comparison[!duplicated(comparison$x_sig_name), , drop = FALSE]
  groups <- split(best$x_sig_name, best$y_sig_name)
  return(groups[lengths(groups) > 1])
}

# Name of the signature merged from the duplicates 'sigs_to_merge'
.merged_signature_name <- function(sigs_to_merge) {
  return(paste0("Merged Signature (", paste(sigs_to_merge, collapse = "."),
                ")"))
}

# Merge each group of duplicates in 'result' (see .duplicate_groups()) into
# one signature, the exposure-weighted mean of the group, whose exposures
# are the sum of the group's
.correct_duplicates <- function(result, comparison) {
  exposures <- exposures(result)
  signatures <- signatures(result)
  groups <- .duplicate_groups(comparison)

  merged_sigs <- NULL
  merged_exposures <- NULL
  for (sigs_to_merge in groups) {
    new_sig_name <- .merged_signature_name(sigs_to_merge)
    group_exposures <- exposures[sigs_to_merge, , drop = FALSE]
    merged <- signatures[, sigs_to_merge, drop = FALSE] %*%
      rowSums(group_exposures)
    merged <- merged / sum(merged)
    colnames(merged) <- new_sig_name
    merged_sigs <- cbind(merged_sigs, merged)
    summed <- matrix(colSums(group_exposures), nrow = 1,
                     dimnames = list(new_sig_name, colnames(exposures)))
    merged_exposures <- rbind(merged_exposures, summed)
  }

  all_merged <- unlist(groups, use.names = FALSE)
  keep_sigs <- !colnames(signatures) %in% all_merged
  keep_rows <- !rownames(exposures) %in% all_merged
  duplicates_corrected <- result
  signatures(duplicates_corrected) <- cbind(
    signatures[, keep_sigs, drop = FALSE], merged_sigs
  )
  exposures(duplicates_corrected) <- rbind(
    exposures[keep_rows, , drop = FALSE], merged_exposures
  )
  return(duplicates_corrected)
}

# Function for addressing composites
.correct_composites <- function(result, compare_cosmic_result, musica, result_true){
  
  Sum <- NULL
  Predicted <- NULL
  True <- NULL
  
  exposures <- exposures(result)
  signatures <- signatures(result)
  num_samples <- dim(exposures)[2]
  
  freq <- table(compare_cosmic_result$x_sig_name)
  composite_signatures <- names(freq[freq > 1])
  #all_sigs_to_merge <- compare_cosmic_result[compare_cosmic_result$y_sig_name == duplicated_signatures, 4]
  
  corrected_sigs <- NULL
  corrected_exposures <- NULL
  
  for (composite_sig in composite_signatures){
    
    composite_exp_all <- NULL
    
    # the full signature to separate
    sig_to_separate <- signatures(result)[, composite_sig]
    exposures_to_separate <- exposures(result)[composite_sig,]
    
    # list of SBS signatures that are the componenets of this composite signature
    sig_components <- compare_cosmic_result[compare_cosmic_result$x_sig_name == composite_sig, 5]
    
    # number of components in this composite sig
    num_components <- length(sig_components)
    
    # data frame of full signatures of components
    component_signatures <- as.data.frame(signatures(result_true)[,sig_components]) # GENERALIZE
    
    separated_sigs <- matrix(ncol = 0, nrow = 96)
    separated_exposures <- matrix(ncol = num_samples, nrow = 0)
    
    # perform nnls   
    nnls_result <- lsei::nnls(as.matrix(component_signatures), as.vector(sig_to_separate))
    
    for (component_index in 1:num_components){
      
      new_signature <- component_signatures[ , component_index] * nnls_result$x[component_index]
      if (sum(new_signature) != 0){
        separated_sigs <- cbind(separated_sigs, new_signature)
      }
      else{
        num_components <- num_components - 1
        sig_components <- sig_components[-component_index]
      }
      
    }
    
    # renormalize
    separated_sigs <- prop.table(separated_sigs,2)
    
    # update exposures
    
    num_samples <- dim(musica@count_tables$SBS96@count_table)[2]
    
    nnls_exposure_results <- data.frame(factor1 = numeric(), factor2 = numeric())
    
    for (sample_index in 1:num_samples){
      
      #A <- as.matrix(signatures(result_true))
      A <- as.matrix(separated_sigs)
      #b <- as.vector(sig_to_separate * exposures_to_separate[sample_index])
      b <- as.vector(musica@count_tables$SBS96@count_table[,sample_index])
      
      nnls_exposure_result <- lsei::nnls(A, b)
      
      nnls_exposure_results[sample_index,1] <- nnls_exposure_result$x[1]
      nnls_exposure_results[sample_index,2] <- nnls_exposure_result$x[2]
      
      
    }
    
    #sig_component_result <- predict_exposure(musica = musica, g = "hg38", table_name = "SBS96", 
    #signature_res = cosmic_v2_sigs, 
    #signatures_to_use =  sig_components, algorithm = "lda")
    
    #tmp_exposures <- exposures(sig_component_result)
    #tmp_exposures <- as.data.frame(t(tmp_exposures))
    
    #tmp_exposures$Sum <- rowSums(tmp_exposures)
    
    # added
    #nnls_exposure_results <- nnls_exposure_results[,c(6,7)]
    #nnls_exposure_results <- nnls_exposure_results[,c(2,1)]
    
    nnls_exposure_results$Sum <- rowSums(nnls_exposure_results)
    
    for (component_index in 1:num_components){
      
      # calculate the percent of the sum that each of the 96 channels contributes (reword this i know it doesnt make sense)
      #tmp_exposures <- transform(tmp_exposures, percent1 = tmp_exposures[component_index] / Sum)
      nnls_exposure_results <- transform(nnls_exposure_results, percent1 = nnls_exposure_results[component_index] / Sum)
      
      #separated_exposures <- rbind(separated_exposures, 
      #exposures_to_separate * tmp_exposures[, component_index + num_components + 1])
      separated_exposures <- rbind(separated_exposures, 
                                   exposures_to_separate * nnls_exposure_results[, component_index + num_components + 1])
      
      
      # for plotting
      tmp_exp <- exposures[composite_sig,]
      tmp_exp_list <- as.data.frame(as.vector(tmp_exp))
      true_exp <- .true_exposure(result_true, sig_components[component_index],
                                 colnames(exposures))
      tmp_exp_all <- cbind(tmp_exp_list, as.data.frame(true_exp)) # GENERALIZE
      colnames(tmp_exp_all) <- c("Predicted", "True")
      tmp_exp_all$Source <- sig_components[component_index]
      
      composite_exp_all <- rbind(composite_exp_all, tmp_exp_all)
      
    }
    
    rownames(composite_exp_all) <- c(1: (num_samples * num_components))
    
    plot <- ggplot(composite_exp_all, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(composite_exp_all$Predicted), max(composite_exp_all$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(composite_exp_all$Predicted), max(composite_exp_all$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of composite signature, ", composite_sig, sep = ""), 
           subtitle = "Before decomposing",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(plot)
    
    colnames <- NULL
    for (component in sig_components){
      #colnames <- c(colnames, paste("Signature", component, "_like", sep = ""))
      colnames <- c(colnames, paste(component, "_like", sep = ""))
      
    }
    
    colnames(separated_sigs) <- colnames
    rownames(separated_exposures) <- colnames
    
    corrected_sigs <- cbind(corrected_sigs, separated_sigs)
    corrected_exposures <- rbind(corrected_exposures, separated_exposures)
    
    # plot exposures after separation
    
    plot_exposures <- composite_exp_all
    plot_exposures$Predicted <- c(t(separated_exposures))
    
    plot <- ggplot(plot_exposures, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(plot_exposures$Predicted), max(plot_exposures$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(plot_exposures$Predicted), max(plot_exposures$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of composite signature, ", composite_sig, sep = ""), 
           subtitle = "After decomposing",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(plot)
    
    
  }
  
  if (is.null(corrected_sigs) == FALSE){
    
    keep_sigs <- !colnames(signatures) %in% composite_signatures
    unchanged_signatures <- as.data.frame(signatures)[, keep_sigs,
                                                      drop = FALSE]
    corrected_sigs <- cbind(unchanged_signatures, corrected_sigs)
    
    keep_rows <- !rownames(exposures) %in% composite_signatures
    unchanged_exposures <- as.data.frame(exposures)[keep_rows, ,
                                                    drop = FALSE]
    corrected_exposures <- rbind(unchanged_exposures, corrected_exposures)
    
  }
  
  else{
    
    keep_sigs <- !colnames(signatures) %in% composite_signatures
    unchanged_signatures <- as.data.frame(signatures)[, keep_sigs,
                                                      drop = FALSE]
    corrected_sigs <- unchanged_signatures
    
    keep_rows <- !rownames(exposures) %in% composite_signatures
    unchanged_exposures <- as.data.frame(exposures)[keep_rows, ,
                                                    drop = FALSE]
    corrected_exposures <- unchanged_exposures
    
  }
  
  composites_corrected <- result
  signatures(composites_corrected) <- as.matrix(corrected_sigs)
  exposures(composites_corrected) <- as.matrix(corrected_exposures)
  
  return(composites_corrected)
  
}

# Matches between the signatures of 'result' and the true signatures of
# 'truth' above 'threshold', without plotting, with low matches of a
# predicted signature that also has a high one removed. No rows if nothing
# matches.
.benchmark_compare <- function(result, truth, threshold,
                               adjustment_threshold) {
  comparison <- .signature_matches(
    signatures(result), signatures(truth), threshold = threshold
  )
  return(.benchmark_comp_adj(comparison, adjustment_threshold))
}

# Function for adjusting comparison 
.benchmark_comp_adj <- function(comparison, adjustment_threshold){
  
  low_threshold_comp <- comparison[comparison$cosine <= adjustment_threshold,]
  high_threshold_comp <- comparison[comparison$cosine > adjustment_threshold,]
  
  indexes_to_keep <- c()
  if (dim(low_threshold_comp)[1] > 0){
    for (index in 1:dim(low_threshold_comp)[1]){
      if (low_threshold_comp[index,4] %in% high_threshold_comp$x_sig_name == FALSE){
        indexes_to_keep <- c(indexes_to_keep, index)
      }
      else{
        existing_cs <- high_threshold_comp[high_threshold_comp$x_sig_name == low_threshold_comp[index,4], 1][1]
        diff <- abs(existing_cs - low_threshold_comp[index,1])
        if (diff < 0.05){
          indexes_to_keep <- c(indexes_to_keep, index)
        }
      }
    }
    
    
    comparison_adj <- rbind(high_threshold_comp, low_threshold_comp[indexes_to_keep,])
    
    return(comparison_adj)
    
  }
  
  else{
    return(comparison)
  }
  
}

# Draw the comparison plot of one stage of a benchmark run, or say there is
# nothing to plot when no signature matched; TRUE if it was drawn
.plot_stage_comparison <- function(full_benchmark, method_id, prediction,
                                   comparison) {
  if (nrow(comparison) == 0) {
    message(
      "No predicted signature matches a true signature at the '",
      prediction, "' stage; skipping its comparison plot."
    )
    return(FALSE)
  }
  benchmark_plot_comparison(full_benchmark, method_id, prediction,
                            same_scale = FALSE)
  return(TRUE)
}

# 'f' (mean, min or max) of cosine similarities 'cs', rounded to 3
# decimals; NA when there are none
.cs_stat <- function(cs, f) {
  if (length(cs) == 0) {
    return(NA_real_)
  }
  return(round(f(cs), 3))
}

# Stop, naming the stage, when 'comparison' has no matches to plot
.check_comparison_has_matches <- function(comparison, prediction) {
  if (nrow(comparison) == 0) {
    stop(
      "No predicted signature matches a true signature at the '",
      prediction, "' stage, so there is nothing to plot."
    )
  }
}

# Function to generate single run summary
.generate_summary <- function(title, result_all, result_true, comparison_results, count_table, final_musica, final_comparison, re_method){
  
  # missing
  
  true_sig_names <- colnames(signatures(result_true))
  num_true <- length(true_sig_names)
  num_missing <- length(true_sig_names[!(true_sig_names %in% comparison_results$y_sig_name)])
  
  # spurious
  
  predicted_sig_names <- colnames(signatures(result_all))
  num_predicted <- length(predicted_sig_names)
  num_spurious <- length(predicted_sig_names[!(predicted_sig_names %in% comparison_results$x_sig_name)])
  
  # duplicate
  
  duplicated_signature_components <- unlist(
    .duplicate_groups(comparison_results), use.names = FALSE
  )
  num_duplicates <- length(duplicated_signature_components)

  # composite
  
  freq <- table(comparison_results$x_sig_name)
  composite_sigs <- names(freq[freq > 1])
  num_composites <- length(composite_sigs)
  
  # dupcomp
  
  dupcomp_count <- 0
  for (sig in composite_sigs) {
    if (sig %in% duplicated_signature_components) {
      dupcomp_count <- dupcomp_count + 1
      num_composites <- num_composites - 1
      num_duplicates <- num_duplicates - 1
    }
  }
  
  # matched
  
  num_direct_matches <- num_predicted - num_composites - num_duplicates - num_spurious - dupcomp_count
  
  # reconstruction error
  
  initial_reconstruction_error <- .get_reconstruction_error(result_all, count_table, re_method)
  initial_reconstruction_error <- round(initial_reconstruction_error, 3)
  end_reconstruction_error <- .get_reconstruction_error(final_musica, count_table, re_method)
  end_reconstruction_error <- round(end_reconstruction_error, 3)
  
  # stability
  
  #  average cosine similarity
  
  cs <- final_comparison$cosine
  avg_cs <- .cs_stat(cs, mean)

  # min and max cosine similarity
  min_cs <- .cs_stat(cs, min)
  max_cs <- .cs_stat(cs, max)
  
  # signatures found
  
  sigs_found <- final_comparison$y_sig_name
  num_found <- length(sigs_found)
  sigs_found <- paste(sigs_found, collapse = ", ") # added when removed below
  #sigs_found_ordered <- paste("SBS", sort(as.numeric(str_remove(sigs_found, "Signature"))), sep = "")
  #sigs_found_ordered <- paste(sigs_found_ordered, collapse = ", ")
  
  
  # one row, named by the method ID, with typed columns
  summary <- data.frame(
    "Num Found" = as.integer(num_found),
    "Num Direct Matched" = as.integer(num_direct_matches),
    "Num Missing" = as.integer(num_missing),
    "Num Spurious" = as.integer(num_spurious),
    "Num Duplicates" = as.integer(num_duplicates),
    "Num Composites" = as.integer(num_composites),
    "Num Dup/Comp" = as.integer(dupcomp_count),
    "Initial RE" = initial_reconstruction_error,
    "Final RE" = end_reconstruction_error,
    "Mean CS" = avg_cs,
    "Min CS" = min_cs,
    "Max CS" = max_cs,
    "Sigs Found" = sigs_found,
    row.names = title, check.names = FALSE, stringsAsFactors = FALSE
  )

  return(summary)
  
}

# One row per true signature (and any other signature found), summarizing
# the Final comparisons of all runs in 'benchmark'
.signature_view_summary <- function(benchmark, prediction_modality) {
  runs <- indv_benchmarks(benchmark)
  truth <- get_model(ground_truth(benchmark), "ground_truth",
                     prediction_modality, "ground_truth")
  comparisons <- lapply(runs, final_comparison)
  methods <- vapply(runs, method_id, character(1))

  found_sigs <- unlist(lapply(comparisons, function(x) x$y_sig_name))
  all_sigs <- unique(c(found_sigs, colnames(signatures(truth))))

  rows <- lapply(all_sigs, function(sig) {
    found <- vapply(comparisons, function(x) sig %in% x$y_sig_name,
                    logical(1))
    cs <- unlist(lapply(comparisons, function(x) {
      x$cosine[x$y_sig_name == sig]
    }))
    data.frame(
      "Times Found" = sum(found),
      "Times Missed" = sum(!found),
      "Mean CS" = .cs_stat(cs, mean),
      "Min CS" = .cs_stat(cs, min),
      "Max CS" = .cs_stat(cs, max),
      "Methods Found" = paste(methods[found], collapse = ", "),
      row.names = sig, check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  return(do.call(rbind, rows))
}

# Function for updating full benchmark object
.update_benchmark <- function(full_benchmark, indv_benchmark, single_summary, prediction_modality){
  
  # update summary, one row per run
  if (nrow(method_view_summary(full_benchmark)) == 0) {
    method_view_summary(full_benchmark) <- single_summary
  } else {
    method_view_summary(full_benchmark) <- rbind(
      method_view_summary(full_benchmark), single_summary
    )
  }
  
  # update single benchmark list
  #NOTE: need to get those setters working
  indv_benchmarks(full_benchmark) <- append(indv_benchmarks(full_benchmark), list(indv_benchmark))
  names(indv_benchmarks(full_benchmark))[length(indv_benchmarks(full_benchmark))] <- method_id(indv_benchmark)
  
  # update sig view summary
  sig_view_summary(full_benchmark) <- .signature_view_summary(
    full_benchmark, prediction_modality
  )
  
  return(full_benchmark)
  
}
