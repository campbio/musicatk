#' @importFrom lsei nnls

#' @title Create a full_benchmark object 
#' @description Initialize a \code{\linkS4class{full_benchmark}} object for benchmarking
#'
#' @param true_signatures A matrix of true signatures by mutational motifs
#' @param true_exposures A matrix of samples by true signature weights
#' @param count_table Summary table with per-sample un-normalized motif counts
#' @param variant_class Mutations are SBS, DBS, or Indel.
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
  
  musica_truth <- create_musica_from_counts(count_table, variant_class)
  
  add_result(as.matrix(true_signatures), as.matrix(true_exposures), musica_truth, "ground_truth",
             variant_class, "ground_truth")
  
  # create benchmark object
  full_benchmark <- new("full_benchmark", ground_truth = musica_truth)
  
  return(full_benchmark)
  
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
#' @param prediction_modality The modality of the prediction to benchmark. Must be
#' "SBS96", "DBS78", or "IND83".
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
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object. ", 
         "Use function 'create_benchmark' to initialize a 'full_benchmark' object")
  }
  
  # check that prediction is a result model object
  if (class(prediction_musica)[1] != "musica"){
    stop("'prediction' must be a 'musica' object.")
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
  
  # extract ground truth musica
  truth <- ground_truth(full_benchmark)
  
  # temporary musica object for benchmarking
  temp_musica <- prediction_musica
  
  # inital comparison of predicted and true signatures
  message("\nComparing to true signatures (initial)...")
  initial_comparison <- compare_results(musica = temp_musica, model_id = prediction_model_id, 
                                        other_model_id = "ground_truth",
                                        modality = prediction_modality, result_name = prediction_result_name, 
                                        other_musica = truth, other_result_name = "ground_truth",
                                        threshold = threshold, result_rename =
                                          paste(prediction_model_id, "- Initial"),
                                        other_result_rename =
                                          "Ground Truth")
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  initial_comparison <- .benchmark_comp_adj(initial_comparison, adjustment_threshold)
  
  # correct duplicate signatures
  message("Correcting duplicates...")
  duplicates_corrected <- .correct_duplicates(get_model(prediction_musica, prediction_result_name,
                                                        prediction_modality, prediction_model_id),
                                              initial_comparison,
                                              get_model(truth, "ground_truth", prediction_modality, "ground_truth"))
  
  # add result_model with duplicates corrected to temp_musica
  add_result(signatures(duplicates_corrected), exposures(duplicates_corrected), temp_musica, prediction_result_name,
             prediction_modality, "duplicates_corrected")
  
  # intermediate comparison of predicted and true signatures
  message("Comparing to true signatures (post duplicates corrected)...")
  duplicates_corrected_comparison <- compare_results(musica = temp_musica, model_id = "duplicates_corrected",
                                                     other_model_id = "ground_truth",
                                                     modality = prediction_modality, result_name = prediction_result_name,
                                                     other_musica = truth, other_result_name = "ground_truth",
                                                     threshold = threshold, result_rename =
                                                       paste(prediction_model_id, "- Intermediate"),
                                                     other_result_rename =
                                                       "Ground Truth")
  
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  duplicates_corrected_comparison <- .benchmark_comp_adj(duplicates_corrected_comparison, adjustment_threshold)
  
  # correct composite signatures
  message("Correcting composites...")
  composites_corrected <- .correct_composites(get_model(temp_musica, prediction_result_name,
                                                        prediction_modality, "duplicates_corrected"),
                                              duplicates_corrected_comparison,
                                              musica = temp_musica,
                                              get_model(truth, "ground_truth", prediction_modality, "ground_truth"))
  
  # add result_model with composites corrected to temp_musica
  add_result(signatures(composites_corrected), exposures(composites_corrected), temp_musica, prediction_result_name,
             prediction_modality, "composites_corrected")
  
  # final comparison between predicted and true signatures
  message("Comparing to true signatures (post composites corrected)...")
  composites_corrected_comparison <- compare_results(musica = temp_musica, model_id = "composites_corrected",
                                                     other_model_id = "ground_truth",
                                                     modality = prediction_modality, result_name = prediction_result_name,
                                                     other_musica = truth, other_result_name = "ground_truth",
                                                     threshold = threshold, result_rename =
                                                       paste(prediction_model_id, "- Final"),
                                                     other_result_rename =
                                                       "Ground Truth")
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  composites_corrected_comparison <- .benchmark_comp_adj(composites_corrected_comparison, adjustment_threshold)
  
  # extract count table
  count_table <- extract_count_tables(temp_musica)
  count_table <- get_count_table(count_table[[prediction_modality]])
  
  # create summary matrix
  message("Creating summary...")
  single_summary <- as.matrix(.generate_summary(method_id, 
                                                get_model(prediction_musica, prediction_result_name,
                                                                     prediction_modality, prediction_model_id), 
                                                get_model(truth, "ground_truth", prediction_modality, "ground_truth"), 
                                                initial_comparison, 
                                                count_table, 
                                                composites_corrected, 
                                                composites_corrected_comparison, re_method))
  
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
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Initial", same_scale = FALSE)

    # Duplicate signature exposures before/afters
    duplicate_plot <- benchmark_plot_duplicate_exposures(full_benchmark, method_id)
    print(duplicate_plot)
    
    # Intermediate signatures plot
    intermediate_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Intermediate")
    print(intermediate_sig_plot)
    # Intermediate comparison plot
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Intermediate", same_scale = FALSE)
    
    # Composite signature exposures before/afters
    composite_plot <- benchmark_plot_composite_exposures(full_benchmark, method_id)
    print(composite_plot)
    
    # Final signatures plot
    final_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Final")
    print(final_sig_plot)
    # Final comparison plot
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Final", same_scale = FALSE)
    # Final exposures plot
    exposure_plot <- benchmark_plot_exposures(full_benchmark, method_id, prediction = "Final")
    print(exposure_plot)
    
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
  if (class(full_benchmark)[1] != "full_benchmark"){
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
  
  # check if prediction is one of Initial, Intermediate, or Final
  valid <- c("Initial", "initial", "Init", "init", "Intermediate", "intermediate",
             "Inter", "inter", "Final", "final", "Fin", "fin")
  if (!(prediction %in% valid)){
    stop("'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'.")
  }
  
  # access musica object for desired prediction
  if (prediction %in% c("Initial", "initial", "Init", "init")){
    result <- initial_pred(indv_benchmark)
  }
  else if (prediction %in% c("Intermediate", "intermediate", "Inter", "inter")){
    result <- intermediate_pred(indv_benchmark)
  }
  else if (prediction %in% c("Final", "final", "Fin", "fin")){
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
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # check if prediction is one of Initial, Intermediate, or Final
  valid <- c("Initial", "initial", "Init", "init", "Intermediate", "intermediate",
             "Inter", "inter", "Final", "final", "Fin", "fin")
  if (!(prediction %in% valid)){
    stop("'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'.")
  }
  
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
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # check if prediction is one of Initial, Intermediate, or Final
  valid <- c("Initial", "initial", "Init", "init", "Intermediate", "intermediate",
             "Inter", "inter", "Final", "final", "Fin", "fin")
  if (!(prediction %in% valid)){
    stop("'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'.")
  }
  
  # access individual benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # access desired prediction
  if (prediction == "Initial"){
    comparison <- initial_comparison(indv_benchmark)
  }
  else if (prediction == "Intermediate"){
    comparison <- intermediate_comparison(indv_benchmark)
  }
  else if (prediction == "Final"){
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
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # check if prediction is one of Initial, Intermediate, or Final
  valid <- c("Initial", "initial", "Init", "init", "Intermediate", "intermediate",
             "Inter", "inter", "Final", "final", "Fin", "fin")
  if (!(prediction %in% valid)){
    stop("'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'.")
  }
  
  comparison <- benchmark_get_comparison(full_benchmark, method_id, prediction)
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
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object.")
  }
  
  # check if this method_id exists
  if (!(method_id %in% names(indv_benchmarks(full_benchmark)))){
    stop("'method_id' ", deparse(substitute(method_id)), " not found in ", deparse(substitute(full_benchmark)))
  }
  
  # check if prediction is one of Initial, Intermediate, or Final
  valid <- c("Initial", "initial", "Init", "init", "Intermediate", "intermediate",
             "Inter", "inter", "Final", "final", "Fin", "fin")
  if (!(prediction %in% valid)){
    stop("'prediction' must be one of: 'Initial', 'Intermediate', or 'Final'.")
  }
  
  # access individual benchmark object
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # access comparison for desired prediction
  comparison <- benchmark_get_comparison(full_benchmark, method_id, prediction)
  
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
    true <- c(true, exposures(truth)[,true_sig])
    sig <- c(sig, rep(true_sig, dim(exposures(truth))[1]))
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
  if (class(full_benchmark)[1] != "full_benchmark"){
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
  num_samples <- dim(exposures(result_true))[1]
  
  # determine duplicate signatures
  freq <- table(comparison$y_sig_name)
  duplicated_signatures <- names(freq[freq > 1])
  
  final_figures <- list()
  
  # loop through duplicated signatures
  for (duplicated_sig in duplicated_signatures){
    
    before_exposures <- NULL
    
    # signatures to combine
    sigs_to_merge <- comparison[comparison$y_sig_name == duplicated_sig & !grepl("like", comparison$x_sig_name), 4]
    
    # get pre-merge exposures
    for (signature in sigs_to_merge){
      
      #tmp_exposures <- exposures(after)[signature,]
      tmp_exposures <- as.data.frame(as.vector(exposures(before)[signature,]))
      tmp_exposures_all <- cbind(tmp_exposures, 
                                 as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
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
    
    # create new name for merged signature
    new_sig_name <- "Merged Signature ("
    for (sig in sigs_to_merge){
      if (sig == sigs_to_merge[1]){
        new_sig_name <- paste(new_sig_name, sig, sep = "")
      }
      else{
        new_sig_name <- paste(new_sig_name, ".", sig, sep = "")
      }
    }
    new_sig_name <- paste(new_sig_name, ")", sep = "")
    
    # get post-merge exposures
    tmp_exposures_list2 <- as.data.frame(as.vector(exposures(after)[new_sig_name,]))
    after_exposures <- cbind(tmp_exposures_list2, as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
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
  if (class(full_benchmark)[1] != "full_benchmark"){
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
  num_samples <- dim(exposures(result_true))[1]
  
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
      tmp_exp_all <- cbind(tmp_exp_list, 
                           as.data.frame(as.numeric(exposures(result_true)[,sig_components[component_index]])))
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
      tmp_exp_all <- cbind(tmp_exp_list, 
                           as.data.frame(as.numeric(exposures(result_true)[,sig_components[component_index]])))
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


# Function for addressing duplicates
.correct_duplicates <- function(result, compare_cosmic_result, result_true){
  
  Sum <- NULL
  Predicted <- NULL
  True <- NULL
  
  exposures <- exposures(result)
  signatures <- signatures(result)
  num_samples <- dim(exposures)[2]
  
  freq <- table(compare_cosmic_result$y_sig_name)
  duplicated_signatures <- names(freq[freq > 1])
  all_sigs_to_merge <- compare_cosmic_result[compare_cosmic_result$y_sig_name %in% duplicated_signatures & !grepl("like", compare_cosmic_result$x_sig_name), 4]
  
  corrected_sigs <- NULL
  corrected_exposures <- NULL
  
  for (duplicated_sig in duplicated_signatures){
    
    duplicate_exposures_all <- NULL
    
    sigs_to_merge <- compare_cosmic_result[compare_cosmic_result$y_sig_name == duplicated_sig & !grepl("like", compare_cosmic_result$x_sig_name), 4]
    
    # add when removed the below
    full_sig_names_to_merge <- sigs_to_merge
    #sigs_to_merge <- str_remove(sigs_to_merge, "Signature")
    
    # convert sig numbers to full names
    #full_sig_names_to_merge <- NULL
    #for (sig_number in sigs_to_merge){
    #  sig_name <- paste("Signature", sig_number, sep = "")
    #  full_sig_names_to_merge <- c(full_sig_names_to_merge, sig_name)
    #}
    
    sum <- 0
    for (signature in full_sig_names_to_merge){
      for (sample_index in 1:num_samples){
        temp <- signatures[, signature] * exposures[signature, sample_index]
        sum <- sum + temp
      }
      
      tmp_exposures <- exposures[signature,]
      tmp_exposures_list <- as.data.frame(as.vector(tmp_exposures))
      tmp_exposures_all <- cbind(tmp_exposures_list, 
                                 as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
      colnames(tmp_exposures_all) <- c("Predicted", "True")
      tmp_exposures_all$Source <- signature
      
      duplicate_exposures_all <- rbind(duplicate_exposures_all, tmp_exposures_all)
    }
    
    ## DOESNT WORK WHEN NOT ACTUALLY DUPICATE
    rownames(duplicate_exposures_all) <- c(1: (num_samples * length(sigs_to_merge)))
    
    plot <- ggplot(duplicate_exposures_all, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      facet_wrap(~Source, scales = "fixed") +
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(duplicate_exposures_all$Predicted), max(duplicate_exposures_all$True)))) +
      scale_y_continuous(labels = scales::comma,limits = c(0, max(max(duplicate_exposures_all$Predicted), max(duplicate_exposures_all$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of duplicated signature, ", duplicated_sig, sep = ""), 
           subtitle = "Before merging",
           x="Predicted Activity", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(plot)
    
    
    total_sum <- sum(sum)
    
    normalized <- sum / total_sum
    
    merged_signature <- as.matrix(normalized)
    
    new_sig_name <- "Merged Signature ("
    for (sig in sigs_to_merge){
      if (sig == sigs_to_merge[1]){
        new_sig_name <- paste(new_sig_name, sig, sep = "")
      }
      else{
        new_sig_name <- paste(new_sig_name, ".", sig, sep = "")
      }
    }
    new_sig_name <- paste(new_sig_name, ")", sep = "")
    
    colnames(merged_signature) <- new_sig_name
    
    corrected_sigs <- cbind(corrected_sigs, merged_signature)
    
    # update exposures
    
    merged_exposures <- 0
    for (signature in full_sig_names_to_merge){
      merged_exposures <- merged_exposures + exposures[signature,]
    }
    
    merged_exposures <- t(as.matrix(merged_exposures))
    
    rownames(merged_exposures) <- new_sig_name
    
    corrected_exposures <- rbind(corrected_exposures, merged_exposures)
    
    # exposure plot after merging
    
    tmp_exposures_list <- as.data.frame(as.vector(merged_exposures))
    tmp_exposures_all <- cbind(tmp_exposures_list, as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
    colnames(tmp_exposures_all) <- c("Predicted", "True")
    
    rownames(tmp_exposures_all) <- c(1:num_samples)
    
    plot <- ggplot(tmp_exposures_all, aes(x = Predicted, y = True)) + 
      geom_point(size = 3) + 
      scale_x_continuous(labels = scales::comma, limits = c(0, max(max(tmp_exposures_all$Predicted), max(tmp_exposures_all$True)))) +
      scale_y_continuous(labels = scales::comma, limits = c(0, max(max(tmp_exposures_all$Predicted), max(tmp_exposures_all$True)))) +
      theme_classic() + 
      labs(title=paste("Exposures of duplicated signature, ", duplicated_sig, sep = ""), 
           subtitle = "After merging",
           x="Predicted Activity (Merged)", y = "True Activity") + 
      theme(text = element_text(size=15), axis.text = element_text(size = 15)) + 
      geom_abline() + 
      geom_smooth(method = "lm") +
      theme(legend.title=element_blank())
    
    #print(plot)
    
  }
  
  if (is.null(corrected_sigs) == FALSE){
    
    unchanged_signatures <- as.data.frame(signatures)[, !colnames(signatures) %in% all_sigs_to_merge]
    corrected_sigs <- cbind(unchanged_signatures, corrected_sigs)
    
    unchanged_exposures <- as.data.frame(exposures)[!rownames(exposures) %in% all_sigs_to_merge,]
    corrected_exposures <- rbind(unchanged_exposures, corrected_exposures)
    
  }
  
  else{
    
    unchanged_signatures <- as.data.frame(signatures)[, !colnames(signatures) %in% all_sigs_to_merge]
    corrected_sigs <- unchanged_signatures
    
    unchanged_exposures <- as.data.frame(exposures)[!rownames(exposures) %in% all_sigs_to_merge,]
    corrected_exposures <- unchanged_exposures
    
  }
  
  duplicates_corrected <- result
  signatures(duplicates_corrected) <- as.matrix(corrected_sigs)
  exposures(duplicates_corrected) <- as.matrix(corrected_exposures)
  
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
      tmp_exp_all <- cbind(tmp_exp_list, 
                           as.data.frame(as.numeric(exposures(result_true)[,sig_components[component_index]]))) # GENERALIZE (was                                                                                                                      Signature not SBS)
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
    
    unchanged_signatures <- as.data.frame(signatures)[, !colnames(signatures) %in% composite_signatures]
    corrected_sigs <- cbind(unchanged_signatures, corrected_sigs)
    
    unchanged_exposures <- as.data.frame(exposures)[!rownames(exposures) %in% composite_signatures,]
    corrected_exposures <- rbind(unchanged_exposures, corrected_exposures)
    
  }
  
  else{
    
    unchanged_signatures <- as.data.frame(signatures)[, !colnames(signatures) %in% composite_signatures]
    corrected_sigs <- unchanged_signatures
    
    unchanged_exposures <- as.data.frame(exposures)[!rownames(exposures) %in% composite_signatures,]
    corrected_exposures <- unchanged_exposures
    
  }
  
  composites_corrected <- result
  signatures(composites_corrected) <- as.matrix(corrected_sigs)
  exposures(composites_corrected) <- as.matrix(corrected_exposures)
  
  return(composites_corrected)
  
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

# Function for claculating RE
.get_reconstruction_error <- function(result, count_table){
  
  # extract exposures and signatures matrices
  expos <- exposures(result)
  sigs <- signatures(result)
  
  # convert count table to probabilities
  count_table_probs <- prop.table(count_table, 2)
  
  # convert exposure matrix to probabilities
  expos_probs <- prop.table(expos, 2)
  
  # multiply signature and exposure matrices
  predicted_count_probs <- sigs %*% expos_probs
  
  # calculate sum of differences (reconstruction error)
  reconstruction_error <- sum(abs(count_table_probs - predicted_count_probs))
  
  # return reconstruction error
  return(reconstruction_error)
  
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
  
  num_duplicates <- 0
  
  freq <- table(comparison_results$y_sig_name)
  
  if(length(names(freq[freq > 1])) != 0){
    duplicated_signatures <- names(freq[freq > 1])
    duplicated_signature_components <- comparison_results[comparison_results$y_sig_name %in% duplicated_signatures, 4]
    
    for (duplicated_sig in duplicated_signatures){
      sigs_to_merge <- comparison_results[comparison_results$y_sig_name == duplicated_sig 
                                          & !grepl("like", comparison_results$x_sig_name), 4]
      
      num_duplicates <- num_duplicates + length(sigs_to_merge)
    }
  }
  
  # composite
  
  freq <- table(comparison_results$x_sig_name)
  composite_sigs <- names(freq[freq > 1])
  num_composites <- length(composite_sigs)
  
  # dupcomp
  
  dupcomp_count <- 0
  if (exists("duplicated_signature_components")){
    for (sig in composite_sigs){
      if (sig %in% duplicated_signature_components){
        dupcomp_count <- dupcomp_count + 1
        num_composites <- num_composites - 1
        num_duplicates <- num_duplicates - 1
      }
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
  avg_cs <- mean(cs)
  avg_cs <- round(avg_cs, 3)
  
  # min cosine similarity
  
  min_cs <- min(cs)
  min_cs <- round(min_cs, 3)
  
  # max cosine similarity
  
  max_cs <- max(cs)
  max_cs <- round(max_cs, 3)
  
  # signatures found
  
  sigs_found <- final_comparison$y_sig_name
  num_found <- length(sigs_found)
  sigs_found <- paste(sigs_found, collapse = ", ") # added when removed below
  #sigs_found_ordered <- paste("SBS", sort(as.numeric(str_remove(sigs_found, "Signature"))), sep = "")
  #sigs_found_ordered <- paste(sigs_found_ordered, collapse = ", ")
  
  
  # put in table
  
  summary <- data.frame(matrix(, nrow=13, ncol=1))
  
  rownames(summary) <- c("Num Found", "Num Direct Matched", "Num Missing", "Num Spurious", "Num Duplicates", "Num Composites", "Num Dup/Comp", "Initial RE", "Final RE", "Mean CS", "Min CS", "Max CS", "Sigs Found") 
  
  colnames(summary) <- title
  
  summary[1,1] <- num_found
  summary[2,1] <- num_direct_matches
  summary[3,1] <- num_missing
  summary[4,1] <- num_spurious
  summary[5,1] <- num_duplicates
  summary[6,1] <- num_composites
  summary[7,1] <- dupcomp_count
  summary[8,1] <- initial_reconstruction_error
  summary[9,1] <- end_reconstruction_error
  summary[10,1] <- avg_cs
  summary[11,1] <- min_cs
  summary[12,1] <- max_cs
  summary[13,1] <- sigs_found
  
  
  return(summary)
  
}

# Function to make sig view summary
.signature_view_summary <- function(benchmark, prediction_modality){
  
  indv_benchmarks <- indv_benchmarks(benchmark)
  truth <- ground_truth(benchmark)
  
  truth <- get_model(truth, "ground_truth", prediction_modality, "ground_truth")
  
  final_comparison_list <- list()
  method_list <- c()
  
  for (index in 1:length(indv_benchmarks)){
    
    final_comparison_list[[index]] <- final_comparison(indv_benchmarks[[index]])
    method_list <- c(method_list, method_id(indv_benchmarks[[index]]))
    
  }
  
  summary_complete <- data.frame(matrix(, nrow=6, ncol=1))
  
  sigs_found <- NULL
  for (final_comparison in final_comparison_list){
    sigs_found <- c(sigs_found, final_comparison$y_sig_name)
  }
  
  sigs_found <- unique(sigs_found)
  sigs_found_ordered <- sigs_found
  #sigs_found_ordered <- paste("Signature", sort(as.numeric(str_remove(sigs_found, "Signature"))), sep = "")
  
  all_sigs <- unique(c(sigs_found_ordered, colnames(signatures(truth))))
  all_sigs_ordered <- all_sigs
  #all_sigs_ordered <- paste("Signature", sort(as.numeric(str_remove(all_sigs, "Signature"))), sep = "")
  
  for (sig in all_sigs_ordered){
    
    summary <- data.frame(matrix(, nrow=6, ncol=1))
    rownames(summary) <- c("Times Found", "Times Missed", "Mean CS", "Min CS", "Max CS", "Methods Found") 
    colnames(summary) <- sig
    
    # times found
    
    count_found <- 0
    index_found <- NULL
    index <- 1
    for (final_comparison in final_comparison_list){
      if (sig %in% final_comparison$y_sig_name == TRUE){
        count_found <- count_found + 1
        index_found <- c(index_found, index)
      }
      index <- index + 1
    }
    
    # times missed
    
    count_missed <- 0
    for (final_comparison in final_comparison_list){
      if (sig %in% final_comparison$y_sig_name == FALSE){
        count_missed <- count_missed + 1
      }
    }
    
    # mean CS
    
    cs <- NULL
    for (final_comparison in final_comparison_list){
      if (sig %in% final_comparison$y_sig_name == TRUE){
        cs <- c(cs, final_comparison[final_comparison$y_sig_name == sig, 1])
      }
    }
    
    mean_cs <- mean(cs)
    mean_cs <- round(mean_cs, 3)
    if (is.null(cs)){
      mean_cs <- NA
    }
    
    # min cs
    
    min_cs <- min(cs)
    min_cs <- round(min_cs, 3)
    if (is.null(cs)){
      min_cs <- NA
    }
    
    # max cs
    
    max_cs <- max(cs)
    max_cs <- round(max_cs, 3)
    if (is.null(cs)){
      max_cs <- NA
    }
    
    # methods
    
    methods <- method_list[c(index_found)]
    methods <- paste(methods, collapse = ", ")
    
    summary[1,1] <- count_found
    summary[2,1] <- count_missed
    summary[3,1] <- mean_cs
    summary[4,1] <- min_cs
    summary[5,1] <- max_cs
    summary[6,1] <- methods
    
    summary_complete <- cbind(summary_complete, summary)
    
  }
  
  summary_complete <- summary_complete[,-1]
  summary_complete[sapply(summary_complete, is.infinite)] <- NA
  
  #summary_complete <- t(summary_complete)
  
  #print(summary_complete)
  return(summary_complete)
  
}

# Function for updating full benchmark object
.update_benchmark <- function(full_benchmark, indv_benchmark, single_summary, prediction_modality){
  
  # update summary
  if (dim(method_view_summary(full_benchmark))[1] == 0 ){
    method_view_summary(full_benchmark) <- single_summary
  }else{
    method_view_summary(full_benchmark) <- cbind(method_view_summary(full_benchmark), single_summary)
  }
  
  # update single benchmark list
  #NOTE: need to get those setters working
  indv_benchmarks(full_benchmark) <- append(indv_benchmarks(full_benchmark), list(indv_benchmark))
  names(indv_benchmarks(full_benchmark))[length(indv_benchmarks(full_benchmark))] <- method_id(indv_benchmark)
  
  # update sig view summary
  sig_view_summary(full_benchmark) <- as.matrix(.signature_view_summary(full_benchmark, prediction_modality))
  
  return(full_benchmark)
  
}
