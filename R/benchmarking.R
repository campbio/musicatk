#' Create a full_benchmark object 
#' 
#' Initialize a \code{\linkS4class{full_benchmark}} object for benchmarking
#'
#' @param true_signatures A matrix of true signatures by mutational motifs
#' @param true_exposures A matrix of samples by true signature weights
#' @param count_table Summary table with per-sample unnormalized motif counts
#'
#' @return A \code{\linkS4class{full_benchmark}} object
#' @export
create_benchmark <- function(true_signatures, true_exposures, count_table){
  
  # create musica result object to hold true exposures and signatures
  truth <- create_musica_result(true_signatures, true_exposures, count_table)
  
  # create benchmark object
  full_benchmark <- new("full_benchmark", ground_truth = truth)
  
  return(full_benchmark)
  
}


#' Run the benchmark framework on a prediction
#' 
#' Perform benchmarking on a signature discovery prediction compared
#' to a ground truth. Potential errors in the predicted signatures, such as
#' composite or duplicate signatures are adjusted, and a summary of the accuracy
#' of the prediction is given.
#'
#' @param full_benchmark An object of class \code{\linkS4class{full_benchmark}}
#' created with the \link{create_benchmark} function or returned from a previous
#' \code{benchmark} run.
#' @param prediction An object of class \code{\linkS4class{musica_result}}
#' containing the predicted signatures and exposures to benchmark.
#' @param method_id An identifier for the prediction being benchmarked. If not
#' supplied, it will be automatically set to the variable name of the prediction
#' provided. Default \code{NULL}.
#' @param threshold Cosine similarity cutoff for comparing preidcted and true
#' signatures. Default \code{0.8}.
#' @param adjustment_threshold Cosine similarity value of high confidence.
#' Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold. Default \code{0.9}.
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
#' @export
benchmark <- function(full_benchmark, prediction, method_id = NULL, 
                      threshold = 0.8, adjustment_threshold = 0.9, description = NULL,
                      plot = TRUE, make_copy = FALSE){
  

  if (make_copy == FALSE){
    var_name <- deparse(substitute(full_benchmark))
  }
  
  # check that full_benchmark is a full_benchmark class object
  if (class(full_benchmark)[1] != "full_benchmark"){
    stop(deparse(substitute(full_benchmark)), " is not a 'full_benchmark' object. ", 
         "Use function 'create_benchmark' to initialize a 'full_benchmark' object")
  }
  
  # check that prediction is a musica result object
  if (class(prediction)[1] != "musica_result"){
    stop("'prediction' must be a 'musica_result' object.")
  }
  
  # if no ID provided, try to make one automatically
  if (is.null(method_id)){
    
    # create ID
    method_id <- deparse(substitute(prediction))
    
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
    message("method_id ", original_id, " already exists. method_id updated to ", method_id)
    
  }
  
  
  if (threshold != 0.8){
    warning("Default threshold overriden. Interpret results with caution if 
            comparing benchmark runs with inconsistent thresholds.")
  }
  
  if (adjustment_threshold != 0.9){
    warning("Default adjustment threshold overriden. Interpret results with caution if 
            comparing benchmark runs with inconsistent thresholds.")
  }
  
  full_benchmark@ground_truth@musica <- get_musica(prediction)
  
  truth <- ground_truth(full_benchmark)
  
  message("\nComparing to true signatures (initial)...")
  
  initial_comparison <- compare_results(prediction, truth, threshold)
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  initial_comparison <- .benchmark_comp_adj(initial_comparison, adjustment_threshold)
  
  message("Correcting duplicates...")
  
  duplicates_corrected <- .correct_duplicates(prediction, initial_comparison, truth)
  
  message("Comparing to true signatures (post duplicates corrected)...")
  
  duplicates_corrected_comparison <- compare_results(duplicates_corrected, truth, 
                                                     threshold = threshold)
  
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  duplicates_corrected_comparison <- .benchmark_comp_adj(duplicates_corrected_comparison, adjustment_threshold)
  
  message("Correcting composites...")
  
  composites_corrected <- .correct_composites(duplicates_corrected, duplicates_corrected_comparison, get_musica(truth), truth)
  
  message("Comparing to true signatures (post composites corrected)...")
  
  composites_corrected_comparison <- compare_results(composites_corrected, truth, threshold = threshold)
  
  # remove anything below 0.9 that already appears with at least 0.05 difference
  composites_corrected_comparison <- .benchmark_comp_adj(composites_corrected_comparison, adjustment_threshold)
  
  # extract count table
  count_table <- extract_count_tables(get_musica(prediction))
  count_table <- count_table$SBS96@count_table
  
  message("Creating summary...")
  
  single_summary <- as.matrix(.generate_summary(method_id, prediction, truth, initial_comparison, 
                                         count_table, composites_corrected, composites_corrected_comparison))
  
  message("Creating individual benchmark object...")
  
  # create single benchmark object
  indv_benchmark <- new("single_benchmark", initial_pred = prediction, intermediate_pred = duplicates_corrected, 
                        final_pred = composites_corrected, initial_comparison = initial_comparison,
                        intermediate_comparison = duplicates_corrected_comparison,
                        final_comparison = composites_corrected_comparison, 
                        single_summary = single_summary, method_id = method_id, threshold = threshold,
                        adjustment_threshold = adjustment_threshold, description = description)
  
  # update full benchmark object
  message("Updating full benchmark object...")
  full_benchmark <- .update_benchmark(full_benchmark, indv_benchmark, single_summary)
  
  # update global variable
  if (make_copy == FALSE){
    assign(var_name, full_benchmark, envir = parent.frame())
  }
  
  if (plot == TRUE){
    
    # plots
    message("Generating plots...")
    
    #message("\nInitial Signatures...\n")
    initial_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Initial", same_scale = FALSE)
    print(initial_sig_plot)
    #message("\nInitial Comparison...\n")
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Initial", same_scale = FALSE)
    #message("\nInitial Exposure Comparison...\n")
    #benchmark_plot_exposures(full_benchmark, method_id, prediction = "Initial")
    
    #message("\nDuplicate signature exposures before/afters...\n")
    duplicate_plot <- benchmark_plot_duplicate_exposures(full_benchmark, method_id)
    print(duplicate_plot)
    
    #message("\nIntermediate Signatures...\n")
    intermediate_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Intermediate", same_scale = FALSE)
    print(intermediate_sig_plot)
    #message("\nIntermediate Comparison...\n")
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Intermediate", same_scale = FALSE)
    #message("\nIntermediate Exposure Comparison...\n")
    #benchmark_plot_exposures(full_benchmark, method_id, prediction = "Intermediate")
    
    #message("\nComposite signature exposures before/afters...\n")
    composite_plot <- benchmark_plot_composite_exposures(full_benchmark, method_id)
    print(composite_plot)
    
    #message("\nFinal Signatures...\n")
    final_sig_plot <- benchmark_plot_signatures(full_benchmark, method_id, prediction = "Final", same_scale = FALSE)
    print(final_sig_plot)
    #message("\nFinal Comparison...\n")
    benchmark_plot_comparison(full_benchmark, method_id, prediction = "Final", same_scale = FALSE)
    #message("\nFinal Exposure Comparison...\n")
    exposure_plot <- benchmark_plot_exposures(full_benchmark, method_id, prediction = "Final")
    print(exposure_plot)
    
  }
  
  print(single_summary)
  
  message("\nDone.\n")
  
  if (make_copy == TRUE){
    return(full_benchmark)
  }
  
  
  
}

#' Predict and benchmark
#' 
#' This function will discover signatures from a musica object using a given
#' algorithm and a range of k values. A new discovery is done for each k value.
#' As each discovery is completed, the prediction is benchmarked.
#' 
#' @param musica A \code{\linkS4class{musica}} object.
#' @param table_name Name of the table to use for signature discovery. Needs
#' to be the same name supplied to the table building functions such as
#' \link{build_standard_table}.
#' @param algorithm Method to use for mutational signature discovery. One of 
#' \code{"lda"} or \code{"nmf"}. Default \code{"lda"}.
#' @param k_min Minimum number of singatures to predict
#' @param k_max Maximum number of signatures to predict
#' @param full_benchmark An object of class \code{\linkS4class{full_benchmark}}
#' created with the \link{create_benchmark} function or returned from a previous
#' \code{benchmark} run.
#' @param threshold Cosine similarity cutoff for comparing preidcted and true
#' signatures. Default \code{0.8}.
#' @param adjustment_threshold Cosine similarity value of high confidence.
#' Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold. Default \code{0.9}.
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
#' @export
predict_and_benchmark <- function(musica, table_name, algorithm = "lda", k_min, 
                                  k_max, full_benchmark, threshold = 0.8,
                                  adjustment_threshold = 0.9, plot = FALSE, 
                                  seed = 1, nstart = 10, par_cores = 1){
  
  
  for (k in c(k_min:k_max)){
    
    message("Discovering signatures for k = ", k, "...")
    
    prediction <- discover_signatures(musica = musica, table_name = table_name,
                                      num_signatures = k, algorithm = algorithm,
                                      seed = seed, nstart = nstart, 
                                      par_cores = par_cores)
    
    message("Benchmarking prediction for k = ", k, "...\n")
    
    full_benchmark <- benchmark(full_benchmark = full_benchmark, 
                                prediction = prediction, 
                                method_id = paste(algorithm, k, sep = ""),
                                threshold = threshold, 
                                adjustment_threshold = adjustment_threshold,
                                description = paste("Prediction from predict_and_benchmark function. Algorithm: ", algorithm, ". K = ", k, ".", sep = ""),
                                plot = plot, make_copy = TRUE)
  }
  
  
  print(full_benchmark@method_view_summary)
  
  return(full_benchmark)
  
}


#' Get a single_benchmark object
#' 
#' Access a \code{\linkS4class{single_benchmark}} object containing information from
#' an individual benchmark run from a \code{\linkS4class{full_benchmark}} object.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object that contains
#' the desired \code{\linkS4class{single_benchmark}} object
#' @param method_id The identifier for the desired \code{\linkS4class{single_benchmark}}
#' object
#'
#' @return A \code{\linkS4class{single_benchmark}} object
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
  
  benchmark <- indv_benchmarks(full_benchmark)[[method_id]]
  
  return(benchmark)
}


#' Get a benchmark prediction
#' 
#' Access a \code{\linkS4class{musica_result}} object containing a particular prediction
#' from a benchmarking analysis.
#'
#' @param indv_benchmark A \code{\linkS4class{single_benchmark}} object containing the
#' desired prediction. This can be accessed using the \link{benchmark_get_entry}
#' function.
#' @param prediction \code{Initial} for the prediction before any benchmarking
#' adjustments have been made, \code{Intermediate} for the prediction after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{Final} for the prediction at the end of the benchmarking adjustments. 
#'
#' @return A \code{\linkS4class{musica_result}} object
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


#' Plot signatures from a benchmarking analysis
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' this function can be used to plot signatures from any step in the benchmarking
#' process. Comparable to the \code{plot_signatures} function but compatible with
#' benchmarking objects.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the signatures to be plotted
#' @param prediction \code{Initial} for the signatures before any benchmarking
#' adjustments have been made, \code{Intermediate} for the signatures after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{Final} for the signatures at the end of the benchmarking adjustments.
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
#' will be adjusted for each signature. Default \code{TRUE}.
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
#' @export
benchmark_plot_signatures <- function(full_benchmark, method_id, prediction,
                                      plotly = FALSE, color_variable = NULL, color_mapping = NULL, text_size = 10,
                                      show_x_labels = TRUE, show_y_labels = TRUE, same_scale = TRUE, y_max = NULL, 
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
  
  # access musica object for desired prediction
  result <- benchmark_get_prediction(indv_benchmark, prediction)
  
  signatures_plot <- plot_signatures(result, plotly = plotly, color_variable = color_variable, color_mapping = color_mapping,
                                     text_size = text_size, show_x_labels = show_x_labels, show_y_labels = show_y_labels,
                                     same_scale = same_scale, y_max = y_max, annotation = annotation, percent = percent)
  
  return(signatures_plot)
  
}


#' Get benchmark comparison table
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' this function can be used to generate the comparison table between true and
#' predicted signatures from any step in the benchmarking process.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the comparison of interest
#' @param prediction \code{Initial} for the comparison before any benchmarking
#' adjustments have been made, \code{Intermediate} for the comparison after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{Final} for the comparison at the end of the benchmarking adjustments.
#'
#' @return A data.frame containing the comparison between true and predicted
#' signatures
#' @export
benchmark_compare_results <- function(full_benchmark, method_id, prediction){
  
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
  
  # access musica object for desired prediction
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


#' Plot a signature comparison from a benchmarking analysis
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' the comparison between the true and predicted signatures at any step of the
#' benchmarking process can be plotted.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object containing the comparison of interest
#' @param prediction \code{Initial} for the comparison before any benchmarking
#' adjustments have been made, \code{Intermediate} for the comparison after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{Final} for the comparison at the end of the benchmarking adjustments.
#' @param decimals Specifies rounding for similarity metric displayed. Default
#' \code{2}.
#' @param same_scale If \code{TRUE}, the scale of the probability for each
#' comparison will be the same. If \code{FALSE}, then the scale of the y-axis
#' will be adjusted for each comparison. Default \code{TRUE}.
#'
#' @return Returns the comparison plot
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
  
  comparison <- benchmark_compare_results(full_benchmark, method_id, prediction)
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  # access musica result object for desired prediction
  if (prediction == "Initial"){
    pred_res <- initial_pred(indv_benchmark)
    res_name <- c("Initial Predicted Signatures")
  }
  else if (prediction == "Intermediate"){
    pred_res <- intermediate_pred(indv_benchmark)
    res_name <- c("Intermediate Predicted Signatures")
  }
  else if (prediction == "Final"){
    pred_res <- final_pred(indv_benchmark)
    res_name <- c("Final Predicted Signatures")
  }
  
  
  comparison_plot <- plot_comparison(comparison, pred_res, ground_truth(full_benchmark),
                                     result_name = res_name, other_result_name = "True Signatures",
                                     decimals = decimals, same_scale = same_scale)
  
  #return(comparison_plot)
  
}


#' Plot exposure comparison from a benchmarking analysis
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' the comparison between the true and predicted exposures at any stage of the
#' benchmarking process can be plotted.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#' @param prediction \code{Initial} for the exposures before any benchmarking
#' adjustments have been made, \code{Intermediate} for the exposures after
#' duplicates have been adjusted but before composites are adjusted, or
#' \code{Final} for the exposures at the end of the benchmarking adjustments.
#'
#' @return Generates a  ggplot object
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
  comparison <- benchmark_compare_results(full_benchmark, method_id, prediction)
  
  truth <- ground_truth(full_benchmark)
  
  # access musica object for desired prediction
  prediction <- benchmark_get_prediction(indv_benchmark, prediction)
  
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
  
  plot_df <- data.frame(Predicted = predicted, True = true, Sig = sig)
  
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


#' Plot effect of duplicate correction
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' the true and predicted exposures can be plotted both before and after the
#' duplicate signature adjustment. The effect of the adjustment can then be
#' observed.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#'
#' @return A ggplot object
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
  
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  result_true <- ground_truth(full_benchmark)
  
  before <- initial_pred(indv_benchmark)
  after <- intermediate_pred(indv_benchmark)
  
  comparison <- initial_comparison(indv_benchmark)
  
  num_samples <- dim(exposures(result_true))[1]
  
  freq <- table(comparison$y_sig_name)
  duplicated_signatures <- names(freq[freq > 1])
  
  count <- 1
  
  for (duplicated_sig in duplicated_signatures){
    
    before_exposures <- NULL
    
    sigs_to_merge <- comparison[comparison$y_sig_name == duplicated_sig & !grepl("like", comparison$x_sig_name), 4]
    
    for (signature in sigs_to_merge){
      
      #tmp_exposures <- exposures(after)[signature,]
      tmp_exposures <- as.data.frame(as.vector(exposures(before)[signature,]))
      tmp_exposures_all <- cbind(tmp_exposures, 
                                 as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
      colnames(tmp_exposures_all) <- c("Predicted", "True")
      tmp_exposures_all$Source <- signature
      
      before_exposures <- rbind(before_exposures, tmp_exposures_all)
    }
    
    ## DOESNT WORK WHEN NOT ACTUALLY DUPICATE
    rownames(before_exposures) <- c(1: (num_samples * length(sigs_to_merge)))
    
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
    
    # exposure plot after merging
    
    tmp_exposures_list2 <- as.data.frame(as.vector(exposures(after)[new_sig_name,]))
    after_exposures <- cbind(tmp_exposures_list2, as.data.frame(as.numeric(exposures(result_true)[,duplicated_sig])))
    colnames(after_exposures) <- c("Predicted", "True")
    
    rownames(after_exposures) <- c(1:num_samples)
    
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
    
    if (count == 1){
      full_figure <- figure
    }
    else{
      full_figure <- ggpubr::ggarrange(full_figure, figure, ncol = 1, heights = c(count - 1,1))
    }
    
    count <- count + 1
    
    #print(figure)
    
    
  }
  
  #print(full_figure)
  if (count != 1){
    return(full_figure)
  }
  
}


#' Plot effect of composite correction
#' 
#' After a prediction has been benchmarked with the \link{benchmark} function,
#' the true and predicted exposures can be plotted both before and after the
#' composite signature adjustment. The effect of the adjustment can then be
#' observed.
#'
#' @param full_benchmark The \code{\linkS4class{full_benchmark}} object for the
#' benchmarking analysis
#' @param method_id The identifier for the \code{\linkS4class{single_benchmark}}
#' object of interest
#'
#' @return A ggplot object
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
  
  indv_benchmark <- benchmark_get_entry(full_benchmark, method_id)
  
  result_true <- ground_truth(full_benchmark)
  
  before <- intermediate_pred(indv_benchmark)
  after <- final_pred(indv_benchmark)
  
  comparison <- intermediate_comparison(indv_benchmark)
  
  num_samples <- dim(exposures(result_true))[1]
  
  freq <- table(comparison$x_sig_name)
  composite_signatures <- names(freq[freq > 1])
  
  count <- 1
  
  for (composite_sig in composite_signatures){
    
    before_exposures <- NULL
    
    # list of SBS signatures that are the components of this composite signature
    sig_components <- comparison[comparison$x_sig_name == composite_sig, 5]
    
    for (component_index in 1:length(sig_components)){
      
      # for plotting
      tmp_exp <- exposures(before)[composite_sig,]
      tmp_exp_list <- as.data.frame(as.vector(tmp_exp))
      tmp_exp_all <- cbind(tmp_exp_list, 
                           as.data.frame(as.numeric(exposures(result_true)[,sig_components[component_index]])))
      colnames(tmp_exp_all) <- c("Predicted", "True")
      tmp_exp_all$Source <- sig_components[component_index]
      
      before_exposures <- rbind(before_exposures, tmp_exp_all)
      
    }
    
    rownames(before_exposures) <- c(1: (num_samples * length(sig_components)))
    
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
    
    print(before_plot)
    
    colnames <- NULL
    for (component in sig_components){
      #colnames <- c(colnames, paste("Signature", component, "_like", sep = ""))
      colnames <- c(colnames, paste(component, "_like", sep = ""))
      
    }
    colnames <- colnames[colnames %in% rownames(exposures(after))]
    sig_components <- gsub("_like", "", colnames)
    
    after_exposures <- NULL
    
    for (component_index in 1:length(sig_components)){
      
      # for plotting
      tmp_exp <- exposures(after)[colnames[component_index],]
      tmp_exp_list <- as.data.frame(as.vector(tmp_exp))
      tmp_exp_all <- cbind(tmp_exp_list, 
                           as.data.frame(as.numeric(exposures(result_true)[,sig_components[component_index]])))
      colnames(tmp_exp_all) <- c("Predicted", "True")
      tmp_exp_all$Source <- sig_components[component_index]
      
      after_exposures <- rbind(after_exposures, tmp_exp_all)
      
    }
    
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
    
    print(after_plot)
    
    figure <- ggpubr::ggarrange(before_plot, after_plot, ncol = 2, nrow = 1)
    
    if (count == 1){
      full_figure <- figure
    }
    else{
      full_figure <- ggpubr::ggarrange(full_figure, figure, ncol = 1, heights = c(count - 1,1))
    }
    
    count <- count + 1
  }
  
  #combined_plot<- gridExtra::grid.arrange(before_plot, after_plot, ncol = 2)
  #return(combined_plot)
  
  if (count != 1){
    return(full_figure)
  }

  
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
.generate_summary <- function(title, result_all, result_true, comparison_results, count_table, final_musica, final_comparison){
  
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
  
  initial_reconstruction_error <- .get_reconstruction_error(result_all, count_table)
  initial_reconstruction_error <- round(initial_reconstruction_error, 3)
  end_reconstruction_error <- .get_reconstruction_error(final_musica, count_table)
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
.signature_view_summary <- function(benchmark){
  
  indv_benchmarks <- indv_benchmarks(benchmark)
  truth <- ground_truth(benchmark)
  
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
.update_benchmark <- function(full_benchmark, indv_benchmark, single_summary){
  
  # update summary
  if (dim(full_benchmark@method_view_summary)[1] == 0 ){
    full_benchmark@method_view_summary <- single_summary
  }else{
    full_benchmark@method_view_summary <- cbind(full_benchmark@method_view_summary, single_summary)
  }
  
  # update single benchmark list
  indv_benchmarks(full_benchmark) <- append(indv_benchmarks(full_benchmark), list(indv_benchmark))
  names(indv_benchmarks(full_benchmark))[length(indv_benchmarks(full_benchmark))] <- method_id(indv_benchmark)
  
  # update sig view summary
  sig_view_summary(full_benchmark) <- as.matrix(.signature_view_summary(full_benchmark))
  
  return(full_benchmark)
  
}


