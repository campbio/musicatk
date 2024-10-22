#' @import ggplot2
#' @importFrom conclust ckmeans
#' @importFrom cluster silhouette
#' @importFrom cluster sortSilhouette

#' @title Compare k values
#' @description Compare the stability and error of various k values to help
#' determine the correct number of signatures (k).
#'
#' @param musica A \code{\linkS4class{musica}} object.
#' @param modality The modality to use, either "SBS96", "DBS78", or "IND83".
#' @param reps Number of times prediction is performed. For each replicate, the
#' count table data is perturbed. Multiple replicates allows for stability
#' analysis by calculating silhouette width on the multiple results. Default
#' \code{100}.
#' @param min_k Lower range of number of signatures for discovery. Default
#' \code{1}.
#' @param max_k Upper range of number of signatures for discovery. Default
#' \code{10}.
#' @param error_type Whether to calculate reconstruction error by proportions
#' ("prop") or raw counts ("raw"). Default \code{"prop"}.
#' @param algorithm Algorithm for signature discovery. Default \code{"nmf"}.
#'
#' @return a data.frame with stats for each k value tested
#' @export
#'
#' @examples
#' data(musica)
#' compare_k_vals(musica, "SBS96", reps = 3, min_k = 1, max_k = 5)
compare_k_vals <- function(musica, modality, reps = 100, min_k = 1, max_k = 10,
                           error_type = "prop", algorithm = "nmf") {
  k_list <- c(min_k:max_k)

  # access count table
  count_table <- extract_count_tables(musica)[[modality]]
  count_table <- get_count_table(count_table)

  # generate bootstraped count tables
  boot_tables_all <- .gen_all_boots(count_table, reps)

  # establish df for sil widths
  results_df <- data.frame(
    k = NULL, sil_width = NULL, min_sil_width = NULL,
    max_sil_width = NULL, sterror_sil_width = NULL,
    error = NULL, min_error = NULL, max_error = NULL,
    sterror_error = NULL
  )

  # establish list for predicted signatures
  signatures_list <- list()

  # loop through the k values to test
  for (index in seq_len(length(k_list))) {
    message("\nk = ", k_list[index], ":\n", sep = "")

    # run discovery
    iterative_results <- .iterative_discovery(
      musica, modality, reps,
      k_list[index], algorithm,
      boot_tables_all
    )

    # result from unpreterbed count table
    result_actual <- iterative_results[[1]]

    # all signatures predicted
    signatures <- iterative_results[[2]]

    # all musica result objects
    predicted_results_list <- iterative_results[[3]]

    # perform clustering on all signatures
    if (k_list[index] > 1) {
      ckmeans <- .run_modified_kmeans(reps, k_list[index], signatures)

      # calculate sil widths
      sil_width <- .get_sil_width(signatures, ckmeans)
      avg_sil_width <- sil_width[[1]]
      min_sil_width <- sil_width[[2]]
      max_sil_width <- sil_width[[3]]
      sd_sil_width <- sil_width[[4]]
      sterror_sil_width <- sd_sil_width / sqrt(k_list[index] * (reps + 1))
    }

    # store signatures
    signatures_list[[paste(k_list[index], "Signatures")]] <- signatures

    # calculate RE
    RE <- .get_average_error(
      result_actual, predicted_results_list, error_type,
      count_table, boot_tables_all
    )
    avg_RE <- RE[[1]]
    min_RE <- RE[[2]]
    max_RE <- RE[[3]]
    sd_RE <- RE[[4]]
    sterror_RE <- sd_RE / sqrt(reps + 1)

    # add to dataframe
    if (k_list[index] > 1) {
      results_df <- rbind(
        results_df,
        data.frame(
          k = k_list[index],
          sil_width = avg_sil_width,
          min_sil_width = min_sil_width,
          max_sil_width = max_sil_width,
          sd_sil_width = sd_sil_width,
          sterror_sil_width = sterror_sil_width,
          error = avg_RE,
          min_error = min_RE,
          max_error = max_RE,
          sd_error = sd_RE,
          sterror_error = sterror_RE
        )
      )
    } else {
      results_df <- rbind(
        results_df,
        data.frame(
          k = k_list[index],
          sil_width = NA,
          min_sil_width = NA,
          max_sil_width = NA,
          sd_sil_width = NA,
          sterror_sil_width = NA,
          error = avg_RE,
          min_error = min_RE,
          max_error = max_RE,
          sd_error = sd_RE,
          sterror_error = sterror_RE
        )
      )
    }
  }

  figure <- plot_k_comparison(results_df)
  print(figure)

  return(results_df)
}


#' @title Plot k comparison
#' @description Plot the results of comparing k values
#'
#' @param k_comparison data.frame with k value comparisons returned from the
#' \code{\link{compare_k_vals}} function.
#'
#' @return a ggplot figure
#' @export
#'
#' @examples
#' data(musica)
#' k_comparison <- compare_k_vals(musica, "SBS96",
#'   reps = 3, min_k = 1,
#'   max_k = 5
#' )
#' plot_k_comparison(k_comparison)
plot_k_comparison <- function(k_comparison) {
  k_vals <- k_comparison$k
  errors <- k_comparison$error
  min_errors <- k_comparison$min_error
  max_errors <- k_comparison$max_error
  sd_errors <- k_comparison$sd_error
  se_errors <- k_comparison$sterror_error
  sws <- k_comparison$sil_width
  min_sws <- k_comparison$min_sil_width
  max_sws <- k_comparison$max_sil_width
  sd_sws <- k_comparison$sd_sil_width
  se_sws <- k_comparison$sterror_sil_width

  errors_to_plot <- (errors - min(min_errors)) /
    (max(max_errors) - min(min_errors))
  se_errors_to_plot <- (se_errors * errors_to_plot) / errors
  sd_errors_to_plot <- (sd_errors * errors_to_plot) / errors

  range <- max(max_errors) - min(min_errors)
  min <- min(min_errors)

  data <- data.frame(
    k_vals = k_vals, errors = errors_to_plot,
    sd_errors = sd_errors_to_plot, sws = sws,
    sd_sws = sd_sws
  )

  re_color <- "red3"
  sw_color <- "blue2"

  sw_top <- 1
  if (max((sws + sd_sws)[c(-1)]) > 1) {
    sw_top <- max((sws + sd_sws)[c(-1)]) + 0.01
  }

  sw_bottom <- 0
  if (min((sws - sd_sws)[c(-1)]) < 0) {
    sw_bottom <- min((sws - sd_sws)[c(-1)]) - 0.01
  }

  if (k_vals[1] == 1) {
    modifier <- 0
  } else {
    modifier <- 1
  }

  figure <- ggplot(data, aes(x = k_vals - modifier)) +
    geom_line(aes(y = sws), color = sw_color, linetype = "dashed") +
    geom_line(aes(y = errors), color = re_color) +
    geom_point(aes(y = sws), color = sw_color) +
    geom_point(aes(y = errors), color = re_color) +
    geom_errorbar(
      aes(
        ymin = (errors - sd_errors),
        ymax = (errors + sd_errors)
      ),
      width = 0.2, color = re_color, alpha = 0.3
    ) +
    geom_errorbar(aes(ymin = sws - sd_sws, ymax = sws + sd_sws),
      width = 0.2, color = sw_color, alpha = 0.3
    ) +
    scale_y_continuous(
      # Features of the first axis
      name = "Silhouette Width (-)",
      limits = c(sw_bottom, sw_top),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~ . * 1,
        name = "Reconstruction Error",
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c(
          round(0 * range + min, 2),
          round(0.25 * range + min, 2),
          round(0.5 * range + min, 2),
          round(0.75 * range + min, 2),
          round(1 * range + min, 2)
        )
      )
    ) +
    scale_x_discrete(limits = factor(k_vals), name = "k value") +
    labs(x = "Number of Signatures (k)") +
    theme_light() +
    theme(
      axis.title.y = element_text(color = sw_color, size = 13),
      axis.title.y.right = element_text(color = re_color, size = 13),
      axis.title.x = element_text(size = 13),
      axis.text = element_text(color = "black", size = 10)
    )

  return(figure)
}

# sample same number of samples, with replacement
.resample_count_table <- function(count_table) {
  colsums <- colSums(count_table)
  prob <- t(t(count_table) / colsums)
  prob[prob == NaN] <- 0
  bootstrapped_count_table <- vapply(seq_len(ncol(count_table)), function(idx) {
    stats::rmultinom(
      n = 1,
      size = colsums[idx],
      prob = prob[, idx]
    )
  }, integer(nrow(count_table)))

  rownames(bootstrapped_count_table) <- rownames(count_table)
  colnames(bootstrapped_count_table) <- colnames(count_table)

  # return bootstrapped count table
  return(bootstrapped_count_table)
}



# generate all bootstrapped count tables that will be used
.gen_all_boots <- function(count_table, reps) {
  boot_tables_all <- NULL

  for (iteration in seq_len(reps)) {
    boot_tables_all[[iteration]] <- .resample_count_table(count_table)
  }

  return(boot_tables_all)
}

# get signatures for set number of iterations on bootstrapped
.iterative_discovery <- function(musica, modality, reps, k, alg,
                                 boot_tables_all) {
  message("Calculating signatures on original count table...\n\n")
  result_discov_actual <- discover_signatures(musica,
    modality = modality,
    num_signatures = k,
    algorithm = alg, nstart = 5,
    model_id = "actual",
    result_name = "k_val_assist",
    make_copy = TRUE
  )

  result_discov_actual <- get_model(
    result_discov_actual, "k_val_assist",
    modality, "actual"
  )

  # create object for storing signatures from all iterations
  result_discovs_boot <- NULL
  signatures <- signatures(result_discov_actual)

  # loop through given number of iterations
  for (iteration in seq_len(reps)) {
    message("Replicate ", iteration, " of ", reps, "...\n")

    # copy musica object
    musica_boot <- musica

    # replace count table with bootstrapped count table
    boot_table <- as.matrix(boot_tables_all[[iteration]])
    class(boot_table) <- "numeric"
    musica_boot@count_tables[[modality]]@count_table <- boot_table

    # find signatures
    result_discov_boot <- discover_signatures(musica_boot,
      modality = modality,
      num_signatures = k,
      algorithm = alg, nstart = 5,
      model_id = "boot",
      result_name = "k_val_assist",
      make_copy = TRUE
    )

    # add result_discov to list
    result_discovs_boot[[iteration]] <- get_model(
      result_discov_boot,
      "k_val_assist", modality,
      "boot"
    )

    # add signatures to matrix of al signatures
    signatures <- cbind(signatures, signatures(
      result_discov_boot,
      "k_val_assist", modality, "boot"
    ))
  }

  # rename columns to be numbered
  colnames <- NULL
  for (i in 0:reps) {
    for (s in seq_len(k)) {
      colnames <- c(colnames, paste(i, ".", s, sep = ""))
    }
  }
  colnames(signatures) <- colnames

  # return the signatures
  return(list(result_discov_actual, signatures, result_discovs_boot))
}

# Run k-means
.run_modified_kmeans <- function(reps, k, signatures) {
  # establish matrix for pairs of signatures that can't be in same cluster
  cant_link <- NULL

  # Take sequence of 1 to total number of signatures (number of iterations times
  # number of signatures generated for each iteration, plus number of signatures
  # per iteration for actual/original signatures) and split it into a list,
  # where each element is a vector listing the signature indexes for one
  # iteration. (List of length (number of iterations + 1) where each element is
  # a vector length (number of signatures generated per iteration))
  sig_codes_for_each_iter <- split(
    seq(1, (reps * k) + k),
    ceiling(seq_along(seq(
      1,
      (reps * k) + k
    )) / k)
  )

  # fill matrix with pairwise signature indexes that must be in separate
  # clusters
  for (i in 1:(reps + 1)) {
    cant_link <- rbind(cant_link, t(combn(sig_codes_for_each_iter[[i]], 2)))
  }

  # must link constraint is empty matrix
  must_link <- matrix(c(NA, NA), nrow = 1)

  ckmeans <- conclust::ckmeans(
    data = t(signatures), k = k,
    mustLink = must_link, cantLink = cant_link
  )

  return(ckmeans)
}

# Calculate silhouette width
.get_sil_width <- function(signatures, ckmeans) {
  # calculate distance matrix
  dists <- as.matrix(dist(t(signatures), method = "manhattan"))

  # compute silhouette widths
  sil <- cluster::silhouette(ckmeans, dists)
  # plot(sil)
  sil_summary <- summary(sil)
  avg_sil <- sil_summary$avg.width
  min_sil <- sil_summary$si.summary[[1]]
  max_sil <- sil_summary$si.summary[[6]]

  sorted <- cluster::sortSilhouette(sil)
  sd_sil <- sd(sorted[, 3])

  return(list(avg_sil, min_sil, max_sil, sd_sil))
}

# calculate proportional reconstruction error
.get_reconstruction_error_prop <- function(result_discov, count_table) {
  # extract exposures and signatures matrices
  expos <- exposures(result_discov)
  sigs <- signatures(result_discov)

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

# calculate raw reconstruction error
.get_reconstruction_error_raw <- function(result_discov, count_table) {
  # extract exposures and signatures matrices
  expos <- exposures(result_discov)
  sigs <- signatures(result_discov)

  # multiply signature and exposure matrices
  predicted_counts <- sigs %*% expos

  # calculate sum of differences (reconstruction error)
  reconstruction_error <- sum(abs(count_table - predicted_counts))

  # return reconstruction error
  return(reconstruction_error)
}

# calculate recontruction error with a given error type
.get_reconstruction_error <- function(result_discov, count_table, error_type) {
  if (error_type %in% c(
    "prop", "Prop", "proportion", "Proportion",
    "proportions", "Proportions", "proportional",
    "Proportional"
  )) {
    error <- .get_reconstruction_error_prop(result_discov, count_table)
  } else if (error_type %in% c("raw", "Raw", "raws", "Raws")) {
    error <- .get_reconstruction_error_raw(result_discov, count_table)
  } else {
    message("Invalid error type. Please enter 'prop' for proprotional RE or
            'raw' for raw RE.")
  }
}

# get average reconstruction error for a k value from all bootstrap iterations
.get_average_error <- function(result_actual, predicted_results_list,
                               error_type, count_table, boot_tables_all) {
  # initialize error vector
  errors <- NULL

  # calculate error for unpreterbed count table prediction
  error <- .get_reconstruction_error(result_actual, count_table, error_type)

  # append to errors
  errors <- c(errors, error)

  # loop through the bootstrapped predictions
  for (index in seq_len(length(predicted_results_list))) {
    # access the musica result object
    prediction <- predicted_results_list[[index]]

    # calcualte error
    error <- .get_reconstruction_error(
      prediction, boot_tables_all[[index]],
      error_type
    )

    # append error
    errors <- c(errors, error)
  }

  # calculate average error
  mean_error <- mean(errors)

  # min and max errors
  min_error <- min(errors)
  max_error <- max(errors)

  # standard deviation
  sd_error <- sd(errors)

  return(list(mean_error, min_error, max_error, sd_error))
}
