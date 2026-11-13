#' @include class_musica.R
NULL

# Full Benchmark object/methods -------------------------------

#' Object that contains information for a full benchmarking analysis where
#' multiple predictions may be benchmarked
#'
#' @slot ground_truth A \code{\linkS4class{musica}} object containing the true
#' signatures and exposures.
#' @slot method_view_summary A data.frame summarizing each benchmark run,
#' with one row per run. See \link{method_view_summary}.
#' @slot sig_view_summary A data.frame summarizing each true signature across
#' all benchmark runs, with one row per signature. See
#' \link{sig_view_summary}.
#' @slot indv_benchmarks A list of
#' \code{\linkS4class{single_benchmark}} objects, each containing information
#' from a single benchmark run.
#' @export
#' @exportClass full_benchmark

setClass(
  "full_benchmark",
  slots = list(
    ground_truth = "musica",
    method_view_summary = "data.frame",
    sig_view_summary = "data.frame",
    indv_benchmarks = "list"
  )
)

# Getter and setter functions

# --------ground_truth --------------
#' @title Retrieve the ground_truth from a full_benchmark object
#' @description  The \code{ground_truth} is a \code{\linkS4class{musica}} object
#' that contains the true signatures and exposures.
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @rdname ground_truth
#' @return A \code{\linkS4class{musica}} object containing the ground truth.
#' @export
setGeneric(
  name = "ground_truth",
  def = function(x) {
    standardGeneric("ground_truth")
  }
)

#' @rdname ground_truth
setMethod(
  f = "ground_truth",
  signature = "full_benchmark",
  definition = function(x) {
    return(x@ground_truth)
  }
)

#' @rdname ground_truth
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @param value A \code{\linkS4class{musica}} object containing the ground
#' truth.
#' @export
setGeneric(
  name = "ground_truth<-",
  def = function(x, value) {
    standardGeneric("ground_truth<-")
  }
)

#' @rdname ground_truth
setReplaceMethod(
  f = "ground_truth",
  signature = c("full_benchmark", "musica"),
  definition = function(x, value) {
    x@ground_truth <- value
    return(x)
  }
)

# --------method_view_summary --------------
#' @title Retrieve the method_view_summary from a full_benchmark object
#' @description  The \code{method_view_summary} is a data.frame with one
#' row per benchmark run (the row names are the method IDs). Its columns are
#' described in \link{benchmark}.
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @rdname method_view_summary
#' @return A data.frame containing the method view summary
#' @export
setGeneric(
  name = "method_view_summary",
  def = function(x) {
    standardGeneric("method_view_summary")
  }
)

#' @rdname method_view_summary
setMethod(
  f = "method_view_summary",
  signature = "full_benchmark",
  definition = function(x) {
    return(x@method_view_summary)
  }
)

#' @rdname method_view_summary
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @param value A data.frame containing the method_view_summary.
#' @export
setGeneric(
  name = "method_view_summary<-",
  def = function(x, value) {
    standardGeneric("method_view_summary<-")
  }
)

#' @rdname method_view_summary
setReplaceMethod(
  f = "method_view_summary",
  signature = c("full_benchmark", "data.frame"),
  definition = function(x, value) {
    x@method_view_summary <- value
    return(x)
  }
)

# --------sig_view_summary --------------
#' @title Retrieve the sig_view_summary from a full_benchmark object
#' @description  The \code{sig_view_summary} is a data.frame with one row
#' per true signature (the row names are the signature names), summarizing
#' the Final stage of every benchmark run: the number of runs that found and
#' missed the signature (\code{"Times Found"}, \code{"Times Missed"}), the
#' mean, minimum and maximum over those runs of its best cosine similarity
#' (\code{"Mean CS"}, \code{"Min CS"}, \code{"Max CS"}; \code{NA} if it was
#' never found), and the runs that found it (\code{"Methods Found"}).
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @rdname sig_view_summary
#' @return A data.frame containing the signature view summary
#' @export
setGeneric(
  name = "sig_view_summary",
  def = function(x) {
    standardGeneric("sig_view_summary")
  }
)

#' @rdname sig_view_summary
setMethod(
  f = "sig_view_summary",
  signature = "full_benchmark",
  definition = function(x) {
    return(x@sig_view_summary)
  }
)

#' @rdname sig_view_summary
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @param value A data.frame containing the sig_view_summary
#' @export
setGeneric(
  name = "sig_view_summary<-",
  def = function(x, value) {
    standardGeneric("sig_view_summary<-")
  }
)

#' @rdname sig_view_summary
setReplaceMethod(
  f = "sig_view_summary",
  signature = c("full_benchmark", "data.frame"),
  definition = function(x, value) {
    x@sig_view_summary <- value
    return(x)
  }
)

# --------indv_benchmarks --------------
#' @title Retrieve the indv_benchmarks list from a full_benchmark object
#' @description  The \code{indv_benchmarks} list is a list of
#' \code{\linkS4class{single_benchmark}} objects, each containing information
#' about a single benchmark run.
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @rdname indv_benchmarks
#' @return A list of \code{\linkS4class{single_benchmark}} objects
#' @export
setGeneric(
  name = "indv_benchmarks",
  def = function(x) {
    standardGeneric("indv_benchmarks")
  }
)

#' @rdname indv_benchmarks
setMethod(
  f = "indv_benchmarks",
  signature = "full_benchmark",
  definition = function(x) {
    return(x@indv_benchmarks)
  }
)

#' @rdname indv_benchmarks
#' @param x A \code{\linkS4class{full_benchmark}} object.
#' @param value A list of \code{\linkS4class{single_benchmark}} objects.
#' @export
setGeneric(
  name = "indv_benchmarks<-",
  def = function(x, value) {
    standardGeneric("indv_benchmarks<-")
  }
)

#' @rdname indv_benchmarks
setReplaceMethod(
  f = "indv_benchmarks",
  signature = c("full_benchmark", "list"),
  definition = function(x, value) {
    x@indv_benchmarks <- value
    return(x)
  }
)
