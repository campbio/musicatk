# Single Benchmark object/methods -------------------------------

#' @title character or NULL class union
#' @description class to allow either NULL or character
#' @keywords internal
#' @noRd
setClassUnion("CharOrNULL", c("character", "NULL"))

#' Object that contains information for an individual benchmark run
#' 
#' @slot initial_pred The initial prediction
#' @slot intermediate_pred The intermediate prediction, after duplicates are
#' corrected, but before composites are corrected
#' @slot final_pred The final prediction
#' @slot initial_comparison The initial comparison between predicted and true
#' signatures
#' @slot intermediate_comparison The intermediate comparison between predicted
#' and true signatures
#' @slot final_comparison The final comparison between predicted and true
#' signatures
#' @slot single_summary A matrix summarizing the benchmark run
#' @slot method_id An identifier for the benchmark run
#' @slot threshold Cosine similarity cutoff for comparing predicted and true
#' signatures
#' @slot adjustment_threshold Cosine similarity value of high confidence.
#' Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold.
#' @slot description Further details about the prediction being benchmarked.
#' @export
#' @exportClass single_benchmark

setClass(
  "single_benchmark",
  slots = list(
    initial_pred = "result_model",
    intermediate_pred = "result_model",
    final_pred = "result_model",
    initial_comparison = "data.frame",
    intermediate_comparison = "data.frame",
    final_comparison = "data.frame",
    single_summary = "matrix",
    method_id = "character",
    threshold = "numeric",
    adjustment_threshold = "numeric",
    description = "CharOrNULL"
  )
)

# Getter and setter functions 

# --------initial_pred --------------
#' @title Retrieve initial prediction from a single_benchmark object
#' @description  The \code{initial_prediction} contains the signatures and
#' exposures before duplicates or composites are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname initial_pred
#' @return An object of class \code{\linkS4class{result_model}} containing
#' signatures and exposures from the initial prediction, before duplicates
#' and composites are corrected.
#' @export
setGeneric(
  name = "initial_pred",
  def = function(x) {
    standardGeneric("initial_pred")
  }
)

#' @rdname initial_pred
setMethod(
  f = "initial_pred",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@initial_pred)
  }
)

#' @rdname initial_pred
#' @param x A \code{\linkS4class{single_benchmark} }object.
#' @param value A \code{\linkS4class{result_model}} object containing the
#' signatures and exposures from the initial prediction.
#' @export
setGeneric(
  name = "initial_pred<-",
  def = function(x, value) {
    standardGeneric("initial_pred<-")
  }
)

#' @rdname initial_pred
setReplaceMethod(
  f = "initial_pred",
  signature = c("single_benchmark", "result_model"),
  definition = function(x, value) {
    x@initial_pred <- value
    return(x)
  }
)

# --------intermediate_pred --------------
#' @title Retrieve intermediate prediction from a single_benchmark object
#' @description  The \code{intermediate_prediction} contains the signatures and
#' exposures after duplicates are corrected but before composites are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname intermediate_pred
#' @return An object of class \code{\linkS4class{result_model}} containing
#' signatures and exposures from the intermediate prediction.
#' @export
setGeneric(
  name = "intermediate_pred",
  def = function(x) {
    standardGeneric("intermediate_pred")
  }
)

#' @rdname intermediate_pred
setMethod(
  f = "intermediate_pred",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@intermediate_pred)
  }
)

#' @rdname intermediate_pred
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A \code{\linkS4class{result_model}} object containing the
#' signatures and exposures from the intermediate prediction.
#' @export
setGeneric(
  name = "intermediate_pred<-",
  def = function(x, value) {
    standardGeneric("intermediate_pred<-")
  }
)

#' @rdname intermediate_pred
setReplaceMethod(
  f = "intermediate_pred",
  signature = c("single_benchmark", "result_model"),
  definition = function(x, value) {
    x@intermediate_pred <- value
    return(x)
  }
)

# --------final_pred --------------
#' @title Retrieve final prediction from a single_benchmark object
#' @description  The \code{final_prediction} contains the signatures and
#' exposures after duplicates and composites are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname final_pred
#' @return An object of class \code{\linkS4class{result_model}} containing
#' signatures and exposures from the final prediction.
#' @export
setGeneric(
  name = "final_pred",
  def = function(x) {
    standardGeneric("final_pred")
  }
)

#' @rdname final_pred
setMethod(
  f = "final_pred",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@final_pred)
  }
)

#' @rdname final_pred
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A \code{\linkS4class{result_model}} object containing the
#' signatures and exposures from the final prediction.
#' @export
setGeneric(
  name = "final_pred<-",
  def = function(x, value) {
    standardGeneric("final_pred<-")
  }
)

#' @rdname final_pred
setReplaceMethod(
  f = "final_pred",
  signature = c("single_benchmark", "result_model"),
  definition = function(x, value) {
    x@final_pred <- value
    return(x)
  }
)

# --------initial_comparison --------------
#' @title Retrieve initial comparison between the prediction and the ground
#' truth from a single_benchmark object
#' @description  The \code{initial_comparison} contains the comparison between
#' the ground truth and the initial prediction, before duplicates or composites
#' are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname initial_comparison
#' @return A data.frame containing the comparisons.
#' @export
setGeneric(
  name = "initial_comparison",
  def = function(x) {
    standardGeneric("initial_comparison")
  }
)

#' @rdname initial_comparison
setMethod(
  f = "initial_comparison",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@initial_comparison)
  }
)

#' @rdname initial_comparison
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A data.frame containing the comparisons between the intial
#' prediction and the ground truth.
#' @export
setGeneric(
  name = "initial_comparison<-",
  def = function(x, value) {
    standardGeneric("initial_comparison<-")
  }
)

#' @rdname initial_comparison
setReplaceMethod(
  f = "initial_comparison",
  signature = c("single_benchmark", "data.frame"),
  definition = function(x, value) {
    x@initial_comparison <- value
    return(x)
  }
)

# --------intermediate_comparison --------------
#' @title Retrieve intermediate comparison between the prediction and the ground
#' truth from a single_benchmark object
#' @description  The \code{intermediate_comparison} contains the comparison
#' between the ground truth and the intermediate prediction, after duplicates
#' are corrected but before composites are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname intermediate_comparison
#' @return A dataframe containing the comparisons.
#' @export
setGeneric(
  name = "intermediate_comparison",
  def = function(x) {
    standardGeneric("intermediate_comparison")
  }
)

#' @rdname intermediate_comparison
setMethod(
  f = "intermediate_comparison",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@intermediate_comparison)
  }
)

#' @rdname intermediate_comparison
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A data.frame containing the comparisons between the intermediate
#' prediction and the ground truth.
#' @export
setGeneric(
  name = "intermediate_comparison<-",
  def = function(x, value) {
    standardGeneric("intermediate_comparison<-")
  }
)

#' @rdname intermediate_comparison
setReplaceMethod(
  f = "intermediate_comparison",
  signature = c("single_benchmark", "data.frame"),
  definition = function(x, value) {
    x@intermediate_comparison <- value
    return(x)
  }
)

# --------final_comparison --------------
#' @title Retrieve final comparison between the prediction and the ground
#' truth from a single_benchmark object
#' @description  The \code{final_comparison} contains the comparison between
#' the ground truth and the final prediction, after duplicates and composites
#' are corrected.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname final_comparison
#' @return A data.frame containing the comparisons.
#' @export
setGeneric(
  name = "final_comparison",
  def = function(x) {
    standardGeneric("final_comparison")
  }
)

#' @rdname final_comparison
setMethod(
  f = "final_comparison",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@final_comparison)
  }
)

#' @rdname final_comparison
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A data.frame containing the comparisons between the intial
#' prediction and the ground truth.
#' @export
setGeneric(
  name = "final_comparison<-",
  def = function(x, value) {
    standardGeneric("final_comparison<-")
  }
)

#' @rdname final_comparison
setReplaceMethod(
  f = "final_comparison",
  signature = c("single_benchmark", "data.frame"),
  definition = function(x, value) {
    x@final_comparison <- value
    return(x)
  }
)

# --------single_summary --------------
#' @title Retrieve the single summary matrix from a single_benchmark object
#' @description  The \code{single_summary} contains the summary information
#' from the benchmark run.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname single_summary
#' @return A matrix containing the summary.
#' @export
setGeneric(
  name = "single_summary",
  def = function(x) {
    standardGeneric("single_summary")
  }
)

#' @rdname single_summary
setMethod(
  f = "single_summary",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@single_summary)
  }
)

#' @rdname single_summary
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value A matrix containing the summary infomation for the run.
#' @export
setGeneric(
  name = "single_summary<-",
  def = function(x, value) {
    standardGeneric("single_summary<-")
  }
)

#' @rdname single_summary
setReplaceMethod(
  f = "single_summary",
  signature = c("single_benchmark", "matrix"),
  definition = function(x, value) {
    x@single_summary <- value
    return(x)
  }
)

# --------method_id --------------
#' @title Retrieve the method id from a single_benchmark object
#' @description  The \code{method_id} is the identifier for the benchmark run.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname method_id
#' @return A matrix containing the summary.
#' @export
setGeneric(
  name = "method_id",
  def = function(x) {
    standardGeneric("method_id")
  }
)

#' @rdname method_id
setMethod(
  f = "method_id",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@method_id)
  }
)

#' @rdname method_id
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value An identifier for the benchmark run.
#' @export
setGeneric(
  name = "method_id<-",
  def = function(x, value) {
    standardGeneric("method_id<-")
  }
)

#' @rdname method_id
setReplaceMethod(
  f = "method_id",
  signature = c("single_benchmark", "character"),
  definition = function(x, value) {
    x@method_id <- value
    return(x)
  }
)

# --------threshold --------------
#' @title Retrieve the threshold from a single_benchmark object
#' @description  The \code{threshold} is the cosine similarity cutoff for
#' comparing predicted and true signatures. 
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname threshold
#' @return The threshold value.
#' @export
setGeneric(
  name = "threshold",
  def = function(x) {
    standardGeneric("threshold")
  }
)

#' @rdname threshold
setMethod(
  f = "threshold",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@threshold)
  }
)

#' @rdname threshold
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value The threshold
#' @export
setGeneric(
  name = "threshold<-",
  def = function(x, value) {
    standardGeneric("threshold<-")
  }
)

#' @rdname threshold
setReplaceMethod(
  f = "threshold",
  signature = c("single_benchmark", "numeric"),
  definition = function(x, value) {
    x@threshold <- value
    return(x)
  }
)

# --------adjustment_threshold --------------
#' @title Retrieve the adjustment_threshold from a single_benchmark object
#' @description  The \code{adjustment_threshold} is a cosine similarity value of
#' high confidence.Comparisons that meet this cutoff are assumed to be likely,
#' while those that fall below the cutoff will be disregarded if the predicted
#' signature is already captured above the threshold.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname adjustment_threshold
#' @return The adjustment_threshold value.
#' @export
setGeneric(
  name = "adjustment_threshold",
  def = function(x) {
    standardGeneric("adjustment_threshold")
  }
)

#' @rdname adjustment_threshold
setMethod(
  f = "adjustment_threshold",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@adjustment_threshold)
  }
)

#' @rdname adjustment_threshold
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value The adjustment_threshold
#' @export
setGeneric(
  name = "adjustment_threshold<-",
  def = function(x, value) {
    standardGeneric("adjustment_threshold<-")
  }
)

#' @rdname adjustment_threshold
setReplaceMethod(
  f = "adjustment_threshold",
  signature = c("single_benchmark", "numeric"),
  definition = function(x, value) {
    x@adjustment_threshold <- value
    return(x)
  }
)

# --------description --------------
#' @title Retrieve the description from a single_benchmark object
#' @description  The \code{description} contains further details about the
#' prediction being benchmarked.
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @rdname description
#' @return The description
#' @export
setGeneric(
  name = "description",
  def = function(x) {
    standardGeneric("description")
  }
)

#' @rdname description
setMethod(
  f = "description",
  signature = "single_benchmark",
  definition = function(x) {
    return(x@description)
  }
)

#' @rdname description
#' @param x A \code{\linkS4class{single_benchmark}} object.
#' @param value The description
#' @export
setGeneric(
  name = "description<-",
  def = function(x, value) {
    standardGeneric("description<-")
  }
)

#' @rdname description
setReplaceMethod(
  f = "description",
  signature = c("single_benchmark", "CharOrNULL"),
  definition = function(x, value) {
    x@description <- value
    return(x)
  }
)






