% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_results.R
\name{compare_results}
\alias{compare_results}
\title{Compare two result files to find similar signatures}
\usage{
compare_results(
  musica,
  model_id,
  other_model_id,
  modality = "SBS96",
  result_name = "result",
  other_musica = NULL,
  other_result_name = "result",
  threshold = 0.9,
  metric = "cosine",
  result_rename = deparse(substitute(model_id)),
  other_result_rename = deparse(substitute(other_model_id)),
  decimals = 2,
  same_scale = FALSE
)
}
\arguments{
\item{musica}{A \code{\linkS4class{musica}} object.}

\item{model_id}{The name of the first model to compare.}

\item{other_model_id}{The name of the second model to compare.}

\item{modality}{Modality of results being compared. Default \code{"SBS96"}.}

\item{result_name}{Name of the result list entry for the first result to
compare. Default \code{"result"}.}

\item{other_musica}{A second \code{\linkS4class{musica}} object. If null, the
provided musica object is used twice. Default \code{NULL}.}

\item{other_result_name}{Name of the result list entry for the second result
to compare. Default \code{"result"}.}

\item{threshold}{threshold for similarity}

\item{metric}{One of \code{"cosine"} for cosine similarity or \code{"jsd"}
for 1 minus the Jensen-Shannon Divergence. Default \code{"cosine"}.}

\item{result_rename}{title for plot of first result signatures}

\item{other_result_rename}{title for plot of second result signatures}

\item{decimals}{Specifies rounding for similarity metric displayed. Default
\code{2}.}

\item{same_scale}{If \code{TRUE}, the scale of the probability for each
signature will be the same. If \code{FALSE}, then the scale of the y-axis
will be adjusted for each signature. Default \code{FALSE}.}
}
\value{
Returns the comparisons
}
\description{
Compare two result files to find similar signatures
}
\examples{
data(res)
compare_results(res,
  model_id = "res", other_model_id = "res",
  modality = "SBS96", threshold = 0.8
)
}
