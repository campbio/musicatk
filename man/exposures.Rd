% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_result_model.R
\name{exposures}
\alias{exposures}
\alias{exposures,musica-method}
\alias{exposures,result_collection-method}
\alias{exposures,result_model-method}
\alias{exposures<-}
\alias{exposures<-,musica,matrix-method}
\alias{exposures<-,result_collection,matrix-method}
\alias{exposures<-,result_model,matrix-method}
\title{Retrieve exposures from a result_model, result_collection, or musica
object}
\usage{
exposures(x, ...)

\S4method{exposures}{musica}(x, result, modality, model_id)

\S4method{exposures}{result_collection}(x, modality, model_id)

\S4method{exposures}{result_model}(x)

exposures(x, ...) <- value

\S4method{exposures}{musica,matrix}(x, result, modality, model_id) <- value

\S4method{exposures}{result_collection,matrix}(x, modality, model_id) <- value

\S4method{exposures}{result_model,matrix}(x) <- value
}
\arguments{
\item{x}{A \code{\linkS4class{musica}},
\code{\linkS4class{result_collection}}, or \code{\linkS4class{result_model}}
object generated by a mutational discovery or prediction tool.}

\item{...}{Other inputs}

\item{result}{Name of result from result_list to assign the exposures.
Used when \code{result} is a \code{\linkS4class{musica}} object.}

\item{modality}{Modality to assign the exposures. Used when
\code{result} is a \code{\linkS4class{musica}} or
\code{\linkS4class{result_collection}} object.}

\item{model_id}{Model identifier to assign the exposures. Used when
\code{result} is a \code{\linkS4class{musica}} or
\code{\linkS4class{result_collection}} object.}

\item{value}{A matrix of samples by signature exposures}
}
\value{
A matrix of exposures
}
\description{
The \code{exposure} matrix contains estimated amount of
each signature for each sample. Rows correspond to each signature and
columns correspond to each sample.
}
\examples{
data(res)
exposures(res, "result", "SBS96", "res")
data(res)
exposures(res, "result", "SBS96", "res") <- matrix()
}
