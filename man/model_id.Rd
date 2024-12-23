% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_result_model.R
\name{model_id}
\alias{model_id}
\alias{model_id,musica-method}
\alias{model_id,result_collection-method}
\alias{model_id,result_model-method}
\alias{model_id<-}
\alias{model_id<-,musica,matrix-method}
\alias{model_id<-,result_collection,matrix-method}
\alias{model_id<-,result_model,matrix-method}
\title{Retrieve model_id from a result_model, result_collection, or musica
object}
\usage{
model_id(x, ...)

\S4method{model_id}{musica}(x, result, modality, model_id)

\S4method{model_id}{result_collection}(x, modality, model_id)

\S4method{model_id}{result_model}(x)

model_id(x, ...) <- value

\S4method{model_id}{musica,matrix}(x, result, modality, model_id) <- value

\S4method{model_id}{result_collection,matrix}(x, modality, model_id) <- value

\S4method{model_id}{result_model,matrix}(x) <- value
}
\arguments{
\item{x}{A \code{\linkS4class{musica}},
\code{\linkS4class{result_collection}}, or \code{\linkS4class{result_model}}
object generated by a mutational discovery or prediction tool.}

\item{...}{Other inputs}

\item{result}{Name of result from result_list to assign the
model_id. Used when \code{result} is a \code{\linkS4class{musica}} object.}

\item{modality}{Modality to assign the model_id. Used when
\code{result} is a \code{\linkS4class{musica}} or
\code{\linkS4class{result_collection}} object.}

\item{model_id}{Model identifier to assign the model_id. Used when
\code{result} is a \code{\linkS4class{musica}} or
\code{\linkS4class{result_collection}} object.}

\item{value}{Model identifier}
}
\value{
The model_id for the model
}
\description{
Model identifier
}
\examples{
data(res)
model_id(res, "result", "SBS96", "res")
}
