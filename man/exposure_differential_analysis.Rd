% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_analysis.R
\name{exposure_differential_analysis}
\alias{exposure_differential_analysis}
\title{Compare exposures of annotated samples}
\usage{
exposure_differential_analysis(
  musica,
  model_name,
  annotation,
  modality = "SBS96",
  result_name = "result",
  method = c("wilcox", "kruskal", "glm.nb"),
  group1 = NULL,
  group2 = NULL,
  ...
)
}
\arguments{
\item{musica}{A \code{\linkS4class{musica}} object.}

\item{model_name}{The name of the model.}

\item{annotation}{Column in the sample_annotations table of the
\code{\linkS4class{musica}} object}

\item{modality}{The modality. Must be "SBS96", "DBS78", or "IND83". Default
\code{"SBS96"}.}

\item{result_name}{Name of the result list entry. Default \code{"result"}.}

\item{method}{Any method in \code{c("wilcox", "kruskal", "glm.nb")}
used to perform differential analysis on signature exposures}

\item{group1}{character vector used in the Wilcox test. Elements in
\code{group1} are compared to elements in \code{group2}. This is
required for \code{annotation} with more than 2 levels.}

\item{group2}{character vector used in the Wilcox test. Elements in
\code{group2} are compared to elements in \code{group1}. This is
required for \code{annotation} with more than 2 levels.}

\item{...}{Additional arguments to be passed to the chosen method}
}
\value{
A matrix containing statistics summarizing the analysis dependent
on the chosen method
}
\description{
\code{exposure_differential_analysis} is used to run
differential analysis on the signature exposures of annotated samples.
}
\examples{
data("res_annot")
exposure_differential_analysis(res_annot,
  model_name = "res_annot",
  annotation = "Tumor_Subtypes", method = "wilcox"
)
}
