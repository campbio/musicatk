% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_exposures.R
\name{plot_exposures}
\alias{plot_exposures}
\title{Display sample exposures with bar, box, or violin plots}
\usage{
plot_exposures(
  result,
  plot_type = c("bar", "box", "violin"),
  proportional = FALSE,
  group_by = "none",
  color_by = c("signature", "annotation"),
  annotation = NULL,
  num_samples = NULL,
  sort_samples = "total",
  threshold = NULL,
  same_scale = FALSE,
  add_points = FALSE,
  point_size = 2,
  label_x_axis = FALSE,
  legend = TRUE,
  plotly = FALSE
)
}
\arguments{
\item{result}{A \code{\linkS4class{musica_result}} object generated by
a mutational discovery or prediction tool.}

\item{plot_type}{One of \code{"bar"}, \code{"box"}, or \code{"violin"}.
Default \code{"bar"}.}

\item{proportional}{If \code{TRUE}, then the exposures will be normalized
to between 0 and 1 by dividing by the total number of counts for each sample.
Default \code{FALSE}.}

\item{group_by}{Determines how to group samples into the subplots
(i.e. facets). One of \code{"none"}, \code{"signature"} or
\code{"annotation"}. If set to \code{"annotation"},
then a sample annotation must be supplied via the 
\code{annotation} parameter. Default \code{"none"}.}

\item{color_by}{Determines how to color the bars or box/violins. One of
\code{"signature"} or \code{"annotation"}. If set to \code{"annotation"},
then a sample annotation must be supplied via the 
\code{annotation} parameter. Default \code{"signature"}.}

\item{annotation}{Sample annotation used to group the subplots and/or
color the bars, boxes, or violins. Default \code{NULL}.}

\item{num_samples}{The top number of sorted samples to display. If
\code{NULL}, then all samples will be displayed. If \code{group_by} is set,
then the top samples will be shown within each group. Default \code{NULL}.}

\item{sort_samples}{This is used to change how samples are sorted in
the barplot from left to right. If set to \code{"total"}, then samples
will be sorted from those with the highest number of mutation counts to the
lowest (regardless of how the parameter \code{"proportional"} is set).
If set to \code{"name"}, then samples are sorted by their name with the
\code{\link[gtools]{mixedsort}} function. If set to 
one or more signature names (e.g. \code{"Signature1"}), then samples will
be sorted from those with the highest level of that signature to the lowest.
If multiple signatures are supplied then, samples will be sorted by each
signature sequentially. Default \code{"total"}.}

\item{threshold}{Exposures less than this threshold will be set to 0.
This is most useful when more than one signature is supplied to 
\code{sort_samples} as samples that are set to zero for the first exposure
will then be sorted by the levels of the second exposure.
Default \code{NULL}.}

\item{same_scale}{If \code{TRUE}, then all subplots will have the 
same scale. Only used when \code{group_by} is set. Default \code{FALSE}.}

\item{add_points}{If \code{TRUE}, then points for individual sample 
exposures will be plotted on top of the violin/box plots. Only used when
\code{plot_type} is set to \code{"violin"} or \code{"box"}. 
Default \code{TRUE}.}

\item{point_size}{Size of the points to be plotted on top of the
violin/box plots. Only used when \code{plot_type} is set to \code{"violin"}
or \code{"box"} and \code{add_points} is set to \code{TRUE}. 
Default \code{2}.}

\item{label_x_axis}{If \code{TRUE}, x-axis labels will be displayed at
the bottom of the plot. Default \code{FALSE}.}

\item{legend}{If \code{TRUE}, the legend will be displayed.
Default \code{TRUE}.}

\item{plotly}{If \code{TRUE}, the the plot will be made interactive
using \code{\link[plotly]{plotly}}. Default \code{FALSE}.}
}
\value{
Generates a ggplot or plotly object
}
\description{
The distributions of mutational signatures can be viewed 
with barplots or box/violin plots. Barplots are most useful for viewing
the proportion of signatures within and across samples. The box/violin plots
are most useful for viewing the distributions of signatures with respect to 
sample annotations. Samples can be grouped using the \code{group_by}
parameter. For barplots, various methods of sorting samples from left
to right can be chosen using the \code{sort_samples} parameter.
}
\examples{
data(res_annot)
plot_exposures(res_annot, plot_type = "bar", annotation = "Tumor_Subtypes")
}
