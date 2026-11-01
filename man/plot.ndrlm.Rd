\name{plot.ndrlm}
\alias{plot.ndrlm}
\title{Plot an NDRLM Path Network}
\description{Creates a layered, printable network of indicators, latent variables, and fitted regression links. A latent variable and its indicators share one community colour.}
\usage{
\method{plot}{ndrlm}(x, sig = 0.05, interactive = FALSE,
    show_edge_labels = TRUE, edge.label.size = 9,
    vertex.label.size = 11, edgescale = 1, shape = "automatic",
    positive.color = "#009E73", negative.color = "#C62828",
    loading.color = "grey60", positive.lty = 1, negative.lty = 1,
    loading.lty = 1, edge.width.range = c(1, 6),
    igraph_args = list(), visnetwork_args = list(), ...)
}
\arguments{
\item{x}{An \code{ndrlm} object.}
\item{sig}{Significance threshold for ordinary least-squares paths. Penalized and LOESS paths do not have p-values and are not filtered by \code{sig}.}
\item{interactive}{Create a \pkg{visNetwork} widget.}
\item{show_edge_labels}{Show loading and path estimates.}
\item{edge.label.size,vertex.label.size}{Label sizes.}
\item{edgescale}{Absolute-weight edge-width multiplier.}
\item{shape}{\code{"automatic"} uses rectangles for observed variables and
circles for latent variables. Community membership is encoded by colour.
Alternatively, supply one shape or one shape per node.}
\item{positive.color,negative.color,loading.color}{Colours for positive paths,
negative paths, and measurement loadings.}
\item{positive.lty,negative.lty,loading.lty}{Line types for the three edge
roles. All are solid by default.}
\item{edge.width.range}{Minimum and maximum widths to which absolute estimates
are rescaled separately within measurement loadings and structural paths.}
\item{igraph_args}{Named static plotting arguments.}
\item{visnetwork_args}{Named \code{network}, \code{options}, and \code{physics} argument lists.}
\item{...}{Additional static \code{plot.igraph} arguments.}
}
\value{An \code{ndrlm_plot} object; printing it renders the graph.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
p <- plot(ndrlm(Y, X, optimize = FALSE))
p
}
\keyword{multivariate}
