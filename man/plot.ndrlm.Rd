\name{plot.ndrlm}
\alias{plot.ndrlm}
\title{Plot an NDRLM Path Network}
\description{Creates a layered, printable network of indicators, latent variables, and fitted regression links.}
\usage{
\method{plot}{ndrlm}(x, sig = 0.05, interactive = FALSE,
    show_edge_labels = TRUE, edge.label.size = 9,
    vertex.label.size = 11, edgescale = 2,
    igraph_args = list(), visnetwork_args = list(), ...)
}
\arguments{
\item{x}{An \code{ndrlm} object.}
\item{sig}{Significance threshold for ordinary least-squares paths. Penalized and LOESS paths do not have p-values and are not filtered by \code{sig}.}
\item{interactive}{Create a \pkg{visNetwork} widget.}
\item{show_edge_labels}{Show loading and path estimates.}
\item{edge.label.size,vertex.label.size}{Label sizes.}
\item{edgescale}{Absolute-weight edge-width multiplier.}
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
