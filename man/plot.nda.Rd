\name{plot.nda}
\alias{plot.nda}
\title{Plot an NDA Association Network}
\description{Creates a printable static \pkg{igraph} or interactive \pkg{visNetwork} network object. Assignment does not draw; printing the result draws it.}
\usage{
\method{plot}{nda}(x, cuts = 0.3, interactive = TRUE, edgescale = 1,
    labeldist = -1.5, show_weights = FALSE, shapes = "circle",
    edge.color.links = TRUE,
    edge.colors = c("darkred", "grey", "darkblue"),
    edge.label.size = 10, edge.label.color = "#333333",
    vertex.label.size = 12, vertex.label.color = "black",
    layout = "nicely", igraph_args = list(), visnetwork_args = list(), ...)
}
\arguments{
\item{x}{An \code{nda} object.}
\item{cuts}{Minimum squared correlation, or a two-value retained range. For association-matrix models a scalar is interpreted relative to the nonzero minimum and maximum.}
\item{interactive}{Create a \pkg{visNetwork} widget instead of a static graph.}
\item{edgescale}{Edge-width multiplier.}
\item{labeldist}{Static vertex-label distance.}
\item{show_weights}{Label edges with unsquared correlations for row-level models and supplied similarities for association-matrix models.}
\item{shapes}{One or more vertex shapes, recycled over indicators.}
\item{edge.color.links}{Color links by signed association or by the association-matrix minimum, median, and maximum.}
\item{edge.colors}{Three-color low/middle/high scale.}
\item{edge.label.size,edge.label.color}{Edge-label appearance.}
\item{vertex.label.size,vertex.label.color}{Vertex-label appearance.}
\item{layout}{An igraph layout suffix such as \code{"nicely"}, \code{"fr"}, or \code{"circle"}.}
\item{igraph_args}{Named static plot arguments.}
\item{visnetwork_args}{Named lists called \code{network}, \code{nodes}, \code{edges}, \code{options}, and \code{layout}, forwarded to matching visNetwork functions.}
\item{...}{Additional static \code{plot.igraph} arguments.}
}
\value{An object of class \code{nda_plot}; print it to render.}
\examples{
fit <- ndr(swiss)
p <- plot(fit, interactive = FALSE, show_weights = TRUE)
p
}
\keyword{multivariate}
