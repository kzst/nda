\name{plot.nda}
\alias{plot.nda}
\title{Plot an NDA Association Network}
\description{Creates a printable static \pkg{igraph} or interactive \pkg{visNetwork} network object. Assignment does not draw; printing the result draws it.}
\usage{
\method{plot}{nda}(x, cuts = 0.3, interactive = FALSE, edgescale = 1,
    labeldist = -1.5, show_weights = FALSE, shapes = NULL,
    shape = "automatic", edge.color.links = TRUE, edge.colors = NULL,
    positive.color = "#009E73", negative.color = "#C62828",
    neutral.color = "grey65", positive.lty = 1, negative.lty = 1,
    neutral.lty = 1, edge.width.range = c(1, 6),
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
\item{shape}{\code{"automatic"} uses one common rectangle for every indicator.
Alternatively, supply one shape or one shape per node. Community membership is
shown by colour, not shape.}
\item{shapes}{Backward-compatible alias for \code{shape}.}
\item{edge.color.links}{Colour links by the sign of the unsquared association.}
\item{edge.colors}{Optional legacy low/middle/high colour scale. When
\code{NULL}, positive, negative, and zero associations use the three explicit
colours below.}
\item{positive.color,negative.color,neutral.color}{Default signed edge colours.}
\item{positive.lty,negative.lty,neutral.lty}{Line types for signed edges. All
are solid by default.}
\item{edge.width.range}{Minimum and maximum edge widths. Absolute association
magnitude is rescaled to this interval.}
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
