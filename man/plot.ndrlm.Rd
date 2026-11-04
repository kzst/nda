\name{plot.ndrlm}
\alias{plot.ndrlm}
\title{Plot an NDRLM Path Network}
\description{Creates a layered, printable network of indicators, latent
variables, and fitted regression links. A latent variable and its indicators
share one community colour. A light qualitative node palette is used by
default. Structural path labels and widths use standardized coefficients by
default, including ridge and elastic-net fits. Legacy unnamed fits recover
their response and standardization data from the fitted model.}
\usage{
\method{plot}{ndrlm}(x, sig = 0.05, interactive = FALSE,
    show_edge_labels = TRUE, edge.label.size = 9,
    vertex.label.size = 11, edgescale = 1, shape = "automatic",
    standardized.paths = TRUE, layer.positions = NULL,
    node.lighten = 0.18, node.alpha = 0.90,
    static.engine = c("grid", "igraph"),
    positive.color = "#009E73", negative.color = "#C62828",
    loading.color = "grey60", positive.lty = 1, negative.lty = 1,
    loading.lty = 2, edge.width.range = c(1, 6),
    edge.label.color = "#333333",
    edge.label.background = "#FFFFFFE6", edge.label.offset = 0.022,
    node.border.color = "#25313C",
    igraph_args = list(), visnetwork_args = list(), ...)
}
\arguments{
\item{x}{An \code{ndrlm} object.}
\item{sig}{Significance threshold for ordinary least-squares paths. Penalized
and LOESS paths do not have p-values and are not filtered by \code{sig}.}
\item{interactive}{Create a \pkg{visNetwork} widget.}
\item{show_edge_labels}{Show loading estimates and structural path
coefficients. Structural path labels are standardized by default.}
\item{edge.label.size,vertex.label.size}{Label sizes.}
\item{edgescale}{Absolute-weight edge-width multiplier.}
\item{shape}{\code{"automatic"} uses rectangles for observed variables and
circles for latent variables. Community membership is encoded by colour.
Alternatively, supply one shape or one shape per node.}
\item{standardized.paths}{If \code{TRUE}, use standardized structural
coefficients. For linear, ridge, and elastic-net fits these are calculated as
\eqn{\beta_j s_{x_j}/s_y} from the complete training observations. Set to
\code{FALSE} to display raw coefficients.}
\item{layer.positions}{Optional horizontal layer positions. Supply four
strictly increasing values for the input-indicator, input-latent,
output-latent, and output-indicator layers, or one value per active layer.
The default closes gaps left by absent latent layers.}
\item{node.lighten,node.alpha}{Amount of white mixed into node fills and
their opacity. Both accept values between zero and one.}
\item{static.engine}{Static renderer. The default \code{"grid"} renderer uses
a centred layered layout and labels with opaque backgrounds; \code{"igraph"}
retains the legacy renderer.}
\item{positive.color,negative.color,loading.color}{Colours for positive paths,
negative paths, and measurement loadings.}
\item{positive.lty,negative.lty,loading.lty}{Line types for the three edge
roles. Structural paths are solid and measurement links are dashed by default.}
\item{edge.width.range}{Minimum and maximum widths to which absolute estimates
are rescaled separately within measurement loadings and structural paths.
Structural widths use the same standardized or raw coefficients as their
labels.}
\item{edge.label.color,edge.label.background,edge.label.offset}{Edge-label
colour, opaque background, and perpendicular displacement from the edge.}
\item{node.border.color}{Node-border colour.}
\item{igraph_args}{Named legacy static plotting arguments, used when
\code{static.engine = "igraph"}.}
\item{visnetwork_args}{Named \code{network}, \code{nodes}, \code{edges},
\code{options}, \code{interaction}, and \code{physics} argument lists.}
\item{...}{Additional legacy static \code{plot.igraph} arguments.}
}
\value{An \code{ndrlm_plot} object; printing it renders the graph.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
p <- plot(ndrlm(Y, X, optimize = FALSE))
p
}
\keyword{multivariate}


