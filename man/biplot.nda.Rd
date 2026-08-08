\name{biplot.nda}
\alias{biplot.nda}
\title{Biplot NDA Scores and Loadings}
\description{Creates a printable two-dimensional biplot or a sequence of three-factor projections.}
\usage{
\method{biplot}{nda}(x, main = NULL, dims = c(1, 2), view_3D = FALSE,
    triples = NULL, interactive = FALSE, angle = 45,
    point_cex = 0.7, label_cex = 0.9, ...)
}
\arguments{
\item{x}{An \code{nda} object fitted to row-level data.}
\item{main}{Optional title.}
\item{dims}{Two factors for the two-dimensional view.}
\item{view_3D}{Display three-factor projections.}
\item{triples}{Optional list of factor triples. By default all triples are used.}
\item{interactive}{Use \pkg{plotly} for rotatable 3D views.}
\item{angle}{Static axonometric projection angle in degrees.}
\item{point_cex,label_cex}{Observation and loading-label size.}
\item{...}{Additional two-dimensional base-graphics arguments.}
}
\value{An \code{nda_biplot} object; printing it renders the selected view or views.}
\examples{
fit <- ndr(swiss)
b <- biplot(fit, dims = c(1, 2))
b
}
\keyword{plot}
