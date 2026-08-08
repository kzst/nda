\name{print.nda}
\alias{print.nda}
\title{Print an NDA Model}
\description{Prints the NDA call, settings, modularity, community sizes, loadings, communalities, and dropped indicators.}
\usage{\method{print}{nda}(x, digits = getOption("digits"), ...)}
\arguments{
\item{x}{An \code{nda} object.}
\item{digits}{Printing precision.}
\item{...}{Additional arguments passed to the summary printer.}
}
\value{\code{x}, invisibly.}
\examples{fit <- ndr(swiss); print(fit)}
\keyword{multivariate}
