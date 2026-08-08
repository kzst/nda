\name{summary.nda}
\alias{summary.nda}
\title{Summarize an NDA Model}
\description{Returns a printable summary without producing output during assignment.}
\usage{\method{summary}{nda}(object, digits = getOption("digits"), ...)}
\arguments{
\item{object}{An \code{nda} object.}
\item{digits}{Preferred printing precision.}
\item{...}{Additional summary arguments.}
}
\value{An object of class \code{summary.nda} containing settings, modularity, community sizes, loadings, communalities, uniquenesses, and inference statistics.}
\examples{fit <- ndr(swiss); s <- summary(fit); s}
\keyword{multivariate}
