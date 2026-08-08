\name{summary.ndrlm}
\alias{summary.ndrlm}
\title{Summarize an NDRLM Model}
\description{Returns a printable regression and reduction summary without printing during assignment.}
\usage{\method{summary}{ndrlm}(object, digits = getOption("digits"), ...)}
\arguments{\item{object}{An \code{ndrlm} object.}\item{digits}{Printing precision.}\item{...}{Additional summary arguments.}}
\value{A \code{summary.ndrlm} object containing settings, variables, response metrics, fitted NDA objects, and regression fits.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
s <- summary(ndrlm(Y, X, optimize = FALSE))
s
}
\keyword{multivariate}
