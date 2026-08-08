\name{print.ndrlm}
\alias{print.ndrlm}
\title{Print an NDRLM Model}
\description{Prints the reduction mode, regression engine, target metrics, and response-specific fits.}
\usage{\method{print}{ndrlm}(x, digits = getOption("digits"), ...)}
\arguments{\item{x}{An \code{ndrlm} object.}\item{digits}{Printing precision.}\item{...}{Additional arguments passed to the summary printer.}}
\value{\code{x}, invisibly.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
print(ndrlm(Y, X, optimize = FALSE))
}
\keyword{multivariate}
