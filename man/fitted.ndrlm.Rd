\name{fitted.ndrlm}
\alias{fitted.ndrlm}
\title{NDRLM Fitted Values}
\description{Extracts fitted values from every response model.}
\usage{\method{fitted}{ndrlm}(object, ...)}
\arguments{\item{object}{An \code{ndrlm} object.}\item{...}{Reserved for method compatibility.}}
\value{A data frame with one column per fitted response.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
fit <- ndrlm(Y, X, optimize = FALSE)
fitted(fit)
}
\keyword{multivariate}
