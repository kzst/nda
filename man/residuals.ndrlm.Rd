\name{residuals.ndrlm}
\alias{residuals.ndrlm}
\title{NDRLM Residuals}
\description{Extracts observed-minus-fitted residuals from every response model.}
\usage{\method{residuals}{ndrlm}(object, ...)}
\arguments{\item{object}{An \code{ndrlm} object.}\item{...}{Reserved for method compatibility.}}
\value{A data frame with one residual column per fitted response.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(y = freeny.y)
fit <- ndrlm(Y, X, optimize = FALSE)
residuals(fit)
}
\keyword{multivariate}
