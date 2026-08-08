\name{predict.ndrlm}
\alias{predict.ndrlm}
\title{Predict from an NDRLM Model}
\description{Transforms new independent variables with the fitted NDA model and predicts every fitted response.}
\usage{
\method{predict}{ndrlm}(object, newdata = NULL, se.fit = FALSE,
    interval = c("none", "confidence", "prediction"), level = 0.95, ...)
}
\arguments{
\item{object}{An \code{ndrlm} object.}
\item{newdata}{New original independent variables. If omitted, fitted values are returned.}
\item{se.fit}{Request standard errors for ordinary least-squares fits.}
\item{interval}{No interval, a confidence interval, or a prediction interval. Intervals are available for ordinary least squares only.}
\item{level}{Confidence level.}
\item{...}{Additional arguments passed to \code{predict.lm}.}
}
\value{A response-named data frame for ordinary predictions, or a response-named list when standard errors or intervals are requested.}
\examples{
X <- as.data.frame(freeny.x); Y <- data.frame(income = freeny.y)
fit <- ndrlm(Y, X, optimize = FALSE)
predict(fit, X[1:3, ])
}
\keyword{multivariate}
