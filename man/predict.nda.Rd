\name{predict.nda}
\alias{predict.nda}
\title{Predict NDA Scores}
\description{Applies training feature standardization, indicator weights, retained centralities, latent-score scaling, and rotation to new rows.}
\usage{\method{predict}{nda}(object, newdata, ...)}
\arguments{
\item{object}{An \code{nda} object fitted to row-level data.}
\item{newdata}{Numeric data containing every original training indicator.}
\item{...}{Reserved for method compatibility.}
}
\value{A matrix with one predicted score column per latent variable. If \code{newdata} is \code{NULL}, training scores are returned.}
\examples{fit <- ndr(swiss); predict(fit, swiss[1:4, ])}
\keyword{multivariate}
