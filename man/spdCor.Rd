\name{spdCor}
\alias{spdCor}
\title{Semipartial Distance-Correlation Matrix}
\description{Forms an asymmetric semipartial distance-correlation matrix from distance covariances.}
\usage{spdCor(x, test = FALSE, adjust = "BH", alpha = NULL,
       parallel = FALSE, cores = 1L)}
\arguments{
\item{x}{Numeric matrix or data frame.}
\item{test}{Return approximate coefficient p-values.}
\item{adjust}{Method passed to \code{p.adjust}.}
\item{alpha}{Optional adjusted-p filtering threshold.}
\item{parallel}{Calculate distance-covariance pairs in parallel where supported.}
\item{cores}{Requested core count.}
}
\details{When the distance-covariance matrix is singular, a Moore-Penrose generalized inverse is used. The optional p-values use a t approximation and should be interpreted as screening diagnostics.}
\value{A matrix, or a list with matrices \code{r} and \code{p} when \code{test = TRUE}.}
\examples{x <- matrix(rnorm(100), 20, 5); spdCor(x)}
\keyword{multivariate}
