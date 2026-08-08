\name{dCov}
\alias{dCov}
\title{Distance-Covariance Matrix}
\description{Calculates distance covariance for two vectors or for every unique pair of matrix columns.}
\usage{dCov(x, y = NULL, test = FALSE, adjust = "BH", alpha = NULL,
     parallel = FALSE, cores = 1L, R = 499L)}
\arguments{
\item{x}{Numeric vector, matrix, or data frame.}
\item{y}{Optional numeric vector paired with \code{x}.}
\item{test}{Run permutation tests and return adjusted p-values.}
\item{adjust}{Method passed to \code{p.adjust}.}
\item{alpha}{Optional adjusted-p threshold used to zero nonsignificant off-diagonal values.}
\item{parallel}{Calculate unique matrix pairs in parallel on supported systems.}
\item{cores}{Requested core count.}
\item{R}{Permutation count passed to \code{energy::dcov.test}.}
}
\value{A covariance value or symmetric matrix when \code{test = FALSE}; otherwise a list with \code{r} and \code{p}.}
\examples{x <- matrix(rnorm(80), 20, 4); dCov(x)}
\keyword{multivariate}
