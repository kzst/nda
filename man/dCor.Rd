\name{dCor}
\alias{dCor}
\title{Distance-Correlation Matrix}
\description{Calculates distance correlation for two vectors or for every unique pair of matrix columns.}
\usage{dCor(x, y = NULL, test = FALSE, adjust = "BH", alpha = NULL,
     parallel = FALSE, cores = 1L, R = 499L)}
\arguments{
\item{x}{Numeric vector, matrix, or data frame.}
\item{y}{Optional numeric vector paired with \code{x}.}
\item{test}{Run permutation tests and return adjusted p-values.}
\item{adjust}{Method passed to \code{p.adjust}.}
\item{alpha}{Optional adjusted-p threshold used to zero nonsignificant coefficients.}
\item{parallel}{Calculate unique matrix pairs in parallel on supported systems.}
\item{cores}{Requested core count.}
\item{R}{Permutation count passed to \code{energy::dcor.test}.}
}
\value{A coefficient or symmetric matrix when \code{test = FALSE}; otherwise a list with \code{r} and \code{p}.}
\examples{
x <- matrix(rnorm(80), 20, 4)
dCor(x)
\dontrun{dCor(x, test = TRUE, R = 999)}
}
\keyword{multivariate}
