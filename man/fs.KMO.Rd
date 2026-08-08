\name{fs.KMO}
\alias{fs.KMO}
\title{KMO-Based Feature Selection}
\description{Iteratively removes the indicator with the lowest measure of sampling adequacy until the requested overall KMO is reached.}
\usage{fs.KMO(data, min_MSA = 0.5, cor.mtx = FALSE,
       max_iter = ncol(data) - 1L, ...)}
\arguments{
\item{data}{Numeric row-level data or a correlation matrix.}
\item{min_MSA}{Target overall Kaiser-Meyer-Olkin value in [0, 1].}
\item{cor.mtx}{Treat \code{data} as a symmetric correlation matrix.}
\item{max_iter}{Maximum number of removals.}
\item{...}{Reserved for backward compatibility.}
}
\value{Filtered data or correlation matrix. The final KMO result and removed variable names are stored in \code{"KMO"} and \code{"removed"} attributes.}
\examples{filtered <- fs.KMO(swiss, min_MSA = 0.6)}
\keyword{multivariate}
