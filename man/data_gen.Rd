\name{data_gen}
\alias{data_gen}
\title{Generate Structured NDA Example Data}
\description{Generates a noisy block matrix without constructing intermediate block-diagonal matrices.}
\usage{data_gen(n, m, nfactors = 2, lambda = 1, sparse = FALSE,
         density = 1, seed = NULL)}
\arguments{
\item{n}{Positive number of rows.}
\item{m}{Positive number of columns.}
\item{nfactors}{Positive number of blocks.}
\item{lambda}{Finite exponential smoothing parameter.}
\item{sparse}{Return a sparse \pkg{Matrix} object.}
\item{density}{Proportion in (0, 1] of block entries retained as nonzero.}
\item{seed}{Optional random seed.}
}
\value{A data frame, or a sparse matrix when \code{sparse = TRUE}.}
\examples{
dense <- data_gen(30, 10, nfactors = 2, seed = 1)
sparse <- data_gen(100, 40, nfactors = 4, sparse = TRUE, density = 0.5)
}
\keyword{array}
