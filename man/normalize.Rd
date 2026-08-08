\name{normalize}
\alias{normalize}
\title{Fast Min-Max Normalization}
\description{Rescales a numeric matrix or data frame globally, by row, or by column.}
\usage{normalize(x, type = c("all", "row", "col"), range = c(0, 1),
          na.rm = TRUE)}
\arguments{
\item{x}{Numeric matrix or data frame.}
\item{type}{Normalize all values together, each row, or each column.}
\item{range}{Two increasing finite output limits.}
\item{na.rm}{Remove missing values when finding limits.}
}
\details{A constant slice is mapped to the midpoint of \code{range}.}
\value{A normalized object with the same matrix or data-frame form as \code{x}.}
\examples{
m <- matrix(rnorm(20), 5, 4)
normalize(m)
normalize(m, type = "row", range = c(-1, 1))
}
\keyword{matrix}
