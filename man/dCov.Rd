\name{dCov}
\alias{dCov}
\title{Calculating distance covariance of two vectors or columns of a matrix
}

\description{Calculating distance covariance of two vectors or columns of a matrix for Generalized Network-based Dimensionality Reduction and Analysis (GNDA).

The calculation is very slow for large matrices!
}

\usage{
dCov(x,y=NULL)
}
\arguments{
  \item{x}{a numeric vector, matrix or data frame.}
  \item{y}{NULL (default) or a vector, matrix or data frame with compatible dimensions to x. The default is equivalent to y = x (but more efficient).}
}
\details{
If \code{x} is a numeric vector, \code{y} must be specified. If \code{x} is a numeric matrix  or numeric data frame, y will be neglected.
}
\value{

Either a distance covariance value of vectors \code{x} and \code{y}, or a distance covariance matrix of \code{x} if \code{x} is a matrix or a dataframe.
}
\references{
Rizzo M, Szekely G (2021). _energy: E-Statistics:
Multivariate Inference via the Energy of Data_. R
package version 1.7-8, <URL:
https://CRAN.R-project.org/package=energy>.

}
\author{
Prof. Zsolt T. Kosztyan, Department of Quantitative Methods, Institute of Management, Faculty of Business and Economics, University of Pannonia, Hungary

e-mail: kosztyan.zsolt@gtk.uni-pannon.hu
}

\examples{
# Specification of distance covariance value of vectors x and y.
x<-rnorm(36)
y<-rnorm(36)
dCov(x,y)
# Specification of distance covariance matrix.
x<-matrix(rnorm(36),nrow=6)
dCov(x)
}

\keyword{array}
\keyword{multivariate}
\concept{correlation matrix}
\concept{distance correlation}

