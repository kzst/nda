\name{dCor}
\alias{dCor}
\title{Calculating distance correlation of two vectors or columns of a matrix
}

\description{Calculating distance correlation of two vectors or columns of a matrix for NDA.

The calculation is very slow for large matrices!
}

\usage{
dCor(x,y=NULL)
}
\arguments{
  \item{x}{a numeric vector, a numeric matrix (in this case y=NULL), or a numeric data frame (in this case y=NULL)}
  \item{y}{a numeric vector (optional)}
}
\details{
If \code{x} is a numeric vector, \code{y} must be specified. If x is a numeric matrix  or numeric data frame, y must be ignored from the parameters.
}
\value{

Either a distance correlation value of vectors \code{x} and \code{y}, or a distance correlation matrix of \code{x}.
}
\references{
Rizzo M, Szekely G (2021). _energy: E-Statistics:
Multivariate Inference via the Energy of Data_. R
package version 1.7-8, <URL:
https://CRAN.R-project.org/package=energy>.

}
\author{
Prof. Zsolt T. Kosztyan, Department of Quantitative Methods, Institute of Management, Faculty of Business and Economics, University of Pannonia, Hungary

e-mail: kzst@gtk.uni-pannon.hu
}

\examples{
# Specification of distance correlation value of vectors x and y.
x<-rnorm(36)
y<-rnorm(36)
dCor(x,y)
# Specification of distance correlaction matrix.
x<-matrix(rnorm(36),nrow=6)
dCor(x)
}

\keyword{array}
\keyword{multivariate}
\concept{correlation matrix}
\concept{distance correlation}

