\name{pdCor}
\alias{pdCor}
\title{Calculating partial distance correlation of columns of a matrix
}

\description{Calculating partial distance correlation of two columns of a matrix for Network-based Dimensionality Reduction and Analysis (NDA).

The calculation is very slow for large matrices!
}

\usage{
pdCor(x)
}
\arguments{
  \item{x}{a a numeric matrix, or a numeric data frame}
}

\value{

Partial distance correlation matrix of \code{x}.
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
# Specification of partial distance correlaction matrix.
x<-matrix(rnorm(36),nrow=6)
pdCor(x)
}

\keyword{array}
\keyword{multivariate}
\concept{correlation matrix}
\concept{distance correlation}

