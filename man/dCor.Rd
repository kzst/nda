\name{dCor}
\alias{dCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculating distance correlation of two vectors or columns of a matrix
}

\description{Calculating distance correlation of two vectors or columns of a matrix for NDA.

The calculation is very slow for large matrices!
}

\usage{
dCor(x,y=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, a numeric matrix (in this case y=NULL), or a numeric data frame (in this case y=NULL)}
  \item{y}{a numeric vector (optional)}
}
\details{
If \code{x} is a numeric vector, \code{y} must be specified. If x is a numeric matrix  or numeric data frame, y must be ignored from the parameters.
}
\value{

Either a distance correlation value of vectors \code{x} and \code{y}, or a distance correlation matrix of \code{x}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~

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
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{correlation matrix}
\concept{distance correlation}
% Use only one concept per line.

