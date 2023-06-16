\name{biplot.nda}
\alias{biplot.nda}
\title{
Biplot function for Generalized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\description{
Biplot function for Generalized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\usage{
\method{biplot}{nda}(x, main=NULL,...)
}
\arguments{
  \item{x}{an object of class 'NDA'.}
  \item{main}{main title of biplot.}
  \item{...}{other graphical parameters.}
}


\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}


\seealso{\code{\link{plot}}, \code{\link{summary}}, \code{\link{ndr}}, \code{\link{data_gen}}.}

\examples{
# Biplot function without feature selection

# Generate 200 x 50 random block matrix with 3 blocks and lambda=0 parameter

df<-data_gen(200,50,3,0)
p<-ndr(df)
biplot(p)
}
\keyword{plot}
