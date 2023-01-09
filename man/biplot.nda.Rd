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


\seealso{\code{\link{plot}}, \code{\link{summary}}, \code{\link{ndr}}.}

\examples{
# Biplot function without feature selection

data("CrimesUSA1990.X")
df<-CrimesUSA1990.X
p<-ndr(df)
biplot(p)

# Biplot function with feature selection
# minimal eigen values (min_evalue) is 0.0065
# minimal communality value (min_communality) is 0.1
# minimal common communality value (com_communalities) is 0.1

p<-ndr(df,min_evalue = 0.0065,min_communality = 0.1,com_communalities = 0.1)
biplot(p)
}
\keyword{plot}
