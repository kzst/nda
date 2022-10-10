\name{plot.nda}
\alias{plot.nda}
\title{
Plot function for Network-based, non-parametric dimensionality reduction and analysis (NDA)
}
\description{
Plot variable network graph
}
\usage{
\method{plot}{nda}(x, cuts=0.3,...)
}
\arguments{
  \item{x}{an object of class 'NDA'.}
  \item{cuts}{minimal correlation square value for an edge in the variable network graph.}
  \item{...}{other graphical parameters.}
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kzst@gtk.uni-pannon.hu
}

\seealso{\code{\link{biplot}}, \code{\link{summary}}, \code{\link{ndr}}.}

\examples{
# Plot function with feature selection

data("CrimesUSA1990.X")
df<-CrimesUSA1990.X
p<-ndr(df)
biplot(p,main="Biplot of CrimesUSA1990 without feature selection")

# Plot function with feature selection
# minimal eigen values (min_evalue) is 0.0065
# minimal communality value (min_communality) is 0.1
# minimal common communality value (com_communalities) is 0.1

p<-ndr(df,min_evalue = 0.0065,min_communality = 0.1,com_communalities = 0.1)

# Plot with default (cuts=0.3)
plot(p)

# Plot with higher cuts
plot(p,cuts=0.6)

}
\keyword{multivariate}
\concept{control chart}
