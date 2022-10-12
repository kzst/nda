\name{summary.nda}
\alias{summary.nda}

\title{
Summary function of Network-based Dimensionality Reduction and Analysis (NDA)
}
\description{
Print summary of Network-based, non-parametric dimensionality reduction and analysis (NDA)
}
\usage{
\method{summary}{nda}(object, digits = getOption("digits"), ...)
}

\arguments{
  \item{object}{an object of class 'nda'.}
  \item{digits}{the number of significant digits to use when \code{add.stats = TRUE}.}
  \item{...}{additional arguments affecting the summary produced.}
}

\value{No return value}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kzst@gtk.uni-pannon.hu
}


\seealso{\code{\link{biplot}}, \code{\link{plot}}, \code{\link{ndr}}.}

\examples{
# Example of summary function of NDA without feature selection

data("CrimesUSA1990.X")
df<-CrimesUSA1990.X
p<-ndr(df)
summary(p)

# Example of summary function of NDA with feature selection
# minimal eigen values (min_evalue) is 0.0065
# minimal communality value (min_communality) is 0.1
# minimal common communality value (com_communalities) is 0.1

p<-ndr(df,min_evalue = 0.0065,min_communality = 0.1,com_communalities = 0.1)
summary(p)

}
\keyword{multivariate}
