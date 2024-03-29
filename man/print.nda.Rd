\name{print.nda}
\alias{print.nda}

\title{
Print function of Generalized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\description{
Print summary of Generalized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\usage{
\method{print}{nda}(x, digits = getOption("digits"), ...)
}

\arguments{
  \item{x}{an object of class 'nda'.}
  \item{digits}{the number of significant digits to use when \code{add.stats = TRUE}.}
  \item{...}{additional arguments affecting the summary produced.}
}
\references{
Kosztyán, Z. T., Katona, A. I., Kurbucz, M. T., & Lantos, Z. (2024). Generalized network-based dimensionality analysis. Expert Systems with Applications, 238, 121779. <URL:
https://doi.org/10.1016/j.eswa.2023.121779>.

}
\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kzst@gtk.uni-pannon.hu
}


\seealso{\code{\link{biplot}}, \code{\link{plot}}, \code{\link{summary}}, \code{\link{ndr}}.}


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
print(p)

}
\keyword{multivariate}
