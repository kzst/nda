\name{data_gen}
\alias{data_gen}

\title{Generate random block matrix for NDA
}

\description{Generate random block matrix for NDA
}

\usage{
data_gen(n,m,nfactors=2,lambda=1)
}

\arguments{
  \item{n}{number of rows}
  \item{m}{number of columns}
  \item{nfactors}{number of blocks (factors, where the default value is 2)}
  \item{lambda}{exponental smoothing, where the default value is 1}
}
\details{
\code{n}, \code{m}, \code{nfactors} must beintegers, and they are not less than 1; lambda should be a positive real number.
}
\value{

  \item{M}{a dataframe of a block matrix}
}
\author{
Prof. Zsolt T. Kosztyan, Department of Quantitative Methods, Institute of Management, Faculty of Business and Economics, University of Pannonia, Hungary

e-mail: kzst@gtk.uni-pannon.hu
}

\examples{
# Specification 30 by 10 random block matrices with 2 blocks/factors
df<-data_gen(30,10)
library(psych)
scree(df)
biplot(ndr(df))
# Specification 40 by 20 random block matrices with 3 blocks/factors
df<-data_gen(40,20,3)
library(psych)
scree(df)
biplot(ndr(df))
plot(ndr(df))

# Specification 50 by 20 random block matrices with 4 blocks/factors
# lambda=0.1
df<-data_gen(50,15,4,0.1)
scree(df)
biplot(ndr(df))
plot(ndr(df))
}

\keyword{array}
\keyword{multivariate}
\concept{random block matrix}

