\name{fs.KMO}
\alias{fs.KMO}
\title{Feature selection for KMO
}
\description{
Drop variables if their MSA_i valus is lower than a threshold, in order to increase the overall KMO (MSA) value.
}
\usage{
fs.KMO(data,min_MSA=0.5,cor.mtx=FALSE)

}
\arguments{
  \item{data}{A numeric data frame}
  \item{min_MSA}{A numeric value. Minimal MSA value for variable i}
  \item{cor.mtx}{Boolean value. The input is either a correlation matrix (cor.mtx=TRUE), or not (cor.mtx=FALSE)}
}
\details{
Low Kaiser-Meyer-Olkin (KMO) measure os sampling adequacy do not suggest to use principal component analysis or factor analysis. Therefore, this function drop variables with low KMO/MSA values.
}
\value{

\item{data}{Cleaned data or the cleaned correlation matrix.}
}

\references{
Abonyi, J., Czvetkó, T., Kosztyán, Z. T., & Héberger, K. (2022). Factor analysis, sparse PCA, and Sum of Ranking Differences-based improvements of the Promethee-GAIA multicriteria decision support technique. Plos one, 17(2), e0264277. doi:10.1371/journal.pone.0264277
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}

\seealso{\code{\link{summary}}.}

\examples{

library(psych)
data(I40_2020)
data<-I40_2020
KMO(fs.KMO(data,min_MSA=0.7,cor.mtx=FALSE))
}

\keyword{multivariate}
\concept{dimensionality}
\concept{reduction}
