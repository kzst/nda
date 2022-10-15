\name{fs.dimred}
\alias{fs.dimred}
\title{Feature selection for PCA and FA
}
\description{
This function drops variables, which have low communality values and/or are common indicators (i.e. correlates more than one latent variable).
}

\usage{
fs.dimred(fn,DF,min_comm=0.25,com_comm=0.25)
}
\arguments{
  \item{fn}{List. Output of a principal (PCA) or an fa (FA) functions.}
  \item{DF}{Numeric data frame, or a numeric matrix of the data table}
  \item{min_comm}{Scalar between 0 to 1. Minimal communality value, which a variable has to be acheved. Default value is 0.25.}
  \item{com_comm}{Scalar between 0 to 1. Minimal difference value between loadings. Default value is 0.25.}
}
\details{
This function only works with principal and fa functions.

This function drops each variable that has a low communality value (under than min_comm value). In other words, that variable does not fit enough of any latent variable.

This function also drops so-called common indicators, which correlate highly with more than one latent variable. And the difference in the correlation is either lower than the com_comm value or the greatest absolute factor loading value is not twice greater than the second greatest factor loading.
}


\value{

\item{dropped_low}{Numeric data frame or numeric matrix. Set of indicators (i.e. variables), which are dropped by their low communalities. This value is NULL if a correlation matrix is used as an input or there is no dropped indicator.}
\item{dropped_com}{Numeric data frame or numeric matrix. Set of dropped common indicators (i.e. common variables). This value is NULL if a correlation matrix is used as an input or there is no dropped indicator.}
\item{remain_DF}{Numeric data frame or numeric matrix. Set of retained indicators}
\item{...}{Other outputs came from \code{\link{principal}}, \code{\link{fa}}, or in \code{\link{ndr}}}
}

\references{
Abonyi, J., Czvetkó, T., Kosztyán, Z. T., & Héberger, K. (2022). Factor analysis, sparse PCA, and Sum of Ranking Differences-based improvements of the Promethee-GAIA multicriteria decision support technique. Plos one, 17(2), e0264277. doi:10.1371/journal.pone.0264277
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}

\seealso{\code{\link{principal}},\code{\link{fa}},\code{\link{ndr}}.}

\examples{

data<-I40_2020

library(psych)

# Principal Component Analysis (PCA)

pca<-principal(data,nfactors=2,covar=TRUE)
pca

# Feature selection with default values

PCA<-fs.dimred(pca,data)
PCA

# List of dropped, low communality value indicators
print(colnames(PCA$dropped_low))

# List of dropped, common communality value indicators
print(colnames(PCA$dropped_com))

# List of retained indicators
print(colnames(PCA$retained_DF))


# Principal Component Analysis (PCA) of correlation matrix

pca<-principal(cor(data,method="spearman"),nfactors=2,covar=TRUE)
pca

# Feature selection
min_comm<-0.25 # Minimal communality value
com_comm<-0.20 # Minimal common communality value

PCA<-fs.dimred(pca,cor(data,method="spearman"),min_comm,com_comm)
PCA

}

\keyword{multivariate}
\concept{dimensionality}
\concept{reduction}
