\name{ndr}
\alias{ndr}
\title{Genearlized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\description{
The main function of Generalized Network-based Dimensionality Reduction and Analysis (GNDA).
}
\usage{
ndr(r,covar=FALSE,cor_method=1,cor_type=1,min_R=0,min_comm=2,Gamma=1,null_model_type=4,
mod_mode=6,min_evalue=0,min_communality=0,com_communalities=0,use_rotation=FALSE,
rotation="oblimin",weight=NULL,seed=NULL)

}
\arguments{
  \item{r}{A numeric data frame}
  \item{covar}{If this value is FALSE (default), it finds the correlation matrix from the raw data. If this value is TRUE, it uses the matrix r as a correlation/similarity matrix.}
  \item{cor_method}{Correlation method (optional).
'1' Pearson's correlation (default),
'2' Spearman's correlation, '3' Kendall's correlation,
'4' Distance correlation}
  \item{cor_type}{Correlation type (optional).
'1' Bivariate correlation (default),
'2' partial correlation, '3' semi-partial correlation}
  \item{min_R}{Minimal square correlation between indicators (default: 0).}
  \item{min_comm}{Minimal number of indicators per community (default: 2).}
  \item{Gamma}{Gamma parameter in multiresolution null modell (default: 1).}
  \item{null_model_type}{'1' Differential Newmann-Grivan's null model, '2' The null model is the mean of square correlations between indicators, '3' The null model is the specified minimal square correlation, '4' Newmann-Grivan's modell (default)}
  \item{mod_mode}{Community-based modularity calculation mode: '1' Louvain modularity, '2' Fast-greedy modularity, '3' Leading Eigen modularity, '4' Infomap modularity, '5' Walktrap modularity, '6' Leiden modularity (default)}
  \item{min_evalue}{Minimal eigenvector centrality value (default: 0)}
  \item{min_communality}{Minimal communality value of indicators (default: 0)}
  \item{com_communalities}{Minimal common communalities (default: 0)}
  \item{use_rotation}{FALSE no rotation (default), TRUE the rotation is used.}
  \item{rotation}{"none", "varimax", "quartimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution. "oblimin" is the default, if use_rotation is TRUE.}
  \item{weight}{The weights of columns. The defalt is NULL (no weights).}
  \item{seed}{default seed value (default=NULL, no seed)}
}
\details{
NDA both works on low and high simple size datasets. If min_evalue=min_communality=com_communalities=0 than there is no feature selection.
}
\value{

\item{communality}{Communality estimates for each item. These are merely the sum of squared factor loadings for that item.  It can be interpreted in correlation matrices.}
\item{loadings}{A standard loading matrix of class “loadings".}
\item{uniqueness}{Uniqueness values of indicators.}
\item{factors}{Number of found factors.}
\item{EVCs}{The list eigenvector centrality value of indicators.}
\item{membership}{The membership value of indicators.}
\item{weight}{The weight of indicators.}
\item{scores}{Estimates of the factor scores are reported (if covar=FALSE).}
\item{centers}{Colum mean of unstandardized score values.}
\item{n.obs}{Number of observations specified or found.}
\item{use_rotation}{FALSE no rotation (default), TRUE the rotation is used.}
\item{rotation}{"none", "varimax", "quartimax", "promax", "oblimin", "simplimax", and "cluster" are possible rotations/transformations of the solution. "oblimin" is the default, if use_rotation is TRUE.}
\item{fn}{Factor name: NDA}
\item{seed}{applied seed value (default=NULL, no seed)}
\item{Call}{Callback function}
}

\references{
Kosztyan, Z. T., Kurbucz, M. T., & Katona, A. I. (2022). Network-based dimensionality reduction of high-dimensional, low-sample-size datasets. Knowledge-Based Systems, 109180. doi:10.1016/j.knosys.2022.109180
}

\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
}

\seealso{\code{\link{plot}}, \code{\link{biplot}}, \code{\link{summary}}.}

\examples{

# Dimension reduction without using any hyperparameters

data(swiss)
df<-swiss
p<-ndr(df)
summary(p)
plot(p)
biplot(p)

# Dimension reduction with using hyperparameters
# min_R=0.1 # The mininal square correlation must be grater than 0.1

p<-ndr(df,min_R = 0.1)
summary(p)
plot(p)

# min_evalue=0.1 # Minimal evector centalities must be greater than 0.1

p<-ndr(df,min_evalue = 0.1)
summary(p)
plot(p)

# minimal and common communality value must be greater than 0.25

p<-ndr(df,min_communality = 0.25,
 com_communalities = 0.25)

# Print factor matrix
cor(p$scores)
plot(p)

# Use factor rotation

p<-ndr(df,min_communality = 0.25,
 com_communalities = 0.25,use_rotation=TRUE)

# Print factor matrix
cor(p$scores)
biplot(p)

# Data reduction - clustering
# Distance is Euclidean's distance
# covar=TRUE means only the distance matrix is considered.

q<-ndr(1-normalize(as.matrix(dist(df))),covar=TRUE)
summary(q)
plot(q)

}

\keyword{multivariate}
\keyword{nonparametric}
\concept{dimensionality}
\concept{reduction}
