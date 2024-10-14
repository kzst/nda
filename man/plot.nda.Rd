\name{plot.nda}
\alias{plot.nda}
\title{
Plot function for Generalized Network-based Dimensionality Reduction and Analysis (GNDA)
}
\description{
Plot variable network graph
}
\usage{
\method{plot}{nda}(x, cuts=0.3, interactive=TRUE,edgescale=1.0,labeldist=-1.5,show_weights=FALSE,...)
}
\arguments{
  \item{x}{an object of class 'NDA'.}
  \item{cuts}{minimal square correlation value for an edge in the correlation network graph (default 0.3).}
  \item{interactive}{Plot interactive visNetwork graph or non-interactive igraph plot (default TRUE).}
  \item{edgescale}{Proportion scale value of edge width.}
  \item{labeldist}{Vertex label distance in non-interactive igraph plot (default value =-1.5).}
  \item{show_weights}{Show edge weights (default FALSE)).}
  \item{...}{other graphical parameters.}
}
\references{
Kosztyán, Z. T., Katona, A. I., Kurbucz, M. T., & Lantos, Z. (2024). Generalized network-based dimensionality analysis. Expert Systems with Applications, 238, 121779. <URL:
https://doi.org/10.1016/j.eswa.2023.121779>.

}
\author{
Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona

e-mail*: kosztyan.zsolt@gtk.uni-pannon.hu
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

# GNDA is used for clustering, where the similarity function is the 1-Euclidean distance
# Data is the swiss data

SIM<-1-normalize(as.matrix(dist(swiss)))
q<-ndr(SIM,covar = TRUE)
plot(q,interactive = FALSE)
}
\keyword{multivariate}
\concept{control chart}
