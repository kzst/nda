\name{fs.dimred}
\alias{fs.dimred}
\title{Loading-Based Feature Selection}
\description{Selects indicators from NDA, PCA, factor-analysis, and other objects that expose loadings, rotations, or variable coordinates.}
\usage{fs.dimred(fn, DF, min_comm = 0.25, com_comm = 0.25, ...)}
\arguments{
\item{fn}{Fitted dimensionality-reduction object. Supported structures include \code{nda}, \code{prcomp}, \code{princomp}, \pkg{psych} PCA/FA, and objects with \code{loadings}, \code{rotation}, or \code{var$coord}.}
\item{DF}{Original numeric data or a symmetric correlation matrix.}
\item{min_comm}{Minimum indicator communality.}
\item{com_comm}{Minimum separation between the two largest squared loadings.}
\item{...}{Reserved for future model-specific extraction.}
}
\value{A list containing \code{dropped_low}, \code{dropped_com}, \code{retained_DF} (also available as \code{remain_DF}), communalities, loadings, and the input model.}
\examples{
pca <- prcomp(swiss, scale. = TRUE)
selected <- fs.dimred(pca, swiss)
colnames(selected$retained_DF)
}
\keyword{multivariate}
