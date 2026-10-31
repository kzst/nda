\name{fs.dimred}
\alias{fs.dimred}
\title{Loading-Based Feature Selection}
\description{Selects indicators from NDA, PCA, factor-analysis, and other objects that expose loadings, rotations, or variable coordinates. When the original model call is available, the model is refitted after every removal so that all later decisions use the current communalities and loadings.}
\usage{fs.dimred(fn, DF, min_comm = 0.25, com_comm = 0.25, ...)}
\arguments{
\item{fn}{Fitted dimensionality-reduction object. Supported structures include \code{nda}, \code{prcomp}, \code{princomp}, \pkg{psych} PCA/FA, and objects with \code{loadings}, \code{rotation}, or \code{var$coord}.}
\item{DF}{Original numeric data or a symmetric correlation matrix.}
\item{min_comm}{Minimum indicator communality.}
\item{com_comm}{Minimum separation between the two largest absolute loadings. An indicator is treated as common only when its largest loading is both less than twice and less than \code{com_comm} above its second-largest loading.}
\item{...}{Reserved for future model-specific extraction.}
}
\value{The final refitted model, preserving its original class, with \code{dropped_low}, \code{dropped_com}, and \code{retained_DF} (also available as \code{remain_DF}) added. The original fitted model is available as \code{model}.}
\examples{
pca <- prcomp(swiss, scale. = TRUE)
selected <- fs.dimred(pca, swiss)
colnames(selected$retained_DF)
}
\keyword{multivariate}
