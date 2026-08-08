\name{ndr}
\alias{ndr}
\title{Generalized Network-Based Dimensionality Reduction and Analysis}
\description{Builds a squared-association network, detects indicator communities, and estimates one centrality-weighted latent variable per retained community.}
\usage{
ndr(r, covar = FALSE, cor_method = 1, cor_type = 1, min_R = 0,
    min_comm = 2, Gamma = 1, null_model_type = 4, mod_mode = 6,
    min_evalue = 0, min_communality = 0, com_communalities = 0,
    use_rotation = FALSE, rotation = "oblimin", weight = NULL,
    seed = NULL, parallel = FALSE, cores = 1L, sparsecalc = FALSE,
    membership = NULL,
    centrality = c("eigenvector", "authority", "hub", "pagerank"),
    standardized = TRUE, P_mat = NULL, test = FALSE,
    p_adjust = "BH", alpha = NULL, bootstrap = FALSE,
    n_boot = 200L, conf = 0.95)
}
\arguments{
\item{r}{Numeric row-level data, or a square association matrix when \code{covar = TRUE}. Sparse \pkg{Matrix} inputs are accepted.}
\item{covar}{Use \code{r} as an association matrix instead of estimating associations from rows.}
\item{cor_method}{Association method. Numeric values 1--6 retain compatibility; textual values are \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, \code{"distance"}, \code{"gaussian_rank"}, and \code{"biweight"}. Gaussian-rank correlation is a fast monotone robust alternative; biweight correlation limits outlier influence.}
\item{cor_type}{Association type: \code{"full"}, \code{"partial"}, or \code{"semipartial"}; 1--3 are accepted aliases.}
\item{min_R}{Minimum squared association retained as an edge.}
\item{min_comm}{Minimum preferred community size. A singleton is retained when it is the only feasible result.}
\item{Gamma}{Positive resolution multiplier used by community detection and null-model subtraction.}
\item{null_model_type}{Null model: 1 configuration-model expectation, 2 mean positive similarity, 3 \code{min_R}, or 4 no matrix modification (default).}
\item{mod_mode}{Community method: Louvain, fast greedy, leading eigenvector, Infomap, Walktrap, or Leiden; either names or values 1--6.}
\item{min_evalue}{Minimum selected centrality weight for periphery filtering.}
\item{min_communality}{Minimum indicator communality.}
\item{com_communalities}{Minimum separation between the two largest squared loadings for common-indicator filtering. Zero disables this filter.}
\item{use_rotation}{Rotate a multi-factor solution.}
\item{rotation}{Rotation understood by \code{psych::principal}.}
\item{weight}{One finite observation-score weight per indicator.}
\item{seed}{Optional random seed.}
\item{parallel}{Use multicore calculation for pairwise distance associations and bootstrap replicates where supported.}
\item{cores}{Number of requested cores.}
\item{sparsecalc}{Store and process the thresholded similarity network as a sparse matrix.}
\item{membership}{Optional supplied community membership, one value per indicator. Zero excludes an indicator.}
\item{centrality}{Latent-variable weighting: eigenvector, authority, hub, or PageRank centrality.}
\item{standardized}{Standardize input indicators using training means and standard deviations before weighting.}
\item{P_mat}{Optional finite expected-edge matrix. When supplied it overrides \code{null_model_type}.}
\item{test}{Calculate association p-values.}
\item{p_adjust}{Method passed to \code{p.adjust} for false-discovery-rate control.}
\item{alpha}{Optional adjusted-p threshold; nonsignificant associations are set to zero.}
\item{bootstrap}{Bootstrap loading intervals from row-level data.}
\item{n_boot}{Number of bootstrap replicates.}
\item{conf}{Bootstrap confidence level.}
}
\details{The method follows three stages: construct a squared-association graph, detect network communities, and form standardized centrality-weighted latent scores. Partial and semipartial networks can be asymmetric. Null-model types 1--3 operate on nonnegative residual similarities; type 4 leaves the similarity matrix unchanged.}
\value{An object of class \code{nda}. Important elements include \code{scores}, \code{loadings}, \code{communality}, \code{uniqueness}, \code{membership}, \code{EVCs}, \code{R}, \code{correlation}, \code{modularity}, \code{reconstruction}, and \code{stats}. The \code{stats} list contains indicator p-values, optional bootstrap intervals, retained and dropped indicators, community sizes, and null-model metadata. Scores are unavailable for association-matrix input, but reconstruction-consistent loadings, communalities, and a loading-based reconstructed association matrix are returned.}
\references{
Kosztyan, Z. T., Kurbucz, M. T., and Katona, A. I. (2022). Network-based dimensionality reduction of high-dimensional, low-sample-size datasets. \emph{Knowledge-Based Systems}, 109180. \doi{10.1016/j.knosys.2022.109180}.

Kosztyan, Z. T., Katona, A. I., Kurbucz, M. T., and Lantos, Z. (2024). Generalized network-based dimensionality analysis. \emph{Expert Systems with Applications}, 238, 121779. \doi{10.1016/j.eswa.2023.121779}.
}
\examples{
fit <- ndr(swiss, cor_method = "spearman", centrality = "pagerank")
summary(fit)
scores <- predict(fit, swiss[1:3, ])

association <- cor(swiss)
fit_matrix <- ndr(association, covar = TRUE, membership = c(1, 1, 1, 2, 2, 2))
fit_matrix$loadings

\dontrun{
fit_boot <- ndr(swiss, cor_method = "distance", test = TRUE,
                bootstrap = TRUE, n_boot = 500, parallel = TRUE, cores = 4)
}
}
\keyword{multivariate}
