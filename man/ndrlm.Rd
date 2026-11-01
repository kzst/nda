\name{ndrlm}
\alias{ndrlm}
\title{Network-Based Dimensionality Reduction and Regression}
\description{
Fits supervised models to original or NDA-reduced input and output variables,
with optional NSGA-II tuning.
}
\usage{
ndrlm(Y, X, latents = "in", dircon = FALSE, optimize = TRUE,
      target = "adj.r.square", rel_weight = FALSE, cor_method = 1,
      cor_type = 1, min_comm = 2, Gamma = 1, null_model_type = 4,
      mod_mode = 1, use_rotation = FALSE, rotation = "oblimin",
      pareto = FALSE, fit_weights = NULL,
      lower.bounds.x = rep(-100, ncol(X)),
      upper.bounds.x = rep(100, ncol(X)),
      lower.bounds.latentx = c(0, 0, 0, 0),
      upper.bounds.latentx = c(0.6, 0.6, 0.6, 0.3),
      lower.bounds.y = rep(-100, ncol(Y)),
      upper.bounds.y = rep(100, ncol(Y)),
      lower.bounds.latenty = c(0, 0, 0, 0),
      upper.bounds.latenty = c(0.6, 0.6, 0.6, 0.3),
      popsize = 20, generations = 30, cprob = 0.7, cdist = 5,
      mprob = 0.2, mdist = 10, seed = NULL,
      regression = c("lm", "ridge", "loess", "elastic_net"),
      model = NULL, alpha = 0.5, lambda = NULL, loess_span = 0.75,
      restarts = 1L, parallel = FALSE, cores = 1L,
      sparsecalc = FALSE, centrality = "eigenvector",
      standardized = TRUE, membership.X = NULL, membership.Y = NULL)
}
\arguments{
\item{Y}{Numeric dependent-variable data.}
\item{X}{Numeric independent-variable data with the same rows as \code{Y}.}
\item{latents}{Use latent variables for inputs (\code{"in"}), outputs
(\code{"out"}), both, or neither.}
\item{dircon}{Retain NDA-excluded variables as direct predictors or responses.}
\item{optimize}{Tune indicator weights and four NDA thresholds with NSGA-II.}
\item{target}{One of adjusted R-squared, R-squared, MAE, MAPE, MASE, MSE,
RMSE, AIC, or BIC.}
\item{rel_weight}{Constrain and normalize indicator weights to be nonnegative
relative weights.}
\item{cor_method,cor_type,min_comm,Gamma,null_model_type,mod_mode,use_rotation,rotation}{Arguments
forwarded to \code{ndr}.}
\item{pareto}{Keep one objective per original dependent variable during
NSGA-II optimization.}
\item{fit_weights}{Optional weights for aggregating response-specific
objectives when their number matches the fitted responses.}
\item{lower.bounds.x,upper.bounds.x}{Lower and upper NSGA-II bounds for
input-indicator weights.}
\item{lower.bounds.latentx,upper.bounds.latentx}{Bounds for input-side
\code{min_evalue}, \code{min_communality}, \code{com_communalities}, and
\code{min_R}.}
\item{lower.bounds.y,upper.bounds.y}{Lower and upper NSGA-II bounds for
output-indicator weights.}
\item{lower.bounds.latenty,upper.bounds.latenty}{Bounds for the four
output-side NDA thresholds.}
\item{popsize,generations,cprob,cdist,mprob,mdist}{NSGA-II population,
generation, crossover, and mutation settings passed to \code{mco::nsga2}.}
\item{seed}{Optional seed; repeated optimization runs use consecutive seeds.}
\item{regression}{Regression engine: ordinary least squares (\code{"lm"}),
ridge (\code{"ridge"}), LOESS (\code{"loess"}), or elastic net
(\code{"elastic_net"}).}
\item{model}{Backward-friendly alias that overrides \code{regression}.}
\item{alpha}{Elastic-net mixing value in [0, 1]. Ridge always uses zero.}
\item{lambda}{Optional penalization value. If omitted, cross-validation
selects lambda.}
\item{loess_span}{LOESS smoothing span.}
\item{restarts}{Number of independent NSGA-II runs; the best feasible
candidate is refitted.}
\item{parallel,cores,sparsecalc,centrality,standardized}{Performance and
latent-score settings forwarded to \code{ndr}.}
\item{membership.X,membership.Y}{Optional fixed NDA community memberships
for the input and output sides.}
}
\details{
NDRLM is predictive and exploratory; its paths do not by themselves establish
causal effects. Ridge and elastic net use \pkg{glmnet}. LOESS is intended for
small predictor sets and does not provide linear-path significance tests. AIC
and BIC for non-likelihood engines are residual-based approximations.

The class and structure of each element of \code{fits} depend on the selected
regression engine. With \code{regression = "lm"}, each \code{fits[[i]]} is a
genuine \code{lm} object and can be passed directly to \code{summary},
\code{coef}, \code{AIC}, \code{BIC}, \code{predict.lm}, and
\code{lm.beta::lm.beta}. The \pkg{lm.beta} package is meaningful only for
ordinary least-squares fits; it must not be applied to ridge, elastic-net, or
LOESS fits.

With \code{regression = "ridge"}, \code{"elastic_net"}, or \code{"loess"},
each \code{fits[[i]]} is an NDRLM engine wrapper. It contains the engine name in
\code{method}, the underlying fitted object in \code{fit}, and common
fitted-value, residual, coefficient, and predictor metadata. Penalized fits
additionally store the selected \code{lambda} and effective \code{alpha}. For
engine-independent code, use \code{summary}, \code{predict}, \code{fitted},
\code{residuals}, and \code{plot} on the enclosing \code{ndrlm} object instead
of relying on these internal fields.

The NDRLM methods retain backward compatibility with both legacy objects that
stored ordinary least-squares fits directly and early 0.3.0 objects that
wrapped ordinary least-squares fits. Newly fitted ordinary least-squares models
use direct \code{lm} objects, while the newer non-\code{lm} engines retain their
wrappers.
}
\value{
An \code{ndrlm} object containing response-specific metrics, fitted NDA
objects, optimized hyperparameters, training designs, NSGA-II results and an
\code{NSGA_trace} evaluation history when optimization is enabled, and a
response-named \code{fits} list. Elements of \code{fits} are \code{lm} objects
for ordinary least squares and NDRLM engine wrappers for ridge, elastic net,
and LOESS.
}
\references{
Kosztyan, Z. T., Kurbucz, M. T., Katona, A. I., and Khan, Z. Network-based
dimensionality reduction and structural modeling.
}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(income = freeny.y)

fit_lm <- ndrlm(Y, X, latents = "none", regression = "lm")
inherits(fit_lm$fits[[1]], "lm")
summary(fit_lm$fits[[1]])
predict(fit_lm, X[1:3, ])

if (requireNamespace("lm.beta", quietly = TRUE)) {
  summary(lm.beta::lm.beta(fit_lm$fits[[1]]))
}

if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_ridge <- ndrlm(Y, X, latents = "none", regression = "ridge",
                     lambda = 0.1)
  fit_ridge$fits[[1]]$method
  predict(fit_ridge, X[1:3, ])
}

if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_en <- ndrlm(Y, X, latents = "none", regression = "elastic_net",
                  alpha = 0.5, lambda = 0.1)
  predict(fit_en, X[1:3, ])
}
}
\keyword{multivariate}
