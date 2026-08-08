\name{ndrlm}
\alias{ndrlm}
\title{Network-Based Dimensionality Reduction and Regression}
\description{Fits supervised models to original or NDA-reduced input and output variables, with optional NSGA-II tuning.}
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
\item{latents}{Use latent variables for inputs (\code{"in"}), outputs (\code{"out"}), both, or neither.}
\item{dircon}{Retain NDA-excluded variables as direct predictors or responses.}
\item{optimize}{Tune indicator weights and four NDA thresholds with NSGA-II.}
\item{target}{One of adjusted R-squared, R-squared, MAE, MAPE, MASE, MSE, RMSE, AIC, or BIC.}
\item{rel_weight}{Constrain and normalize indicator weights to be nonnegative relative weights.}
\item{cor_method,cor_type,min_comm,Gamma,null_model_type,mod_mode,use_rotation,rotation}{Arguments forwarded to \code{ndr}.}
\item{pareto}{Keep one objective per original dependent variable during NSGA-II optimization.}
\item{fit_weights}{Optional weights for aggregating response-specific objectives when their number matches the fitted responses.}
\item{lower.bounds.x,upper.bounds.x}{Lower and upper NSGA-II bounds for input-indicator weights.}
\item{lower.bounds.latentx,upper.bounds.latentx}{Bounds for input-side \code{min_evalue}, \code{min_communality}, \code{com_communalities}, and \code{min_R}.}
\item{lower.bounds.y,upper.bounds.y}{Lower and upper NSGA-II bounds for output-indicator weights.}
\item{lower.bounds.latenty,upper.bounds.latenty}{Bounds for the four output-side NDA thresholds.}
\item{popsize,generations,cprob,cdist,mprob,mdist}{NSGA-II population, generation, crossover, and mutation settings passed to \code{mco::nsga2}.}
\item{seed}{Optional seed; repeated optimization runs use consecutive seeds.}
\item{regression}{Regression engine: ordinary least squares, ridge, LOESS, or elastic net.}
\item{model}{Backward-friendly alias that overrides \code{regression}.}
\item{alpha}{Elastic-net mixing value in [0, 1]. Ridge always uses zero.}
\item{lambda}{Optional penalization value. If omitted, cross-validation selects lambda.}
\item{loess_span}{LOESS smoothing span.}
\item{restarts}{Number of independent NSGA-II runs; the best feasible candidate is refitted.}
\item{parallel,cores,sparsecalc,centrality,standardized}{Performance and latent-score settings forwarded to \code{ndr}.}
\item{membership.X,membership.Y}{Optional fixed NDA community memberships for the input and output sides.}
}
\details{NDRLM is predictive and exploratory; its paths do not by themselves establish causal effects. Ridge and elastic net use \pkg{glmnet}. LOESS is intended for small predictor sets and does not provide linear-path significance tests. AIC and BIC for non-likelihood engines are residual-based approximations.}
\value{An \code{ndrlm} object containing fitted regression models, response-specific metrics, fitted NDA objects, optimized hyperparameters, training designs, and NSGA-II results.}
\references{Kosztyan, Z. T., Kurbucz, M. T., Katona, A. I., and Khan, Z. Network-based dimensionality reduction and structural modeling.}
\examples{
X <- as.data.frame(freeny.x)
Y <- data.frame(income = freeny.y)
fit <- ndrlm(Y, X, optimize = FALSE, regression = "lm")
summary(fit)
predict(fit, X[1:3, ])

\dontrun{
fit_en <- ndrlm(Y, X, regression = "elastic_net", alpha = 0.5,
                popsize = 20, generations = 50, restarts = 3)
}
}
\keyword{multivariate}
