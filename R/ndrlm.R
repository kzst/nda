#-----------------------------------------------------------------------------#
#                                                                             #
#  GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (GNDA)     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona,      #
#              Zahid Khan                                                     #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: August 2026                                                #
#-----------------------------------------------------------------------------#


# Network-based dimensionality reduction and regression.

.ndrlm_relative_weights <- function(x) {
  x[!is.finite(x)] <- 0
  if (sum(x) <= .Machine$double.eps) x[] <- 1
  x * length(x) / sum(x)
}

.ndrlm_unpack <- function(hyperparams, latents, nx, ny, rel_weight) {
  result <- list(weight.X = rep(1, nx), weight.Y = rep(1, ny),
                 params.X = rep(0, 4L), params.Y = rep(0, 4L))
  if (latents == "in") {
    result$weight.X <- hyperparams[seq_len(nx)]
    result$params.X <- hyperparams[nx + seq_len(4L)]
  } else if (latents == "out") {
    result$weight.Y <- hyperparams[seq_len(ny)]
    result$params.Y <- hyperparams[ny + seq_len(4L)]
  } else if (latents == "both") {
    result$weight.X <- hyperparams[seq_len(nx)]
    result$params.X <- hyperparams[nx + seq_len(4L)]
    offset <- nx + 4L
    result$weight.Y <- hyperparams[offset + seq_len(ny)]
    result$params.Y <- hyperparams[offset + ny + seq_len(4L)]
  }
  if (rel_weight) {
    result$weight.X <- .ndrlm_relative_weights(result$weight.X)
    result$weight.Y <- .ndrlm_relative_weights(result$weight.Y)
  }
  result
}

.ndrlm_reduce <- function(data, params, weights, membership, settings) {
  ndr(data,
    min_evalue = params[1L], min_communality = params[2L],
    com_communalities = params[3L], min_R = params[4L], weight = weights,
    cor_method = settings$cor_method, cor_type = settings$cor_type,
    min_comm = settings$min_comm, Gamma = settings$Gamma,
    null_model_type = settings$null_model_type, mod_mode = settings$mod_mode,
    use_rotation = settings$use_rotation, rotation = settings$rotation,
    parallel = settings$parallel, cores = settings$cores,
    sparsecalc = settings$sparsecalc, membership = membership,
    centrality = settings$centrality, standardized = settings$standardized,
    bootstrap = FALSE, seed = settings$seed)
}

.ndrlm_design <- function(X, Y, latents, dircon, NDAin = NULL, NDAout = NULL) {
  indep <- X
  dep <- Y
  extra_x <- extra_y <- FALSE
  dropped_x <- dropped_y <- character()
  if (latents %in% c("in", "both")) {
    dropped_x <- colnames(X)[NDAin$membership == 0L]
    extra_x <- dircon && length(dropped_x) > 0L
    indep <- as.data.frame(NDAin$scores)
    colnames(indep) <- paste0("NDAin", seq_len(NDAin$factors))
    if (extra_x) indep <- cbind(indep, X[, dropped_x, drop = FALSE])
  }
  if (latents %in% c("out", "both")) {
    dropped_y <- colnames(Y)[NDAout$membership == 0L]
    extra_y <- dircon && length(dropped_y) > 0L
    dep <- as.data.frame(NDAout$scores)
    colnames(dep) <- paste0("NDAout", seq_len(NDAout$factors))
    if (extra_y) dep <- cbind(dep, Y[, dropped_y, drop = FALSE])
  }
  list(indep = as.data.frame(indep), dep = as.data.frame(dep),
       extra_vars.X = extra_x, extra_vars.Y = extra_y,
       dropped_X = dropped_x, dropped_Y = dropped_y)
}

# Accessors keep both the CRAN-era lm fits and the 0.3.0 wrapped fits usable.
.ndrlm_fit_metadata <- function(fit) {
  metadata <- attr(fit, "ndrlm_metadata", exact = TRUE)
  if (is.null(metadata)) list() else metadata
}

.ndrlm_fit_method <- function(fit) {
  if (inherits(fit, "lm")) return("lm")
  fit$method %||% NULL
}

.ndrlm_underlying_fit <- function(fit) {
  if (inherits(fit, "lm")) fit else fit$fit
}

.ndrlm_fitted_values <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  if (!is.null(metadata$fitted)) return(metadata$fitted)
  if (!inherits(fit, "lm")) return(fit$fitted)
  stats::napredict(fit$na.action, stats::fitted(fit))
}

.ndrlm_residual_values <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  if (!is.null(metadata$residuals)) return(metadata$residuals)
  if (!inherits(fit, "lm")) return(fit$residuals)
  stats::naresid(fit$na.action, stats::residuals(fit))
}

.ndrlm_fit_coefficients <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  if (!is.null(metadata$coefficients)) return(metadata$coefficients)
  if (!inherits(fit, "lm")) return(fit$coefficients %||% numeric())
  coefficients <- stats::coef(fit)
  coefficients[names(coefficients) != "(Intercept)"]
}

.ndrlm_standardize_coefficients <- function(coefficients, predictors, response) {
  coefficient_names <- names(coefficients)
  coefficients <- as.numeric(coefficients)
  names(coefficients) <- coefficient_names
  if (!length(coefficients) || is.null(coefficient_names)) return(coefficients)
  predictors <- as.data.frame(predictors, check.names = FALSE)
  response <- as.numeric(response)
  response_sd <- stats::sd(response)
  result <- stats::setNames(rep(NA_real_, length(coefficients)),
                            coefficient_names)
  if (!is.finite(response_sd) || response_sd <= .Machine$double.eps) {
    return(result)
  }
  predictor_sd <- vapply(predictors, stats::sd, numeric(1L))
  matched <- intersect(coefficient_names, names(predictor_sd))
  result[matched] <- coefficients[matched] * predictor_sd[matched] / response_sd
  result
}

.ndrlm_fit_standardized_coefficients <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  if (!is.null(metadata$standardized_coefficients)) {
    return(metadata$standardized_coefficients)
  }
  if (!inherits(fit, "lm") && !is.null(fit$standardized_coefficients)) {
    return(fit$standardized_coefficients)
  }
  coefficients <- .ndrlm_fit_coefficients(fit)
  underlying <- .ndrlm_underlying_fit(fit)
  if (!inherits(underlying, "lm")) return(coefficients)
  model_frame <- tryCatch(stats::model.frame(underlying),
                          error = function(e) NULL)
  design <- tryCatch(stats::model.matrix(underlying),
                     error = function(e) NULL)
  if (is.null(model_frame) || is.null(design)) return(coefficients)
  response <- tryCatch(stats::model.response(model_frame),
                       error = function(e) NULL)
  if (is.null(response)) return(coefficients)
  intercept <- colnames(design) == "(Intercept)"
  design <- design[, !intercept, drop = FALSE]
  if (!all(names(coefficients) %in% colnames(design)) &&
      ncol(design) == length(coefficients)) {
    colnames(design) <- names(coefficients)
  }
  standardized <- .ndrlm_standardize_coefficients(
    coefficients, design, response)
  if (all(!is.finite(standardized))) coefficients else standardized
}

.ndrlm_fit_p_values <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  if (!is.null(metadata$p_values)) return(metadata$p_values)
  if (!inherits(fit, "lm")) return(fit$p_values %||% numeric())
  coefficient_table <- tryCatch(summary(fit)$coefficients,
                                error = function(e) NULL)
  coefficients <- .ndrlm_fit_coefficients(fit)
  p_values <- stats::setNames(rep(NA_real_, length(coefficients)),
                              names(coefficients))
  if (!is.null(coefficient_table) && ncol(coefficient_table) >= 4L) {
    matched <- intersect(names(coefficients), rownames(coefficient_table))
    p_values[matched] <- coefficient_table[matched, 4L]
  }
  p_values
}

.ndrlm_fit_response <- function(fit) {
  metadata <- .ndrlm_fit_metadata(fit)
  response <- metadata$response %||% fit$response %||% NULL
  if (length(response) && !is.na(response[1L]) &&
      nzchar(as.character(response[1L]))) {
    return(as.character(response[1L]))
  }
  underlying <- .ndrlm_underlying_fit(fit)
  if (inherits(underlying, "lm")) {
    formula <- tryCatch(stats::formula(underlying),
                        error = function(e) NULL)
    if (!is.null(formula) && length(formula) >= 3L) {
      response <- all.vars(formula[[2L]])
      if (length(response)) return(response[1L])
    }
  }
  NULL
}

.ndrlm_lm_formula <- function(response, predictors, env = parent.frame()) {
  terms <- lapply(predictors, as.name)
  rhs <- if (length(terms) == 1L) terms[[1L]] else
    Reduce(function(left, right) call("+", left, right), terms)
  formula <- call("~", as.name(response), rhs)
  class(formula) <- "formula"
  environment(formula) <- env
  formula
}

.ndrlm_fit_one <- function(y, x, response, regression, alpha, lambda, loess_span) {
  original_names <- colnames(x)
  safe_names <- make.unique(c(".response", make.names(original_names)))[-1L]
  safe_x <- x
  colnames(safe_x) <- safe_names
  data <- data.frame(.response = as.numeric(y), safe_x, check.names = FALSE)
  keep <- stats::complete.cases(data)
  train <- data[keep, , drop = FALSE]
  if (nrow(train) < 3L) stop("Too few complete observations for regression.", call. = FALSE)
  formula <- stats::reformulate(safe_names, response = ".response")
  if (regression == "lm") {
    lm_data <- data.frame(.response = as.numeric(y), x, check.names = FALSE)
    lm_train <- lm_data[keep, , drop = FALSE]
    lm_formula <- .ndrlm_lm_formula(".response", original_names)
    fit <- stats::lm(lm_formula, data = lm_train)
    fitted <- rep(NA_real_, nrow(data)); fitted[keep] <- stats::fitted(fit)
    coefficients <- stats::coef(fit)[-1L]
    names(coefficients) <- original_names
    standardized_coefficients <- .ndrlm_standardize_coefficients(
      coefficients, lm_train[, original_names, drop = FALSE],
      lm_train$.response)
    p_values <- rep(NA_real_, length(original_names))
    names(p_values) <- original_names
    coefficient_table <- summary(fit)$coefficients
    rows <- intersect(original_names, rownames(coefficient_table))
    p_values[match(rows, original_names)] <- coefficient_table[rows, 4L]
    attr(fit, "ndrlm_metadata") <- list(
      method = regression, fitted = fitted,
      residuals = data$.response - fitted, response = response,
      x_names = original_names, safe_names = original_names,
      coefficients = coefficients,
      standardized_coefficients = standardized_coefficients,
      p_values = p_values
    )
    return(fit)
  }
  if (regression == "loess") {
    fit <- stats::loess(formula, data = train, span = loess_span,
                        control = stats::loess.control(surface = "direct"))
    fitted <- rep(NA_real_, nrow(data)); fitted[keep] <- stats::fitted(fit)
    weights <- vapply(train[, -1L, drop = FALSE], function(z)
      stats::cor(z, stats::fitted(fit), use = "complete.obs"), numeric(1))
    names(weights) <- original_names
    loess_p <- rep(NA_real_, length(weights)); names(loess_p) <- original_names
    return(list(method = regression, fit = fit, fitted = fitted,
                residuals = data$.response - fitted, response = response,
                x_names = original_names, safe_names = safe_names, coefficients = weights,
                standardized_coefficients = weights,
                p_values = loess_p))
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for ridge and elastic-net regression.", call. = FALSE)
  }
  matrix_x <- stats::model.matrix(~ . - 1, data = train[, -1L, drop = FALSE])
  model_alpha <- if (regression == "ridge") 0 else alpha
  if (is.null(lambda)) {
    folds <- min(10L, max(3L, floor(nrow(train) / 5L)))
    fit <- glmnet::cv.glmnet(matrix_x, train$.response, alpha = model_alpha,
                            nfolds = folds)
    selected_lambda <- fit$lambda.min
    pred <- as.numeric(stats::predict(fit, newx = matrix_x, s = "lambda.min"))
    coefficients <- as.matrix(stats::coef(fit, s = "lambda.min"))[-1L, 1L]
  } else {
    fit <- glmnet::glmnet(matrix_x, train$.response, alpha = model_alpha,
                          lambda = lambda)
    selected_lambda <- lambda[1L]
    pred <- as.numeric(stats::predict(fit, newx = matrix_x, s = selected_lambda))
    coefficients <- as.matrix(stats::coef(fit, s = selected_lambda))[-1L, 1L]
  }
  fitted <- rep(NA_real_, nrow(data)); fitted[keep] <- pred
  names(coefficients) <- original_names
  standardization_data <- train[, -1L, drop = FALSE]
  colnames(standardization_data) <- original_names
  standardized_coefficients <- .ndrlm_standardize_coefficients(
    coefficients, standardization_data, train$.response)
  list(method = regression, fit = fit, fitted = fitted,
       residuals = data$.response - fitted, response = response,
       x_names = original_names, safe_names = safe_names,
       design_names = colnames(matrix_x),
       coefficients = coefficients,
       standardized_coefficients = standardized_coefficients,
       p_values = stats::setNames(rep(NA_real_, length(coefficients)), original_names),
       lambda = selected_lambda, alpha = model_alpha)
}

.ndrlm_predict_one <- function(model, newdata, se.fit = FALSE,
                               interval = "none", level = 0.95, ...) {
  if (inherits(model, "lm")) {
    metadata <- .ndrlm_fit_metadata(model)
    if (length(metadata)) {
      x <- as.data.frame(newdata[, metadata$x_names, drop = FALSE])
      colnames(x) <- metadata$safe_names
    } else {
      x <- as.data.frame(newdata)
    }
    return(stats::predict(model, newdata = x, se.fit = se.fit,
                          interval = interval, level = level, ...))
  }
  x <- as.data.frame(newdata[, model$x_names, drop = FALSE])
  colnames(x) <- model$safe_names
  if (model$method == "lm") {
    return(stats::predict(model$fit, newdata = x, se.fit = se.fit,
                          interval = interval, level = level, ...))
  }
  if (se.fit || interval != "none") {
    warning("Standard errors and prediction intervals are only available for lm fits.",
            call. = FALSE)
  }
  if (model$method == "loess") return(stats::predict(model$fit, newdata = x))
  matrix_x <- stats::model.matrix(~ . - 1, data = x)
  missing_columns <- setdiff(model$design_names, colnames(matrix_x))
  if (length(missing_columns)) {
    matrix_x <- cbind(matrix_x, matrix(0, nrow(matrix_x), length(missing_columns),
      dimnames = list(NULL, missing_columns)))
  }
  matrix_x <- matrix_x[, model$design_names, drop = FALSE]
  as.numeric(stats::predict(model$fit, newx = matrix_x,
    s = if (inherits(model$fit, "cv.glmnet")) "lambda.min" else model$lambda))
}

.ndrlm_fit_set <- function(dep, indep, regression, alpha, lambda,
                           loess_span, target, fit_weights, keep_fits = TRUE) {
  fits <- vector("list", ncol(dep))
  metric <- numeric(ncol(dep))
  for (i in seq_len(ncol(dep))) {
    fits[[i]] <- .ndrlm_fit_one(dep[[i]], indep, colnames(dep)[i], regression,
                                alpha, lambda, loess_span)
    metric[i] <- .nda_metric(dep[[i]], .ndrlm_fitted_values(fits[[i]]),
                             fits[[i]], target)
  }
  names(fits) <- colnames(dep)
  names(metric) <- colnames(dep)
  aggregate <- if (is.null(fit_weights)) mean(metric) else
    stats::weighted.mean(metric, fit_weights)
  list(fits = if (keep_fits) fits else NULL, metric = metric, aggregate = aggregate)
}

#' Network-Based Dimensionality Reduction and Regression
#'
#' Fits supervised models to original or NDA-reduced input and output
#' variables, with optional NSGA-II tuning.
#'
#' @param Y Numeric dependent-variable data.
#' @param X Numeric independent-variable data with the same rows as \code{Y}.
#' @param latents Use latent variables for inputs (\code{"in"}), outputs
#'   (\code{"out"}), both, or neither.
#' @param dircon Retain NDA-excluded variables as direct predictors or responses.
#' @param optimize Tune indicator weights and four NDA thresholds with NSGA-II.
#' @param target One of adjusted R-squared, R-squared, MAE, MAPE, MASE, MSE,
#'   RMSE, AIC, or BIC.
#' @param rel_weight Constrain and normalize indicator weights to be
#'   nonnegative relative weights.
#' @param cor_method,cor_type,min_comm,Gamma,null_model_type,mod_mode,use_rotation,rotation
#'   Arguments forwarded to \code{ndr}.
#' @param pareto Keep one objective per original dependent variable during
#'   NSGA-II optimization.
#' @param fit_weights Optional weights for aggregating response-specific
#'   objectives when their number matches the fitted responses.
#' @param lower.bounds.x,upper.bounds.x Lower and upper NSGA-II bounds for
#'   input-indicator weights.
#' @param lower.bounds.latentx,upper.bounds.latentx Bounds for input-side
#'   \code{min_evalue}, \code{min_communality}, \code{com_communalities}, and
#'   \code{min_R}.
#' @param lower.bounds.y,upper.bounds.y Lower and upper NSGA-II bounds for
#'   output-indicator weights.
#' @param lower.bounds.latenty,upper.bounds.latenty Bounds for the four
#'   output-side NDA thresholds.
#' @param popsize,generations,cprob,cdist,mprob,mdist NSGA-II population,
#'   generation, crossover, and mutation settings passed to \code{mco::nsga2}.
#' @param seed Optional seed; repeated optimization runs use consecutive seeds.
#' @param regression Regression engine: ordinary least squares (\code{"lm"}),
#'   ridge (\code{"ridge"}), LOESS (\code{"loess"}), or elastic net
#'   (\code{"elastic_net"}).
#' @param model Backward-friendly alias that overrides \code{regression}.
#' @param alpha Elastic-net mixing value in [0, 1]. Ridge always uses zero.
#' @param lambda Optional penalization value. If omitted, cross-validation
#'   selects lambda.
#' @param loess_span LOESS smoothing span.
#' @param restarts Number of independent NSGA-II runs; the best feasible
#'   candidate is refitted.
#' @param parallel,cores,sparsecalc,centrality,standardized Performance and
#'   latent-score settings forwarded to \code{ndr}.
#' @param membership.X,membership.Y Optional fixed NDA community memberships
#'   for the input and output sides.
#'
#' @details
#' NDRLM is predictive and exploratory; its paths do not by themselves
#' establish causal effects. Ridge and elastic net use \pkg{glmnet}. LOESS is
#' intended for small predictor sets and does not provide linear-path
#' significance tests. AIC and BIC for non-likelihood engines are
#' residual-based approximations.
#'
#' The class and structure of each element of \code{fits} depend on the selected
#' regression engine. With \code{regression = "lm"}, each \code{fits[[i]]} is a
#' genuine \code{lm} object and can be passed directly to \code{summary},
#' \code{coef}, \code{AIC}, \code{BIC}, \code{predict.lm}, and
#' \code{lm.beta::lm.beta}. The \pkg{lm.beta} package is meaningful only for
#' ordinary least-squares fits; it must not be applied to ridge, elastic-net,
#' or LOESS fits.
#'
#' With \code{regression = "ridge"}, \code{"elastic_net"}, or \code{"loess"},
#' each \code{fits[[i]]} is an NDRLM engine wrapper. It contains the engine name
#' in \code{method}, the underlying fitted object in \code{fit}, and common
#' fitted-value, residual, coefficient, and predictor metadata. Penalized fits
#' additionally store the selected \code{lambda}, effective \code{alpha}, and
#' standardized coefficients calculated on the complete training observations.
#' The plot method applies the same standardization to ordinary least-squares
#' fits, including compatible legacy \code{lm} objects. For engine-independent
#' code, use \code{summary}, \code{predict},
#' \code{fitted}, \code{residuals}, and \code{plot} on the enclosing
#' \code{ndrlm} object instead of relying on these internal fields.
#'
#' The NDRLM methods retain backward compatibility with both legacy objects
#' that stored ordinary least-squares fits directly and early 0.3.0 objects
#' that wrapped ordinary least-squares fits. Newly fitted ordinary
#' least-squares models use direct \code{lm} objects, while the newer
#' non-\code{lm} engines retain their wrappers.
#'
#' @return An \code{ndrlm} object containing response-specific metrics, fitted
#'   NDA objects, optimized hyperparameters, training designs, NSGA-II results,
#'   an \code{NSGA_trace} evaluation history when optimization is enabled, and
#'   a response-named \code{fits} list. Elements of \code{fits} are
#'   \code{lm} objects for ordinary least squares and NDRLM engine wrappers for
#'   ridge, elastic net, and LOESS.
#'
#' @references Kosztyan, Z. T., Kurbucz, M. T., Katona, A. I., and Khan, Z.
#'   Network-based dimensionality reduction and structural modeling.
#'
#' @examples
#' X <- as.data.frame(freeny.x)
#' Y <- data.frame(income = freeny.y)
#'
#' fit_lm <- ndrlm(Y, X, latents = "none", regression = "lm")
#' inherits(fit_lm$fits[[1]], "lm")
#' summary(fit_lm$fits[[1]])
#' predict(fit_lm, X[1:3, ])
#'
#' if (requireNamespace("lm.beta", quietly = TRUE)) {
#'   summary(lm.beta::lm.beta(fit_lm$fits[[1]]))
#' }
#'
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   fit_ridge <- ndrlm(Y, X, latents = "none", regression = "ridge",
#'                      lambda = 0.1)
#'   fit_ridge$fits[[1]]$method
#'   predict(fit_ridge, X[1:3, ])
#' }
#'
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   fit_en <- ndrlm(Y, X, latents = "none", regression = "elastic_net",
#'                   alpha = 0.5, lambda = 0.1)
#'   predict(fit_en, X[1:3, ])
#' }
#'
#' @keywords multivariate
#' @export
ndrlm <- function(Y, X, latents = "in", dircon = FALSE, optimize = TRUE,
                  target = "adj.r.square", rel_weight = FALSE,
                  cor_method = 1, cor_type = 1, min_comm = 2, Gamma = 1,
                  null_model_type = 4, mod_mode = 1, use_rotation = FALSE,
                  rotation = "oblimin", pareto = FALSE, fit_weights = NULL,
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
                  model = NULL, alpha = 0.5, lambda = NULL,
                  loess_span = 0.75, restarts = 1L, parallel = FALSE,
                  cores = 1L, sparsecalc = FALSE,
                  centrality = "eigenvector", standardized = TRUE,
                  membership.X = NULL, membership.Y = NULL) {
  call <- match.call()
  X <- as.data.frame(as.matrix(.nda_numeric_matrix(X, "X")))
  Y <- as.data.frame(as.matrix(.nda_numeric_matrix(Y, "Y")))
  if (nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows.", call. = FALSE)
  latents <- match.arg(tolower(latents), c("in", "out", "both", "none"))
  if (!is.null(model)) regression <- model
  regression <- match.arg(tolower(gsub("[- ]", "_", regression)),
                          c("lm", "ridge", "loess", "elastic_net"))
  targets <- c("adj.r.square", "r.square", "MAE", "MAPE", "MASE", "MSE",
               "RMSE", "AIC", "BIC")
  if (!target %in% targets) stop("Unknown target metric: ", target, call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0 || alpha > 1) {
    stop("alpha must be between zero and one.", call. = FALSE)
  }
  if (!is.numeric(restarts) || length(restarts) != 1L || restarts < 1L) {
    stop("restarts must be a positive integer.", call. = FALSE)
  }
  if (latents == "none") optimize <- pareto <- FALSE
  if (optimize && (popsize < 4L || popsize %% 4L != 0L)) {
    stop("popsize must be a positive multiple of four for NSGA-II.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(fit_weights) && length(fit_weights) != ncol(Y)) {
    stop("fit_weights must have one value per original dependent variable.", call. = FALSE)
  }
  nx <- ncol(X); ny <- ncol(Y)
  initial <- switch(latents,
    `in` = c(rep(1, nx), rep(0, 4L)),
    `out` = c(rep(1, ny), rep(0, 4L)),
    both = c(rep(1, nx), rep(0, 4L), rep(1, ny), rep(0, 4L)),
    none = numeric())
  lower <- switch(latents,
    `in` = c(lower.bounds.x, lower.bounds.latentx),
    `out` = c(lower.bounds.y, lower.bounds.latenty),
    both = c(lower.bounds.x, lower.bounds.latentx,
             lower.bounds.y, lower.bounds.latenty), none = numeric())
  upper <- switch(latents,
    `in` = c(upper.bounds.x, upper.bounds.latentx),
    `out` = c(upper.bounds.y, upper.bounds.latenty),
    both = c(upper.bounds.x, upper.bounds.latentx,
             upper.bounds.y, upper.bounds.latenty), none = numeric())
  if (rel_weight) {
    lower[seq_along(lower) <= length(initial)] <- pmax(lower[seq_along(lower) <= length(initial)], 0)
  }
  if (length(lower) != length(initial) || length(upper) != length(initial)) {
    stop("Optimization bounds do not match the number of hyperparameters.", call. = FALSE)
  }
  settings <- list(cor_method = cor_method, cor_type = cor_type,
    min_comm = min_comm, Gamma = Gamma, null_model_type = null_model_type,
    mod_mode = mod_mode, use_rotation = use_rotation, rotation = rotation,
    parallel = parallel, cores = cores, sparsecalc = sparsecalc,
    centrality = centrality, standardized = standardized, seed = seed)

  build <- function(hyperparams, keep_fits = TRUE) {
    parsed <- .ndrlm_unpack(hyperparams, latents, nx, ny, rel_weight)
    NDAin <- NDAout <- NULL
    if (latents %in% c("in", "both")) {
      NDAin <- .ndrlm_reduce(X, parsed$params.X, parsed$weight.X,
                             membership.X, settings)
    }
    if (latents %in% c("out", "both")) {
      NDAout <- .ndrlm_reduce(Y, parsed$params.Y, parsed$weight.Y,
                              membership.Y, settings)
    }
    design <- .ndrlm_design(X, Y, latents, dircon, NDAin, NDAout)
    selected_fit_weights <- fit_weights
    if (!is.null(selected_fit_weights) && length(selected_fit_weights) != ncol(design$dep)) {
      selected_fit_weights <- NULL
    }
    fitted_set <- .ndrlm_fit_set(design$dep, design$indep, regression, alpha,
      lambda, loess_span, target, selected_fit_weights, keep_fits)
    c(list(parsed = parsed, NDAin = NDAin, NDAout = NDAout, design = design),
      fitted_set)
  }

  maximize <- target %in% c("adj.r.square", "r.square")
  trace_state <- new.env(parent = emptyenv())
  trace_state$recording <- FALSE
  trace_state$restart <- 0L
  trace_state$evaluation <- 0L
  trace_state$rows <- list()
  objective <- function(hyperparams) {
    value <- tryCatch(build(hyperparams, keep_fits = FALSE)$metric,
      error = function(e) rep(if (maximize) -1e12 else 1e12,
                              if (pareto) max(1L, ncol(Y)) else 1L))
    if (any(!is.finite(value))) value <- rep(if (maximize) -1e12 else 1e12,
                                             length(value))
    if (pareto && length(value) != ncol(Y)) value <- rep(mean(value), ncol(Y))
    if (!pareto) value <- mean(value)
    cost <- if (maximize) -value else value
    if (isTRUE(trace_state$recording)) {
      trace_state$evaluation <- trace_state$evaluation + 1L
      scalar_cost <- mean(cost)
      target_value <- if (maximize) -scalar_cost else scalar_cost
      trace_state$rows[[length(trace_state$rows) + 1L]] <- data.frame(
        restart = trace_state$restart,
        evaluation = trace_state$evaluation,
        generation = max(0L, ceiling(trace_state$evaluation / popsize) - 1L),
        target_value = target_value,
        cost = scalar_cost
      )
    }
    cost
  }

  optimized <- FALSE
  NSGA <- NULL
  NSGA_trace <- NULL
  hyperparams <- initial
  if (optimize && length(initial)) {
    if (!requireNamespace("mco", quietly = TRUE)) {
      stop("Package 'mco' is required when optimize = TRUE.", call. = FALSE)
    }
    runs <- vector("list", as.integer(restarts))
    candidates <- list()
    candidate_costs <- numeric()
    for (run in seq_len(as.integer(restarts))) {
      if (!is.null(seed)) set.seed(seed + run - 1L)
      trace_state$recording <- TRUE
      trace_state$restart <- run
      trace_state$evaluation <- 0L
      runs[[run]] <- mco::nsga2(objective, idim = length(initial),
        odim = if (pareto) ncol(Y) else 1L, lower.bounds = lower,
        upper.bounds = upper, popsize = popsize, generations = generations,
        cprob = cprob, cdist = cdist, mprob = mprob, mdist = mdist,
        vectorized = FALSE)
      trace_state$recording <- FALSE
      for (row in seq_len(nrow(runs[[run]]$par))) {
        candidate <- runs[[run]]$par[row, ]
        value <- objective(candidate)
        candidates[[length(candidates) + 1L]] <- candidate
        candidate_costs <- c(candidate_costs, mean(value))
      }
    }
    best <- which.min(candidate_costs)
    if (length(best) && is.finite(candidate_costs[best])) {
      hyperparams <- candidates[[best]]
      optimized <- TRUE
    } else warning("Optimization did not find a feasible solution; initial values were used.",
                   call. = FALSE)
    NSGA <- if (length(runs) == 1L) runs[[1L]] else runs
    if (length(trace_state$rows)) {
      NSGA_trace <- do.call(rbind, trace_state$rows)
      split_index <- split(seq_len(nrow(NSGA_trace)), NSGA_trace$restart)
      NSGA_trace$best_target <- NA_real_
      for (index in split_index) {
        NSGA_trace$best_target[index] <- if (maximize) {
          cummax(NSGA_trace$target_value[index])
        } else cummin(NSGA_trace$target_value[index])
      }
      rownames(NSGA_trace) <- NULL
    }
  }
  final <- tryCatch(build(hyperparams, keep_fits = TRUE), error = function(e) {
    if (!identical(hyperparams, initial)) {
      warning("The optimized solution failed during refitting; initial values were used.",
              call. = FALSE)
      optimized <<- FALSE
      hyperparams <<- initial
      build(initial, keep_fits = TRUE)
    } else stop(e)
  })

  parsed <- final$parsed
  result <- list(Call = call, target = target,
    fval = if (pareto) final$metric else final$aggregate,
    fit_stats = data.frame(response = names(final$metric), metric = final$metric,
                           row.names = NULL),
    hyperparams = hyperparams, pareto = pareto, X = X, Y = Y,
    latents = latents, fits = final$fits, regression = regression,
    alpha = alpha, lambda = lambda, loess_span = loess_span,
    optimized = optimized, extra_vars.X = final$design$extra_vars.X,
    extra_vars.Y = final$design$extra_vars.Y,
    dircon_X = final$design$dropped_X, dircon_Y = final$design$dropped_Y,
    indep = final$design$indep, dep = final$design$dep,
    indep_names = colnames(final$design$indep),
    dep_names = colnames(final$design$dep), seed = seed, fn = "NDRLM")
  if (latents %in% c("in", "both")) {
    result$NDAin <- final$NDAin
    result$NDAin_weight <- parsed$weight.X
    result$NDAin_min_evalue <- parsed$params.X[1L]
    result$NDAin_min_communality <- parsed$params.X[2L]
    result$NDAin_com_communalities <- parsed$params.X[3L]
    result$NDAin_min_R <- parsed$params.X[4L]
  }
  if (latents %in% c("out", "both")) {
    result$NDAout <- final$NDAout
    result$NDAout_weight <- parsed$weight.Y
    result$NDAout_min_evalue <- parsed$params.Y[1L]
    result$NDAout_min_communality <- parsed$params.Y[2L]
    result$NDAout_com_communalities <- parsed$params.Y[3L]
    result$NDAout_min_R <- parsed$params.Y[4L]
  }
  if (optimize) {
    result$NSGA <- NSGA
    result$NSGA_trace <- NSGA_trace
  }
  class(result) <- c("ndrlm", "list")
  result
}


