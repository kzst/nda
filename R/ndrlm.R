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

.ndrlm_fit_one <- function(y, x, response, regression, alpha, lambda, loess_span) {
  original_names <- colnames(x)
  safe_names <- paste0(".x", seq_len(ncol(x)))
  safe_x <- x
  colnames(safe_x) <- safe_names
  data <- data.frame(.response = as.numeric(y), safe_x, check.names = FALSE)
  keep <- stats::complete.cases(data)
  train <- data[keep, , drop = FALSE]
  if (nrow(train) < 3L) stop("Too few complete observations for regression.", call. = FALSE)
  formula <- stats::reformulate(safe_names, response = ".response")
  if (regression == "lm") {
    fit <- stats::lm(formula, data = train)
    fitted <- rep(NA_real_, nrow(data)); fitted[keep] <- stats::fitted(fit)
    coefficients <- stats::coef(fit)[-1L]
    names(coefficients) <- original_names
    p_values <- rep(NA_real_, length(original_names))
    names(p_values) <- original_names
    coefficient_table <- summary(fit)$coefficients
    rows <- intersect(safe_names, rownames(coefficient_table))
    p_values[match(rows, safe_names)] <- coefficient_table[rows, 4L]
    return(list(method = regression, fit = fit, fitted = fitted,
                residuals = data$.response - fitted, response = response,
                x_names = original_names, safe_names = safe_names,
                coefficients = coefficients, p_values = p_values))
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
  list(method = regression, fit = fit, fitted = fitted,
       residuals = data$.response - fitted, response = response,
       x_names = original_names, safe_names = safe_names,
       design_names = colnames(matrix_x),
       coefficients = coefficients,
       p_values = stats::setNames(rep(NA_real_, length(coefficients)), original_names),
       lambda = selected_lambda, alpha = model_alpha)
}

.ndrlm_predict_one <- function(model, newdata, se.fit = FALSE,
                               interval = "none", level = 0.95, ...) {
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
    metric[i] <- .nda_metric(dep[[i]], fits[[i]]$fitted, fits[[i]], target)
  }
  names(fits) <- colnames(dep)
  names(metric) <- colnames(dep)
  aggregate <- if (is.null(fit_weights)) mean(metric) else
    stats::weighted.mean(metric, fit_weights)
  list(fits = if (keep_fits) fits else NULL, metric = metric, aggregate = aggregate)
}

#' Network-based dimensionality reduction and regression
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
  objective <- function(hyperparams) {
    value <- tryCatch(build(hyperparams, keep_fits = FALSE)$metric,
      error = function(e) rep(if (maximize) -1e12 else 1e12,
                              if (pareto) max(1L, ncol(Y)) else 1L))
    if (any(!is.finite(value))) value <- rep(if (maximize) -1e12 else 1e12,
                                             length(value))
    if (pareto && length(value) != ncol(Y)) value <- rep(mean(value), ncol(Y))
    if (!pareto) value <- mean(value)
    if (maximize) -value else value
  }

  optimized <- FALSE
  NSGA <- NULL
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
      runs[[run]] <- mco::nsga2(objective, idim = length(initial),
        odim = if (pareto) ncol(Y) else 1L, lower.bounds = lower,
        upper.bounds = upper, popsize = popsize, generations = generations,
        cprob = cprob, cdist = cdist, mprob = mprob, mdist = mdist,
        vectorized = FALSE)
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
  if (optimize) result$NSGA <- NSGA
  class(result) <- c("ndrlm", "list")
  result
}
