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



# Internal utilities shared by NDA and NDRLM.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.nda_numeric_matrix <- function(x, name = "x", min_cols = 1L) {
  if (inherits(x, "Matrix")) {
    if (!is.numeric(x)) stop(name, " must be numeric.", call. = FALSE)
    out <- x
  } else {
    out <- as.matrix(x)
    if (!is.numeric(out)) {
      stop(name, " must be a numeric matrix or data frame.", call. = FALSE)
    }
  }
  if (length(dim(out)) != 2L || nrow(out) < 1L || ncol(out) < min_cols) {
    stop(name, " has invalid dimensions.", call. = FALSE)
  }
  if (is.null(colnames(out))) colnames(out) <- paste0("V", seq_len(ncol(out)))
  out
}

.nda_choice <- function(x, choices, numeric_choices = NULL,
                        name = deparse(substitute(x))) {
  if (length(x) != 1L || is.na(x)) {
    stop(name, " must have length one and cannot be missing.", call. = FALSE)
  }
  if (is.numeric(x) && !is.null(numeric_choices)) {
    key <- as.character(as.integer(x))
    if (!key %in% names(numeric_choices)) stop("Unknown ", name, ": ", x, call. = FALSE)
    return(unname(numeric_choices[[key]]))
  }
  x <- tolower(gsub("[- ]", "_", as.character(x)))
  aliases <- c(bivariate = "full", semi_partial = "semipartial",
               semi = "semipartial", distance_correlation = "distance",
               gaussian = "gaussian_rank", page_rank = "pagerank")
  if (x %in% names(aliases)) x <- aliases[[x]]
  match.arg(x, choices)
}

.nda_scale <- function(x, center = NULL, scale = NULL, enabled = TRUE) {
  x <- as.matrix(x)
  if (!enabled) {
    return(list(x = x, center = rep(0, ncol(x)), scale = rep(1, ncol(x))))
  }
  center <- center %||% colMeans(x, na.rm = TRUE)
  scale <- scale %||% apply(x, 2L, stats::sd, na.rm = TRUE)
  center[!is.finite(center)] <- 0
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  z[!is.finite(z)] <- 0
  list(x = z, center = center, scale = scale)
}

.nda_apply_pairs <- function(p, fun, parallel = FALSE, cores = 1L) {
  pairs <- which(upper.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
  work <- lapply(seq_len(nrow(pairs)), function(k) pairs[k, ])
  if (!length(work)) return(list())
  if (parallel && cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(work, fun, mc.cores = min(as.integer(cores), length(work)))
  } else {
    if (parallel && cores > 1L && .Platform$OS.type == "windows") {
      warning("Pairwise parallelism uses a single core on Windows.", call. = FALSE)
    }
    lapply(work, fun)
  }
}

.nda_pairwise_matrix <- function(x, statistic, test_fun = NULL, test = FALSE,
                                 adjust = "BH", alpha = NULL,
                                 parallel = FALSE, cores = 1L) {
  x <- as.matrix(.nda_numeric_matrix(x))
  p <- ncol(x)
  r <- diag(1, p)
  dimnames(r) <- list(colnames(x), colnames(x))
  pv <- matrix(0, p, p, dimnames = dimnames(r))
  vals <- .nda_apply_pairs(p, function(ij) {
    a <- x[, ij[1L]]
    b <- x[, ij[2L]]
    keep <- stats::complete.cases(a, b)
    if (sum(keep) < 3L) return(c(r = NA_real_, p = NA_real_))
    estimate <- statistic(a[keep], b[keep])
    prob <- if (test) test_fun(a[keep], b[keep]) else NA_real_
    c(r = estimate, p = prob)
  }, parallel, cores)
  pairs <- which(upper.tri(r), arr.ind = TRUE)
  for (k in seq_along(vals)) {
    i <- pairs[k, 1L]
    j <- pairs[k, 2L]
    r[i, j] <- r[j, i] <- vals[[k]]["r"]
    if (test) pv[i, j] <- pv[j, i] <- vals[[k]]["p"]
  }
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  if (!test) return(r)
  upper <- upper.tri(pv)
  pv[upper] <- stats::p.adjust(pv[upper], method = adjust)
  pv[lower.tri(pv)] <- t(pv)[lower.tri(pv)]
  diag(pv) <- 0
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
      stop("alpha must be a number between zero and one.", call. = FALSE)
    }
    r[pv > alpha] <- 0
    diag(r) <- 1
  }
  list(r = r, p = pv)
}

.nda_inverse <- function(x) {
  out <- tryCatch(solve(x), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out)) || any(abs(diag(out)) <= .Machine$double.eps)) {
    warning("A generalized inverse was used for a singular association matrix.",
            call. = FALSE)
    out <- MASS::ginv(x)
  }
  out
}

.nda_pvalues <- function(r, n, controls = 0L, adjust = "BH", alpha = NULL) {
  df <- max(1, n - controls - 2L)
  z <- abs(r) * sqrt(df / pmax(1 - r^2, .Machine$double.eps))
  p <- 2 * stats::pt(-z, df = df)
  diag(p) <- 0
  upper <- upper.tri(p)
  p[upper] <- stats::p.adjust(p[upper], method = adjust)
  p[lower.tri(p)] <- t(p)[lower.tri(p)]
  if (!is.null(alpha)) {
    r[p > alpha] <- 0
    diag(r) <- 1
  }
  list(r = r, p = p)
}

.nda_cor_method <- function(x) {
  .nda_choice(x,
    c("pearson", "spearman", "kendall", "distance", "gaussian_rank", "biweight"),
    c(`1` = "pearson", `2` = "spearman", `3` = "kendall", `4` = "distance",
      `5` = "gaussian_rank", `6` = "biweight"), "cor_method")
}

.nda_cor_type <- function(x) {
  .nda_choice(x, c("full", "partial", "semipartial"),
              c(`1` = "full", `2` = "partial", `3` = "semipartial"),
              "cor_type")
}

.nda_biweight <- function(x) {
  x <- as.matrix(x)
  med <- apply(x, 2L, stats::median, na.rm = TRUE)
  mad <- apply(x, 2L, stats::mad, na.rm = TRUE, constant = 1)
  mad[!is.finite(mad) | mad == 0] <- 1
  u <- sweep(sweep(x, 2L, med, "-"), 2L, 9 * mad, "/")
  w <- (1 - pmin(u^2, 1))^2
  z <- sweep(x, 2L, med, "-") * w
  stats::cor(z, use = "pairwise.complete.obs")
}

.nda_association <- function(x, cor_method = "pearson", cor_type = "full",
                             test = FALSE, adjust = "BH", alpha = NULL,
                             parallel = FALSE, cores = 1L) {
  method <- .nda_cor_method(cor_method)
  type <- .nda_cor_type(cor_type)
  x <- as.matrix(.nda_numeric_matrix(x))
  if (ncol(x) == 1L) {
    one <- matrix(1, 1, 1, dimnames = list(colnames(x), colnames(x)))
    return(if (test) list(r = one, p = matrix(0, 1, 1, dimnames = dimnames(one))) else one)
  }
  if (method == "distance") {
    return(switch(type,
      full = dCor(x, test = test, adjust = adjust, alpha = alpha,
                  parallel = parallel, cores = cores),
      partial = pdCor(x, test = test, adjust = adjust, alpha = alpha,
                      parallel = parallel, cores = cores),
      semipartial = spdCor(x, test = test, adjust = adjust, alpha = alpha,
                           parallel = parallel, cores = cores)))
  }
  transformed <- x
  base_method <- method
  if (method == "gaussian_rank") {
    transformed <- apply(x, 2L, function(z) {
      rnk <- rank(z, na.last = "keep", ties.method = "average")
      stats::qnorm((rnk - 0.5) / sum(!is.na(rnk)))
    })
    base_method <- "pearson"
  }
  full <- if (method == "biweight") .nda_biweight(transformed) else
    stats::cor(transformed, use = "pairwise.complete.obs", method = base_method)
  full[!is.finite(full)] <- 0
  diag(full) <- 1
  r <- switch(type,
    full = full,
    partial = {
      inv <- .nda_inverse(full)
      den <- sqrt(pmax(abs(diag(inv)), .Machine$double.eps))
      out <- -inv / outer(den, den)
      diag(out) <- 1
      out
    },
    semipartial = {
      if (method %in% c("pearson", "spearman", "kendall") &&
          requireNamespace("ppcor", quietly = TRUE)) {
        ppcor::spcor(transformed, method = base_method)$estimate
      } else {
        inv <- .nda_inverse(full)
        residual_precision <- diag(inv) - t(t(inv^2) / diag(inv))
        den <- sqrt(pmax(diag(full), .Machine$double.eps)) *
          sqrt(pmax(abs(residual_precision), .Machine$double.eps))
        out <- -stats::cov2cor(inv) / den
        diag(out) <- 1
        out
      }
    })
  dimnames(r) <- list(colnames(x), colnames(x))
  if (!test) return(r)
  .nda_pvalues(r, nrow(x), controls = if (type == "full") 0L else ncol(x) - 2L,
               adjust = adjust, alpha = alpha)
}

.nda_metric <- function(observed, predicted, fit = NULL, target = "adj.r.square") {
  keep <- stats::complete.cases(observed, predicted)
  y <- observed[keep]
  p <- predicted[keep]
  rss <- sum((y - p)^2)
  tss <- sum((y - mean(y))^2)
  r2 <- if (tss > 0) 1 - rss / tss else 0
  engine <- if (is.list(fit) && !is.null(fit$method)) fit$method else NULL
  lm_fit <- if (identical(engine, "lm")) fit$fit else if (inherits(fit, "lm")) fit else NULL
  k <- if (!is.null(lm_fit)) {
    sum(!is.na(stats::coef(lm_fit))) - 1L
  } else if (!is.null(engine) && engine == "loess") {
    max(1, fit$fit$trace.hat)
  } else if (!is.null(fit$coefficients)) {
    max(1L, sum(abs(fit$coefficients) > sqrt(.Machine$double.eps), na.rm = TRUE))
  } else 1L
  n <- length(y)
  switch(target,
    "adj.r.square" = 1 - (1 - r2) * (n - 1) / max(1, n - k - 1),
    "r.square" = r2,
    "MAE" = mean(abs(y - p)),
    "MAPE" = mean(abs((y - p) / ifelse(y == 0, NA, y)), na.rm = TRUE),
    "MASE" = mean(abs(y - p)) / mean(abs(diff(y))),
    "MSE" = mean((y - p)^2),
    "RMSE" = sqrt(mean((y - p)^2)),
    "AIC" = if (!is.null(lm_fit)) stats::AIC(lm_fit) else n * log(rss / n) + 2 * (k + 1),
    "BIC" = if (!is.null(lm_fit)) stats::BIC(lm_fit) else n * log(rss / n) + log(n) * (k + 1))
}
