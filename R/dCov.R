#' Matrix-based distance covariance
#' @export
dCov <- function(x, y = NULL, test = FALSE, adjust = "BH", alpha = NULL,
                 parallel = FALSE, cores = 1L, R = 499L) {
  statistic <- function(a, b) {
    if (length(a) > 50L && exists("dcov2d", envir = asNamespace("energy"),
                                  inherits = FALSE)) {
      sqrt(pmax(energy::dcov2d(a, b, type = "V"), 0))
    } else energy::dcov(a, b)
  }
  test_fun <- function(a, b) energy::dcov.test(a, b, R = R)$p.value
  if (!is.null(y)) {
    keep <- stats::complete.cases(x, y)
    estimate <- statistic(x[keep], y[keep])
    if (!test) return(estimate)
    return(list(r = estimate, p = test_fun(x[keep], y[keep])))
  }
  out <- .nda_pairwise_matrix(x, statistic, test_fun, test, adjust, alpha,
                              parallel, cores)
  m <- as.matrix(x)
  diagonal <- vapply(seq_len(ncol(m)), function(i) statistic(m[, i], m[, i]),
                     numeric(1))
  if (test) {
    diag(out$r) <- diagonal
  } else {
    diag(out) <- diagonal
  }
  out
}
