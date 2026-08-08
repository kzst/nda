#' Matrix-based partial distance correlation
#' @export
pdCor <- function(x, test = FALSE, adjust = "BH", alpha = NULL,
                  parallel = FALSE, cores = 1L) {
  x <- as.matrix(.nda_numeric_matrix(x))
  cvx <- dCov(x, parallel = parallel, cores = cores)
  inv <- .nda_inverse(cvx)
  den <- sqrt(pmax(abs(diag(inv)), .Machine$double.eps))
  out <- -inv / outer(den, den)
  diag(out) <- 1
  dimnames(out) <- list(colnames(x), colnames(x))
  if (!test) return(out)
  .nda_pvalues(out, nrow(x), controls = ncol(x) - 2L,
               adjust = adjust, alpha = alpha)
}
