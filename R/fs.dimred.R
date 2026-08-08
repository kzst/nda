#' Feature selection from dimensionality-reduction loadings
#' @export
fs.dimred <- function(fn, DF, min_comm = 0.25, com_comm = 0.25, ...) {
  data <- as.matrix(.nda_numeric_matrix(DF, "DF"))
  loadings <- if (inherits(fn, "nda")) {
    fn$loadings
  } else if (inherits(fn, "prcomp")) {
    fn$rotation
  } else if (inherits(fn, "princomp")) {
    unclass(fn$loadings)
  } else if (!is.null(fn$loadings)) {
    unclass(fn$loadings)
  } else if (!is.null(fn$rotation)) {
    fn$rotation
  } else if (!is.null(fn$var$coord)) {
    fn$var$coord
  } else {
    stop("fn does not expose loadings, a rotation matrix, or variable coordinates.",
         call. = FALSE)
  }
  loadings <- as.matrix(loadings)
  if (is.null(rownames(loadings))) rownames(loadings) <- colnames(data)[seq_len(nrow(loadings))]
  squared <- loadings^2
  communality <- if (!is.null(fn$communality)) as.numeric(fn$communality) else
    apply(squared, 1L, max, na.rm = TRUE)
  names(communality) <- rownames(loadings)
  low <- names(communality)[communality < min_comm]
  common <- character()
  if (ncol(squared) > 1L) {
    ordered <- t(apply(squared, 1L, sort, decreasing = TRUE))
    common <- rownames(loadings)[ordered[, 1L] < ordered[, 2L] + com_comm |
                                  ordered[, 1L] <= 2 * ordered[, 2L]]
  }
  common <- setdiff(common, low)
  removed <- union(low, common)
  correlation_input <- nrow(data) == ncol(data) && isSymmetric(data) &&
    setequal(rownames(data), colnames(data))
  if (correlation_input) {
    retained <- data[setdiff(colnames(data), removed),
                     setdiff(colnames(data), removed), drop = FALSE]
    dropped_low <- dropped_com <- NULL
  } else {
    retained <- data[, setdiff(colnames(data), removed), drop = FALSE]
    dropped_low <- if (length(low)) data[, low, drop = FALSE] else NULL
    dropped_com <- if (length(common)) data[, common, drop = FALSE] else NULL
    if (is.data.frame(DF)) {
      retained <- as.data.frame(retained)
      if (!is.null(dropped_low)) dropped_low <- as.data.frame(dropped_low)
      if (!is.null(dropped_com)) dropped_com <- as.data.frame(dropped_com)
    }
  }
  list(dropped_low = dropped_low, dropped_com = dropped_com,
       retained_DF = retained, remain_DF = retained,
       communality = communality, loadings = loadings, model = fn)
}
