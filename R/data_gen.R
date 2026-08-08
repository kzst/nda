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


#' Generate structured data for NDA examples
#' @export
data_gen <- function(n, m, nfactors = 2, lambda = 1, sparse = FALSE,
                     density = 1, seed = NULL) {
  args <- c(n = n, m = m, nfactors = nfactors)
  if (any(!is.finite(args)) || any(args < 1) || any(args != floor(args))) {
    stop("n, m, and nfactors must be positive integers.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("lambda must be a finite number.", call. = FALSE)
  }
  if (!is.numeric(density) || length(density) != 1L || density <= 0 || density > 1) {
    stop("density must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  row_group <- pmin(nfactors, ceiling(seq_len(n) * nfactors / n))
  col_group <- pmin(nfactors, ceiling(seq_len(m) * nfactors / m))
  columns_by_group <- split(seq_len(m), col_group)
  columns <- unlist(columns_by_group[as.character(row_group)], use.names = FALSE)
  rows <- rep(seq_len(n), lengths(columns_by_group[as.character(row_group)]))
  if (density < 1) {
    keep <- stats::runif(length(rows)) <= density
    rows <- rows[keep]
    columns <- columns[keep]
  }
  generated <- 1 - stats::runif(length(rows)) / exp(lambda)
  names_out <- list(paste0("case", seq_len(n)), paste0("V", seq_len(m)))
  if (sparse) {
    Matrix::sparseMatrix(i = rows, j = columns, x = generated,
                         dims = c(n, m), dimnames = names_out)
  } else {
    values <- matrix(0, n, m, dimnames = names_out)
    values[cbind(rows, columns)] <- generated
    as.data.frame(values)
  }
}
