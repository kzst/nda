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


#' KMO-based feature selection
#' @export
fs.KMO <- function(data, min_MSA = 0.5, cor.mtx = FALSE,
                   max_iter = ncol(data) - 1L, ...) {
  x <- as.matrix(.nda_numeric_matrix(data, "data"))
  if (!is.numeric(min_MSA) || length(min_MSA) != 1L || min_MSA < 0 || min_MSA > 1) {
    stop("min_MSA must be between zero and one.", call. = FALSE)
  }
  if (cor.mtx && (nrow(x) != ncol(x) || !isSymmetric(x))) {
    stop("With cor.mtx = TRUE, data must be a symmetric square matrix.", call. = FALSE)
  }
  correlation <- if (cor.mtx) x else stats::cor(x, use = "pairwise.complete.obs")
  removed <- character()
  iter <- 0L
  repeat {
    result <- psych::KMO(correlation)
    if (result$MSA >= min_MSA || ncol(x) <= 2L || iter >= max_iter) break
    drop_index <- which.min(result$MSAi)
    drop_name <- names(result$MSAi)[drop_index]
    if (!length(drop_name) || is.na(drop_name) || !nzchar(drop_name)) {
      drop_name <- colnames(x)[drop_index]
    }
    removed <- c(removed, drop_name)
    correlation <- correlation[-drop_index, -drop_index, drop = FALSE]
    if (cor.mtx) x <- correlation else x <- x[, -drop_index, drop = FALSE]
    iter <- iter + 1L
  }
  out <- if (is.data.frame(data) && !cor.mtx) as.data.frame(x) else x
  attr(out, "KMO") <- result
  attr(out, "removed") <- removed
  out
}
