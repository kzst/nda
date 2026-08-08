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


#' Matrix-based semi-partial distance correlation
#' @export
spdCor <- function(x, test = FALSE, adjust = "BH", alpha = NULL,
                   parallel = FALSE, cores = 1L) {
  x <- as.matrix(.nda_numeric_matrix(x))
  cvx <- dCov(x, parallel = parallel, cores = cores)
  inv <- .nda_inverse(cvx)
  residual_precision <- diag(inv) - t(t(inv^2) / diag(inv))
  den <- sqrt(pmax(diag(cvx), .Machine$double.eps)) *
    sqrt(pmax(abs(residual_precision), .Machine$double.eps))
  out <- -stats::cov2cor(inv) / den
  diag(out) <- 1
  dimnames(out) <- list(colnames(x), colnames(x))
  if (!test) return(out)
  .nda_pvalues(out, nrow(x), controls = ncol(x) - 2L,
               adjust = adjust, alpha = alpha)
}
