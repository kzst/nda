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


#' Matrix-based distance correlation
#' @export
dCor <- function(x, y = NULL, test = FALSE, adjust = "BH", alpha = NULL,
                 parallel = FALSE, cores = 1L, R = 499L) {
  statistic <- function(a, b) {
    if (length(a) > 50L && exists("dcor2d", envir = asNamespace("energy"),
                                  inherits = FALSE)) {
      sqrt(pmax(energy::dcor2d(a, b, type = "V"), 0))
    } else energy::dcor(a, b)
  }
  test_fun <- function(a, b) energy::dcor.test(a, b, R = R)$p.value
  if (!is.null(y)) {
    keep <- stats::complete.cases(x, y)
    estimate <- statistic(x[keep], y[keep])
    if (!test) return(estimate)
    return(list(r = estimate, p = test_fun(x[keep], y[keep])))
  }
  .nda_pairwise_matrix(x, statistic, test_fun, test, adjust, alpha,
                       parallel, cores)
}
