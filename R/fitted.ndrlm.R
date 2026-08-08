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


#' Fitted values from an NDRLM model
#' @export
fitted.ndrlm <- function(object, ...) {
  if (!inherits(object, "ndrlm")) stop("object must inherit from 'ndrlm'.", call. = FALSE)
  out <- vapply(object$fits, function(fit) fit$fitted, numeric(nrow(object$dep)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  dimnames(out) <- list(rownames(object$dep), names(object$fits))
  as.data.frame(out)
}
