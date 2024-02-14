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
# Last modified: February 2024                                                #
#-----------------------------------------------------------------------------#
#SUMMARY FUNCTION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#
#' @export
summary.nda <- function(object,  digits =  getOption("digits"), ...) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (methods::is(object,"nda")){
    communality <- object$communality
    loadings <- object$loadings
    uniqueness <- object$uniqueness
    factors <- object$factors
    scores <- object$scores
    n.obs <- object$n.obs
    factors <- object$factors
    if (!is.null(scores)){
      results<-list(cummunality = communality, loadings = loadings,
                    uniqueness = uniqueness,
                    factors = factors,
                    scores = scores,
                    n.obs = n.obs)
    }else{
      results<-list(cummunality = communality, loadings = loadings,
                    uniqueness = uniqueness,
                    factors = factors,
                    n.obs = n.obs)
    }
    return(results)
    print.nda(object)
  }else{
    summary(object,...)
  }
}
