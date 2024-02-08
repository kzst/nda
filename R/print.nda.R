#-----------------------------------------------------------------------------#
#                                                                             #
#  GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (GNDA)     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona       #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: February 2023                                                #
#-----------------------------------------------------------------------------#
#' @export
print.nda <- function(x,  digits =  getOption("digits"), ...) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (methods::is(x,"nda")){
    communality <- x$communality
    loadings <- x$loadings
    uniqueness <- x$uniqueness
    factors <- x$factors
    scores <- x$scores
    n.obs <- x$n.obs
    factors <- x$factors
    cat("\nPrint of the NDA:\n")
    cat("\nNumber of latent variables: ",factors)
    cat("\nNumber of observations: ",n.obs)
    cat("\nCommunalities:\n")
    print(communality,digits = digits, ...)
    cat("\nFactor loadings:\n")
    print(loadings,digits = digits, ...)
    if (!is.null(scores)){
      cat("\nFactor scores:\n")
      print(scores,digits = digits, ...)
      cat("\n\nCorrelation matrix of factor scores:\n")
      print(stats::cor(scores),digits = digits, ...)
    }
  }else{
    print(x,...)
  }
}
