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

spdCor<-function(x){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package \"MASS\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))

  # sample number
  n <- dim(x)[1]

  # given variables' number
  gp <- dim(x)[2]-2

  # covariance matrix
  cvx <- dCov(x)

  # inverse covariance matrix
  if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
    icvx <- MASS::ginv(cvx)
  }else
    icvx <- solve(cvx)

  # semi-partial correlation
  spcor <- -stats::cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
  diag(spcor) <- 1
  spdCor<-spcor
  spdCor
}

