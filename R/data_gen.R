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
#DATA GENERATION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#
#' @export
data_gen<-function(n,m,nfactors=2,lambda=1){
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  M<-NA
  if (n>=1)
    {
      if (m>=1)
        {
          M<-matrix(0,nrow=n,ncol=m)
            if (nfactors>=1)
              {
                L<-replicate(nfactors,matrix(1,ceiling(n/nfactors),
                                    ceiling(m/nfactors)),simplify=FALSE)
                M<-Matrix::bdiag(L)
                M<-as.matrix(M[1:n,1:m])
                N<-matrix(stats::runif(n*m),n,m)
                M<-M-N*M/exp(lambda)
              }
            else
              {
                warning("nfactors must be equal to or greater than 1!")
              }
      }
        else
          {
            warning("m must be equal to or greater than 1!")
          }
    }
      else
        {
          warning("n must be equal to or greater than 1!")
        }
  return(as.data.frame(M))
}
