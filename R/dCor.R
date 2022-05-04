#-----------------------------------------------------------------------------#
#                                                                             #
#  NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)                  #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona       #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: May 2022                                                     #
#-----------------------------------------------------------------------------#

#' @export

dCor<-function(x,y=NULL){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.null(y)){
    if (is.data.frame(x)|is.matrix(x)){
      dC<-matrix(0,nrow=ncol(x),ncol=ncol(x))
      for (i in c(1:ncol(x))){
        for (j in c(1:ncol(x))){
          dC[i,j]<-energy::dcor(x[,i],x[,j])
        }
      }
      rownames(dC)<-colnames(x)
      colnames(dC)<-colnames(x)
      dCor<-dC
      dCor
    }else{
      stop("Error: x must be a matrix or a dataframe!")
      dCor<-NULL
    }
  }else{
    dCor<-energy::dcor(x,y)
    dCor
  }
}

