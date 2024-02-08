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

dCor<-function(x,y=NULL){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(x)|is.matrix(x)){
    if (min(dim(x))>1){
      if (!is.null(y)){
        warning("x is a matrix or a data.frame with at least two columns, y is neglected.")
      }
      y<-NULL
    }
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
      dCor<-NULL
      stop("Error: x must be a matrix or a dataframe!")
    }
  }else{
    dCor<-energy::dcor(x,y)
    dCor
  }
}

