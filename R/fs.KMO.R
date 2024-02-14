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
######### Feature selection for KMO #######
#' @export
fs.KMO<-function(data,min_MSA=0.5,cor.mtx=FALSE){
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(data)|is.matrix(data)){
    if (ncol(data)>=2){
      x<-data
      loop=TRUE
      while(loop==TRUE){
        kmo<-psych::KMO(x)
        if (min(kmo$MSAi)>min_MSA){loop=FALSE}else{
          i<-which.min(kmo$MSAi)
          if (cor.mtx==TRUE){
            x<-x[-i,-i]
          }else{
            x<-x[,-i]
          }
        }
        if (ncol(x)<2){
          loop=FALSE
        }
      }
      return(x)
    }else{
      stop("Error: data must contain at least 2 columns!")
      step.KMO<-NULL
      return(step.KMO)
    }
  }else{
    stop("Error: data must be a matrix or a dataframe!")
    step.KMO<-NULL
    return(step.KMO)
  }
}
