#-----------------------------------------------------------------------------#
#                                                                             #
#  NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)                  #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona       #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: October 2022                                                 #
#-----------------------------------------------------------------------------#

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
<<<<<<< HEAD
      loop=TRUE
      while(loop==TRUE){
        kmo<-psych::KMO(x)
        if (min(kmo$MSAi)>min_MSA){loop=FALSE}else{
=======
      repeat{
        kmo<-psych::KMO(x)
        if (min(kmo$MSAi)>min_MSA){break}else{
>>>>>>> f16cf31a08b36536137b862845015c4923264d69
          i<-which.min(kmo$MSAi)
          if (cor.mtx==TRUE){
            x<-x[-i,-i]
          }else{
            x<-x[,-i]
          }
        }
<<<<<<< HEAD
        if (ncol(x)<2){
          loop=FALSE
        }
=======
        if (ncol(x)<2){break}
>>>>>>> f16cf31a08b36536137b862845015c4923264d69
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
