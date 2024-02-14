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
######### Normalize entire data, row, or column #######
#' @export
normalize <- function(x,type="all")
{
  results<-NULL
  if ((is.data.frame(x))|(is.matrix(x))|
      (is.array(x))){
    results<-((x - min(x)) / (max(x) - min(x)))
    if ("row" %in% type){
      for (i in 1:nrow(x)){
        results[i,]<-((x[i,] - min(x[i,])) / (max(x[i,]) - min(x[i,])))
      }
    }else{
      if ("col" %in% type){
        for (i in 1:ncol(x)){
          results[,i]<-((x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i])))
        }
      }
    }
  }else{
    stop(
      "Only matrix, array or data.frame can be used in this function!",
      call. = FALSE
    )
  }
  return(results)
}
