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
### PREDICT SCORES NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ##
#' @export
predict.nda <- function(object,  newdata,...) {
  if (methods::is(object,"nda")){
    Call<-object$Call
    LOADING<-object$loadings
    SCORES<-object$scores
    EVCs<-object$EVCs
    center<-object$center
    membership<-object$membership
    weight<-object$weight
    factors<-object$factors
    use_rotation<-object$use_rotation
    rotation<-object$rotation
    seed<-object$seed
    if (!is.null(seed)){
      set.seed(seed)
    }
    if (length(membership)!=ncol(newdata)){
      stop(
        "The columns of newdata and the original date must be same.",
        call. = FALSE
      )
    }

    Coords<-1:length(membership)
    L<-as.data.frame(matrix(0,nrow = nrow(newdata),ncol=factors))
    colnames(L)<-colnames(SCORES)
    rownames(L)<-rownames(newdata)
    if (is.null(weight)){
      weight=rep(1,ncol(r))
    }
    r<-t(t(newdata)*weight)
    DATA<-r
    X<-r

    for (i in 1:factors){
      EVC<-EVCs[[i]]
      Coordsi<-Coords[membership==i]
      result<-NA
      try(result <- as.matrix(rowSums(r[,Coordsi] %*% EVC)),silent=TRUE)
      if (is.null(nrow(is.nan(result)))){
        try(result <- as.matrix(rowSums(r[,Coordsi] * EVC)),silent=TRUE)
      }
      L[,i]<-result
    }
    if (ncol(L)>1 && use_rotation==TRUE){
      L<-psych::principal(L,nfactors = dim(L)[2],
                          rotate = rotation)$scores
    }else{
      L<-scale(L,center = center)
    }
    return(L)
  }
}
