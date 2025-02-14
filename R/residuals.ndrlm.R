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
# RESIDUALS FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) #
#' @export
residuals.ndrlm <- function(object,  ...) {
  if (methods::is(object,"ndrlm")){
    Call<-object$Call
    fval<-object$fval
    pareto<-object$pareto
    X<-object$X
    Y<-object$Y
    latents<-object$latents
    if (latents %in% c("in","both")){
      NDAin<-object$NDAin
      NDAin_weight<-object$NDAin_weight
      NDAin_min_evalue<-object$NDAin_min_evalue
      NDAin_min_communality<-object$NDAin_min_communality
      NDAin_com_communalities<-object$NDAin_com_communalities
      NDAin_min_R<-object$NDAin_com_communalities
    }
    if (latents %in% c("out","both")){
      NDAout<-object$NDAout
      NDAout_weight<-object$NDAout_weight
      NDAout_min_evalue<-object$NDAout_min_evalue
      NDAout_min_communality<-object$NDAout_min_communality
      NDAout_com_communalities<-object$NDAout_com_communalities
      NDAout_min_R<-object$NDAout_com_communalities
    }
    fits<-object$fits
    optimized<-object$optimized
    if (optimized==TRUE){
      NSGA<-object$NSGA
    }
    extra_vars.X<-object$extra_vars.X
    extra_vars.Y<-object$extra_vars.Y
    if (latents %in% c("in","both")){
      if (extra_vars.X==TRUE){
        dircon_X<-object$dircon_X
      }
    }
    if (latents %in% c("out","both")){
      if (extra_vars.Y==TRUE){
        dircon_Y<-object$dircon_Y
      }
    }
    fn<-object$fn
    dep<-Y
    if (latents %in% c("out","both")){
      if (extra_vars.Y==TRUE){
        dep<-cbind(NDAout$scores,Y[,NDAout$membership==0])
        dep<-as.data.frame(dep)
        colnames(dep)<-c(paste("NDAout",1:NDAout$factors,sep=""),
                         colnames(Y)[NDAout$membership==0])
      }else{
        dep<-NDAout$scores
        colnames(dep)<-paste("NDAout",1:NDAout$factors,sep="")
      }
    }
    indep<-X
    if (latents %in% c("in","both")){
      if (extra_vars.X==TRUE){
        indep<-cbind(NDAin$scores,X[,NDAin$membership==0])
        indep<-as.data.frame(indep)
        colnames(indep)<-c(paste("NDAin",1:NDAin$factors,sep=""),
                           colnames(X)[NDAin$membership==0])
      }else{
        indep<-NDAin$scores
        colnames(indep)<-paste("NDAin",1:NDAin$factors,sep="")
      }
    }
    RESIDUALS<-as.data.frame(matrix(0,nrow=nrow(dep),ncol=ncol(dep)))
    colnames(RESIDUALS)<-colnames(dep)
    rownames(RESIDUALS)<-rownames(dep)
    for (i in 1:length(fits)){
      RESIDUALS[,i]<-stats::fitted(fits[[i]])
    }
    RESIDUALS<-RESIDUALS[,1:length(fits)]
    return(RESIDUALS)
  }
}
