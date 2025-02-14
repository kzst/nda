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
# Last modified: December 2024                                                #
#-----------------------------------------------------------------------------#
## PRINT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
#' @export
print.ndrlm <- function(x, digits = getOption("digits"), ...) {
  if (methods::is(x,"ndrlm")){
    Call<-x$Call
    target<-x$target
    fval<-x$fval
    pareto<-x$pareto
    X<-x$X
    Y<-x$Y
    latents<-x$latents
    if (latents %in% c("in","both")){
      NDAin<-x$NDAin
      NDAin_weight<-x$NDAin_weight
      NDAin_min_evalue<-x$NDAin_min_evalue
      NDAin_min_communality<-x$NDAin_min_communality
      NDAin_com_communalities<-x$NDAin_com_communalities
      NDAin_min_R<-x$NDAin_min_R
    }
    if (latents %in% c("out","both")){
      NDAout<-x$NDAout
      NDAout_weight<-x$NDAout_weight
      NDAout_min_evalue<-x$NDAout_min_evalue
      NDAout_min_communality<-x$NDAout_min_communality
      NDAout_com_communalities<-x$NDAout_com_communalities
      NDAout_min_R<-x$NDAout_min_R
    }
    fits<-x$fits
    optimized<-x$optimized
    if (optimized==TRUE){
      NSGA<-x$NSGA
    }
    extra_vars.X<-x$extra_vars.X
    extra_vars.Y<-x$extra_vars.Y
    if (latents %in% c("in","both")){
      if (extra_vars.X==TRUE){
        dircon_X<-x$dircon_X
      }
    }
    if (latents %in% c("out","both")){
      if (extra_vars.Y==TRUE){
        dircon_Y<-x$dircon_Y
      }
    }
    fn<-x$fn
    cat("\nBrief summary of NDRLM:\n")
    cat("\nFunction call: ")
    print(Call)

    cat("\nNumber of independent variables: ",ncol(X))
    cat("\nNumber of dependent variables: ",ncol(Y))
    if (latents %in% c("in","both")){
      cat("\nNumber of latent-independent variables: ",NDAin$factors)
    }
    if (latents %in% c("out","both")){
      cat("\nNumber of latent-dependent variables: ",NDAout$factors)
    }
    if (latents %in% c("in","both")){
      if (extra_vars.X==TRUE){
        cat("\nNumber of dropped independent variables: ",sum((NDAin$membership==0)))
      }
    }
    if (latents %in% c("out","both")){
      if (extra_vars.Y==TRUE){
        cat("\nNumber of dropped dependent variables: ",sum((NDAout$membership==0)))
      }
    }
    if (latents %in% c("in","both")){
      cat("\n\nSummary of dimensionality reduction for independent variables\n")
      print.nda(NDAin,digits = digits)
    }
    if (latents %in% c("out","both")){
      cat("\n\nSummary of dimensionality reduction for dependent variables\n")
      print.nda(NDAout,digits = digits)
    }
    cat("\n\nSummary of fitting\n")
    if (optimized==TRUE){
      cat("\nOptimized fittings\n")
      cat("\nTarget performance measure: ",target)
    }else{
      cat("\nNon-optimized fittings\n")
    }

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

    cat("\nList of dependent variables: ",toString(colnames(dep)))
    cat("\nList of independent variables: ",toString(colnames(indep)))
    if (latents %in% c("in","both")){
      cat("\nList of latent-independent variables: ",toString(paste("NDAin",1:NDAin$factors,sep="")))
      if (extra_vars.X==TRUE){
        cat("\nList of non-groupped independent variables: ",toString(dircon_X))
      }
    }
    if (latents %in% c("out","both")){
      cat("\nList of latent-dependent variables: ",toString(paste("NDAout",1:NDAout$factors,sep="")))
      if (extra_vars.Y==TRUE){
        cat("\nList of non-groupped independent variables: ",toString(dircon_Y))
      }
    }

    for (i in 1:length(fits)){
      cat("\nFitting for variable ",colnames(fits[[i]]$model)[1])
      print(lm.beta::summary.lm.beta(lm.beta::lm.beta(fits[[i]])))
    }
  }
}
