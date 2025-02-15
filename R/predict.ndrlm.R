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
predict.ndrlm <- function(object,  newdata, se.fit = FALSE, scale = NULL, df = Inf,
                          interval = c("none", "confidence", "prediction"),
                          level = 0.95, type = c("response", "terms"),
                          terms = NULL, na.action = stats::na.pass,
                          pred.var = 1/weights, weights = 1,...) {
  if (methods::is(object,"ndrlm")){
    Call<-object$Call
    fval<-object$fval
    seed<-object$seed
    if (!is.null(seed)){
      set.seed(seed)
    }
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
    prediction<-list()

    if (is.null(newdata)){
      for (i in 1:length(fits)){
        prediction[[i]]<-stats::fitted(fits[[i]])
      }
    }else{
      for (i in 1:length(fits)){
        newdata<-as.data.frame(newdata)
        newdata.X<-as.data.frame(newdata[,colnames(X)])
        colnames(newdata.X)<-colnames(X)
        rownames(newdata.X)<-rownames(newdata)
        newdata.indep<-newdata.X
        newdata.Y<-newdata[,colnames(Y)]
        newdata.Y<-as.data.frame(newdata[,colnames(Y)])
        colnames(newdata.Y)<-colnames(Y)
        rownames(newdata.Y)<-rownames(newdata)
        newdata.dep<-newdata.Y
        if (latents %in% c("in","both")){
          newdata.NDAin<-predict.nda(NDAin,newdata.X)
          colnames(newdata.NDAin)<-paste("NDAin",1:NDAin$factors,sep="")
          newdata.indep<-cbind(newdata.X,newdata.NDAin)[,colnames(indep)]
          newdata.indep<-as.data.frame(newdata.indep)
          colnames(newdata.indep)<-colnames(indep)
          rownames(newdata.indep)<-rownames(newdata)
        }
        if (latents %in% c("out","both")){
          newdata.NDAout<-predict.nda(NDAout,newdata.Y)
          colnames(newdata.NDAout)<-paste("NDAout",1:NDAout$factors,sep="")
          newdata.dep<-cbind(newdata.Y,newdata.NDAout)[,colnames(dep)]
          newdata.dep<-as.data.frame(newdata.dep)
          colnames(newdata.dep)<-colnames(dep)
          rownames(newdata.dep)<-rownames(newdata)
        }

        newdata.final<-cbind(newdata.dep[,i],newdata.indep)
        colnames(newdata.final)[1]<-colnames(newdata.dep)[i]
        colnames(newdata.final)[-1]<-colnames(newdata.indep)
        newdata.final<-as.data.frame(newdata.final)
        rownames(newdata.final)<-rownames(newdata)
        prediction[[i]]<-stats::predict.lm(fits[[i]],newdata = newdata.final,
                                    se.fit = se.fit, scale = scale, df = df,
                                    interval = interval,
                                    level = level, type = type,
                                    terms = terms, na.action = na.action,
                                    pred.var = pred.var, weights = weights)
      }
    }
    return(prediction)
  }
}
