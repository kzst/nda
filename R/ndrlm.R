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
# Last modified: November 2024                                                #
#-----------------------------------------------------------------------------#
### GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (GNDR) ##
#' @export
ndrlm<-function(Y,X,latents="in",dircon=FALSE,optimize=TRUE,
                target="adj.r.square",rel_weight=FALSE,
                cor_method=1,
                cor_type=1,min_comm=2,Gamma=1,
                null_model_type=4,mod_mode=1,use_rotation=FALSE,
                rotation="oblimin",pareto=FALSE,fit_weights=NULL,
                lower.bounds.x = c(rep(-100,ncol(X))),
                upper.bounds.x = c(rep(100,ncol(X))),
                lower.bounds.latentx = c(0,0,0,0),
                upper.bounds.latentx = c(0.6,0.6,0.6,0.3),
                lower.bounds.y = c(rep(-100,ncol(Y))),
                upper.bounds.y = c(rep(100,ncol(Y))),
                lower.bounds.latenty = c(0,0,0,0),
                upper.bounds.latenty = c(0.6,0.6,0.6,0.3),
                popsize = 20, generations = 30, cprob = 0.7, cdist = 5,
                mprob = 0.2, mdist=10, seed=NULL){
  cl<-match.call()
  if (!requireNamespace("mco", quietly = TRUE)) {
    stop(
      "Package \"mco\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!(target %in% c("adj.r.square","r.sqauare","MAE","MAPE","MASE","MSE","RMSE"))){
    stop(
      "Target must be either adj.r.square, r.sqauare, MAE, MAPE, MASE, MSE, or RMSE",
      call. = FALSE
    )
  }
  ERROR<-TRUE
  Y<-as.data.frame(Y)
  X<-as.data.frame(X)
  extra_vars.X=dircon
  extra_vars.Y=dircon
  weight.X<-rep(1,ncol(X))
  weight.Y<-rep(1,ncol(Y))
  latent.X<-c(0,0,0,0)
  latent.Y<-c(0,0,0,0)
  if (rel_weight){
    lower.bounds.x[lower.bounds.x<0]<-0
    lower.bounds.y[lower.bounds.y<0]<-0
    upper.bounds.x[upper.bounds.x<0]<-0
    upper.bounds.y[upper.bounds.y<0]<-0
  }
  if (("in" %in% latents)==FALSE){ # Pareto-optimiality can be found,
    pareto=FALSE         #  if there are only latent-independent variables
  }
  if (("none" %in% latents)==TRUE){ # If there are no latent variables,
    optimize=FALSE # there is no way to optimize fittings.
  }

  hyperparams=switch(
    latents,
    "in"=c(weight.X,latent.X),
    "out"=c(weight.Y,latent.Y),
    "both"=c(weight.X,latent.X,weight.Y,latent.Y),
    "none"=NULL
  )
  lower.bounds=switch(
    latents,
    "in"=c(lower.bounds.x,lower.bounds.latentx),
    "out"=c(lower.bounds.y,lower.bounds.latenty),
    "both"=c(lower.bounds.x,lower.bounds.latentx,
             lower.bounds.y,lower.bounds.latenty),
    "none"=NULL
  )
  upper.bounds=switch(
    latents,
    "in"=c(upper.bounds.x,upper.bounds.latentx),
    "out"=c(upper.bounds.y,upper.bounds.latenty),
    "both"=c(upper.bounds.x,upper.bounds.latentx,
             upper.bounds.y,upper.bounds.latenty),
    "none"=NULL
  )
  tmp_hyper<-hyperparams
  cost<-function(hyperparams){ # Cost function
    hyperparams[is.na(hyperparams)]<-0
    if ("in" %in% latents){
      weight.X<-hyperparams[1:ncol(X)]
      params.X<-hyperparams[-c(1:ncol(X))]
    }else{
      if ("out" %in% latents){
        weight.Y<-hyperparams[1:ncol(Y)]
        params.Y<-hyperparams[-c(1:ncol(Y))]
      }else{
        if ("both" %in% latents){
          weight.X<-hyperparams[1:ncol(X)]
          params.X<-hyperparams[(ncol(X)+1):(ncol(X)+4)]
          weight.Y<-hyperparams[(ncol(X)+5):(ncol(X)+4+ncol(Y))]
          params.Y<-hyperparams[(ncol(X)+5+ncol(Y)):(ncol(X)+4+ncol(Y)+4)]
        }
      }
    }
    if (latents %in% c("in","both")){ # For latent-independent variables
      NDA_in<-try(ndr(X,min_evalue = params.X[1],
                      min_communality=params.X[2],
                      com_communalities = params.X[3],
                      min_R = params.X[4],weight=weight.X,covar=FALSE,
                      cor_method=cor_method,
                      cor_type=cor_type,
                      min_comm=min_comm,
                      Gamma=Gamma,
                      null_model_type=null_model_type,
                      mod_mode=mod_mode,
                      use_rotation=use_rotation,
                      rotation=rotation),silent=TRUE)
    }

    if (latents %in% c("out","both")){ # For latent-dependent variables
      NDA_out<-try(ndr(Y,min_evalue = params.Y[1],
                       min_communality=params.Y[2],
                       com_communalities = params.Y[3],
                       min_R = params.Y[4],weight=weight.Y,covar=FALSE,
                       cor_method=cor_method,
                       cor_type=cor_type,
                       min_comm=min_comm,
                       Gamma=Gamma,
                       null_model_type=null_model_type,
                       mod_mode=mod_mode,
                       use_rotation=use_rotation,
                       rotation=rotation),silent=TRUE)
    }
    errorvalue<-switch(target,
                       "adj.r.square" = 0,
                       "r.sqauare" = 0,
                       "MAE" =Inf,
                       "MAPE" = Inf,
                       "MASE" = Inf,
                       "MSE" = Inf,
                       "RMSE" = Inf
    )
    if (pareto==TRUE){
      res<-rep(errorvalue,ncol(Y))
    }else{
      res<-errorvalue
    }
    error<-FALSE
    if (latents %in% c("out","both")){
      if (inherits(NDA_out,"try-error")){
        error<-TRUE
        return(res)
      }
    }

    if (latents %in% c("in","both")){
      if (inherits(NDA_in,"try-error")){
        error<-TRUE
        return(res)
      }
    }
    if (error==FALSE){
      extra_vars.X<-FALSE
      extra_vars.Y<-FALSE
      dropped_X<-NULL
      dropped_Y<-NULL
      if (latents %in% c("in")){
        if (!inherits(NDA_in,"try-error")){
          if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
            extra_vars.X<-TRUE
            dropped_X<-X[,NDA_in$membership==0]
          }
        }else{return(res)}
      }else{
        if (latents %in% c("out")){
          if (!inherits(NDA_out,"try-error")){
            if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
              extra_vars.Y<-TRUE
              dropped_Y<-Y[,NDA_out$membership==0]
            }
          }else{return(res)}
        }else{
          if (latents %in% c("both")){
            if (!inherits(NDA_in,"try-error")){
              if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
                extra_vars.X<-TRUE
                dropped_X<-X[,NDA_in$membership==0]
              }
            }else{return(res)}
            if (!inherits(NDA_out,"try-error")){
              if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
                extra_vars.Y<-TRUE
                dropped_Y<-Y[,NDA_out$membership==0]
              }
            }else{return(res)}
          }
        }
      }
      dep<-Y
      if (latents %in% c("out","both")){
        if (extra_vars.Y==TRUE){
          dep<-cbind(as.data.frame(NDA_out$scores),dropped_Y)
          dep<-as.data.frame(dep)
          colnames(dep)<-c(paste("NDAout",1:NDA_out$factors,sep=""),
                           colnames(Y)[NDA_out$membership==0])
        }else{
          dep<-as.data.frame(NDA_out$scores)
          colnames(dep)<-paste("NDAout",1:NDA_out$factors,sep="")
        }
      }
      indep<-X
      if (latents %in% c("in","both")){
        if (extra_vars.X==TRUE){
          indep<-cbind(as.data.frame(NDA_in$scores),dropped_X)
          indep<-as.data.frame(indep)
          colnames(indep)<-c(paste("NDAin",1:NDA_in$factors,sep=""),
                             colnames(X)[NDA_in$membership==0])
        }else{
          indep<-as.data.frame(NDA_in$scores)
          colnames(indep)<-paste("NDAin",1:NDA_in$factors,sep="")
        }
      }

      for (i in 1:ncol(dep))
      {
        data<-cbind(dep[,i],indep)
        colnames(data)[1]<-colnames(dep)[i]
        colnames(data)[-1]<-colnames(indep)
        data<-as.data.frame(data)
        fit<-stats::lm(str2lang(paste(paste("`",colnames(data)[1],"`",sep=""),"~",
                               gsub(", ","+",
                                    toString(
                                      paste('`',
                                            colnames(data)[-1],'`',sep=""))))),data)
        res[i]<-switch(target,
                       "adj.r.square" = stats::summary.lm(fit)$adj.r.squared,
                       "r.sqauare" = stats::summary.lm(fit)$r.squared,
                       "MAE" = Metrics::mae(data[,1],stats::fitted(fit)),
                       "MAPE" = Metrics::mape(data[,1],stats::fitted(fit)),
                       "MASE" = Metrics::mase(data[,1],stats::fitted(fit)),
                       "MSE" = Metrics::mse(data[,1],stats::fitted(fit)),
                       "RMSE" = Metrics::rmse(data[,1],stats::fitted(fit))
        )
      }
      if (pareto==TRUE){
        return(res)
      }else{
        if (is.null(fit_weights)){
          return(mean(res))
        }else{
          return(stats::weighted.mean(res,fit_weights))
        }
      }
    }
  }
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  if (target %in% c("adj.r.square","r.square")){
    costmin <- function(hyperparams) -cost(hyperparams)
  }else{
    costmin <- function(hyperparams) cost(hyperparams)
  }

  if (pareto==TRUE){
    ODIM<-ncol(Y)
  }else{
    ODIM<-1
  }
  if ((optimize==TRUE)&&(latents %in% c("in","out","both"))){
    NSGA <- mco::nsga2(fn=costmin,idim=length(hyperparams),odim=ODIM,
                       lower.bounds = lower.bounds,
                       upper.bounds = upper.bounds,
                       popsize = popsize,
                       generations = generations, cprob = cprob, cdist = cdist,
                       mprob = mprob, mdist=mdist,vectorized = FALSE)
    HYPERPARAMS<-NSGA$par
    HYPERPARAMS[is.na(HYPERPARAMS)]<-0
    I<-1 ## Try to find a feasible solution
    while (I<=nrow(HYPERPARAMS)&&(ERROR==TRUE)){
      hyperparams<-HYPERPARAMS[I,]
      hyperparams[is.na(hyperparams)]<-0
      if ("in" %in% latents){
        weight.X<-hyperparams[1:ncol(X)]
        if (rel_weight){
          weight.X<-weight.X*ncol(X)/sum(weight.X)
        }
        params.X<-hyperparams[-c(1:ncol(X))]
      }else{
        if ("out" %in% latents){
          weight.Y<-hyperparams[1:ncol(Y)]
          if (rel_weight){
            weight.Y<-weight.Y*ncol(Y)/sum(weight.Y)
          }
          params.Y<-hyperparams[-c(1:ncol(Y))]
        }else{
          if ("both" %in% latents){
            weight.X<-hyperparams[1:ncol(X)]
            params.X<-hyperparams[(ncol(X)+1):(ncol(X)+4)]
            weight.Y<-hyperparams[(ncol(X)+5):(ncol(X)+4+ncol(Y))]
            params.Y<-hyperparams[(ncol(X)+5+ncol(Y)):(ncol(X)+4+ncol(Y)+4)]
            if (rel_weight){
              weight.X<-weight.X*ncol(X)/sum(weight.X)
              weight.Y<-weight.Y*ncol(Y)/sum(weight.Y)
            }
          }
        }
      }
      error<-FALSE
      if (latents %in% c("in","both")){ # For latent-independent variables
        NDA_in<-try(ndr(X,min_evalue = params.X[1],
                        min_communality=params.X[2],
                        com_communalities = params.X[3],
                        min_R = params.X[4],weight=weight.X,covar=FALSE,
                        cor_method=cor_method,
                        cor_type=cor_type,
                        min_comm=min_comm,
                        Gamma=Gamma,
                        null_model_type=null_model_type,
                        mod_mode=mod_mode,
                        use_rotation=use_rotation,
                        rotation=rotation),silent=TRUE)
        if (inherits(NDA_in,"try-error")){
          error<-TRUE
        }
      }

      if (latents %in% c("out","both")){ # For latent-dependent variables
        NDA_out<-try(ndr(Y,min_evalue = params.Y[1],
                         min_communality=params.Y[2],
                         com_communalities = params.Y[3],
                         min_R = params.Y[4],weight=weight.Y,covar=FALSE,
                         cor_method=cor_method,
                         cor_type=cor_type,
                         min_comm=min_comm,
                         Gamma=Gamma,
                         null_model_type=null_model_type,
                         mod_mode=mod_mode,
                         use_rotation=use_rotation,
                         rotation=rotation),silent=TRUE)
        if (inherits(NDA_in,"try-error")){
          error<-TRUE
        }
      }
      if (error==FALSE){
        ERROR<-FALSE
        fits<-list()
        extra_vars.X<-FALSE
        extra_vars.Y<-FALSE
        dropped_X<-NULL
        dropped_Y<-NULL
        if (latents %in% c("in")){
          if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
            extra_vars.X<-TRUE
            dropped_X<-X[,NDA_in$membership==0]
          }
        }else{
          if (latents %in% c("out")){
            if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
              extra_vars.Y<-TRUE
              dropped_Y<-Y[,NDA_out$membership==0]
            }
          }else{
            if (latents %in% c("both")){
              if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
                extra_vars.X<-TRUE
                dropped_X<-X[,NDA_in$membership==0]
              }
              if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
                extra_vars.Y<-TRUE
                dropped_Y<-Y[,NDA_out$membership==0]
              }
            }
          }
        }

        dep<-Y
        if (latents %in% c("out","both")){
          if ((extra_vars.Y==TRUE)&&(!is.null(dropped_Y))){
            dep<-cbind(as.data.frame(NDA_out$scores),dropped_Y)
            dep<-as.data.frame(dep)
            colnames(dep)<-c(paste("NDAout",1:NDA_out$factors,sep=""),
                             colnames(Y)[NDA_out$membership==0])
          }else{
            dep<-as.data.frame(NDA_out$scores)
            colnames(dep)<-paste("NDAout",1:NDA_out$factors,sep="")
          }
        }
        indep<-X
        if (latents %in% c("in","both")){
          if ((extra_vars.X==TRUE)&&(!is.null(dropped_X))){
            indep<-cbind(as.data.frame(NDA_in$scores),dropped_X)
            indep<-as.data.frame(indep)
            colnames(indep)<-c(paste("NDAin",1:NDA_in$factors,sep=""),
                               colnames(X)[NDA_in$membership==0])
          }else{
            indep<-as.data.frame(NDA_in$scores)
            colnames(indep)<-paste("NDAin",1:NDA_in$factors,sep="")
          }
        }

        for (i in 1:ncol(dep))
        {
          data<-cbind(dep[,i],indep)
          colnames(data)[1]<-colnames(dep)[i]
          colnames(data)[-1]<-colnames(indep)
          data<-as.data.frame(data)
          fit<-stats::lm(str2lang(paste(paste("`",colnames(data)[1],"`",sep=""),"~",
                                 gsub(", ","+",
                                      toString(
                                        paste('`',
                                              colnames(data)[-1],'`',sep=""))))),data)

          fits[[i]]<-fit

        }
      }else{I<-I+1}
    }
    if (ERROR==TRUE){
      warning(
        "The NSGA has failed, the hyper-parameters are restored to the initial values"
      )
      hyperparams<-tmp_hyper
    }
  }

  if (ERROR==TRUE){ # If not optimized, or cannot be optimized.


    hyperparams[is.na(hyperparams)]<-0

    if ("in" %in% latents){
      weight.X<-hyperparams[1:ncol(X)]
      if (rel_weight){
        weight.X<-weight.X*ncol(X)/sum(weight.X)
      }
      params.X<-hyperparams[-c(1:ncol(X))]
    }else{
      if ("out" %in% latents){
        weight.Y<-hyperparams[1:ncol(Y)]
        if (rel_weight){
          weight.Y<-weight.Y*ncol(Y)/sum(weight.Y)
        }
        params.Y<-hyperparams[-c(1:ncol(Y))]
      }else{
        if ("both" %in% latents){
          weight.X<-hyperparams[1:ncol(X)]
          params.X<-hyperparams[(ncol(X)+1):(ncol(X)+4)]
          weight.Y<-hyperparams[(ncol(X)+5):(ncol(X)+4+ncol(Y))]
          params.Y<-hyperparams[(ncol(X)+5+ncol(Y)):(ncol(X)+4+ncol(Y)+4)]
          if (rel_weight){
            weight.X<-weight.X*ncol(X)/sum(weight.X)
            weight.Y<-weight.Y*ncol(Y)/sum(weight.Y)
          }
        }
      }
    }
    if (latents %in% c("in","both")){ # For latent-independent variables
      NDA_in<-try(ndr(X,min_evalue = params.X[1],
                      min_communality=params.X[2],
                      com_communalities = params.X[3],
                      min_R = params.X[4],weight=weight.X,covar=FALSE,
                      cor_method=cor_method,
                      cor_type=cor_type,
                      min_comm=min_comm,
                      Gamma=Gamma,
                      null_model_type=null_model_type,
                      mod_mode=mod_mode,
                      use_rotation=use_rotation,
                      rotation=rotation),silent=TRUE)
    }

    if (latents %in% c("out","both")){ # For latent-dependent variables
      NDA_out<-try(ndr(Y,min_evalue = params.Y[1],
                       min_communality=params.Y[2],
                       com_communalities = params.Y[3],
                       min_R = params.Y[4],weight=weight.Y,covar=FALSE,
                       cor_method=cor_method,
                       cor_type=cor_type,
                       min_comm=min_comm,
                       Gamma=Gamma,
                       null_model_type=null_model_type,
                       mod_mode=mod_mode,
                       use_rotation=use_rotation,
                       rotation=rotation),silent=TRUE)
    }
    fits<-list()
    extra_vars.X<-FALSE
    extra_vars.Y<-FALSE
    dropped_X<-NULL
    dropped_Y<-NULL
    if (latents %in% c("in")){
      if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
        extra_vars.X<-TRUE
        dropped_X<-X[,NDA_in$membership==0]
      }
    }else{
      if (latents %in% c("out")){
        if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
          extra_vars.Y<-TRUE
          dropped_Y<-Y[,NDA_out$membership==0]
        }
      }else{
        if (latents %in% c("both")){
          if ((dircon==TRUE)&&(sum(NDA_in$membership==0)>0)){
            extra_vars.X<-TRUE
            dropped_X<-X[,NDA_in$membership==0]
          }
          if ((dircon==TRUE)&&(sum(NDA_out$membership==0)>0)){
            extra_vars.Y<-TRUE
            dropped_Y<-Y[,NDA_out$membership==0]
          }
        }
      }
    }

    dep<-Y
    if (latents %in% c("out","both")){
      if ((extra_vars.Y==TRUE)&&(!is.null(dropped_Y))){
        dep<-cbind(as.data.frame(NDA_out$scores),dropped_Y)
        dep<-as.data.frame(dep)
        colnames(dep)<-c(paste("NDAout",1:NDA_out$factors,sep=""),
                         colnames(Y)[NDA_out$membership==0])
      }else{
        dep<-as.data.frame(NDA_out$scores)
        colnames(dep)<-paste("NDAout",1:NDA_out$factors,sep="")
      }
    }
    indep<-X
    if (latents %in% c("in","both")){
      if ((extra_vars.X==TRUE)&&(!is.null(dropped_X))){
        indep<-cbind(as.data.frame(NDA_in$scores),dropped_X)
        indep<-as.data.frame(indep)
        colnames(indep)<-c(paste("NDAin",1:NDA_in$factors,sep=""),
                           colnames(X)[NDA_in$membership==0])
      }else{
        indep<-as.data.frame(NDA_in$scores)
        colnames(indep)<-paste("NDAin",1:NDA_in$factors,sep="")
      }
    }

    for (i in 1:ncol(dep))
    {
      data<-cbind(dep[,i],indep)
      colnames(data)[1]<-colnames(dep)[i]
      colnames(data)[-1]<-colnames(indep)
      data<-as.data.frame(data)
      fit<-stats::lm(str2lang(paste(paste("`",colnames(data)[1],"`",sep=""),"~",
                                    gsub(", ","+",
                                         toString(
                                           paste('`',
                                                 colnames(data)[-1],'`',sep=""))))),data)

      fits[[i]]<-fit
    }

  }

  P<-list()
  P$Call<-cl
  P$target<-target
  P$fval<-cost(hyperparams)
  P$hyperparams<-hyperparams
  P$pareto<-pareto
  P$X<-X
  P$Y<-Y
  P$latents<-latents
  if (latents %in% c("in","both")){
    P$NDAin<-NDA_in
    P$NDAin_weight<-weight.X
    P$NDAin_min_evalue<-params.X[1]
    P$NDAin_min_communality<-params.X[2]
    P$NDAin_com_communalities<-params.X[3]
    P$NDAin_min_R <- params.X[4]
  }
  if (latents %in% c("out","both")){
    P$NDAout<-NDA_out
    P$NDAout_weight<-weight.Y
    P$NDAout_min_evalue<-params.Y[1]
    P$NDAout_min_communality<-params.Y[2]
    P$NDAout_com_communalities<-params.Y[3]
    P$NDAout_min_R <- params.Y[4]
  }
  P$fits<-fits
  P$optimized<-optimize
  if (optimize==TRUE){
    P$NSGA<-NSGA
  }
  P$extra_vars.X<-extra_vars.X
  P$extra_vars.Y<-extra_vars.Y
  if (latents %in% c("in","both")){
    if (extra_vars.X==TRUE){
      P$dircon_X<-colnames(dropped_X)
    }
  }
  if (latents %in% c("out","both")){
    if (extra_vars.Y==TRUE){
      P$dircon_Y<-colnames(dropped_Y)
    }
  }
  P$seed<-seed
  P$fn<-"NDRLM"
  class(P)<-c("ndrlm","list")
  return(P)
}

