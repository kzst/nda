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
## SUMMARY FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
#' @export
summary.ndrlm <- function(object,  digits =  getOption("digits"), ...) {
  if (methods::is(object,"ndrlm")){
    Call<-object$Call
    target<-object$target
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
    results<-list(Call=Call,
                    fval=fval,
                    target=target,
                    pareto=pareto,
                    X = X,
                    Y = Y,
                    latents = latents,
                    NDAin=unlist(ifelse(latents %in% c("in","both"),
                                     list(NDAin),
                                     list(NULL))),
                    NDAin_weight=unlist(ifelse(latents %in% c("in","both"),
                                      list(NDAin_weight),
                                      list(NULL))),
                    NDAin_min_evalue=unlist(ifelse(latents %in% c("in","both"),
                                             list(NDAin_min_evalue),
                                             list(NULL))),
                    NDAin_min_communality=unlist(ifelse(latents %in% c("in","both"),
                                                 list(NDAin_min_communality),
                                                 list(NULL))),
                    NDAin_com_communalities=unlist(ifelse(latents %in% c("in","both"),
                                                      list(NDAin_com_communalities),
                                                      list(NULL))),
                    NDAin_min_R=unlist(ifelse(latents %in% c("in","both"),
                                                      list(NDAin_min_R),
                                                      list(NULL))),
                    NDAout=unlist(ifelse(latents %in% c("out","both"),
                                     list(NDAout),
                                     list(NULL))),
                    NDAout_weight=unlist(ifelse(latents %in% c("out","both"),
                                      list(NDAout_weight),
                                      list(NULL))),
                    NDAout_min_evalue=unlist(ifelse(latents %in% c("out","both"),
                                             list(NDAout_min_evalue),
                                             list(NULL))),
                    NDAout_min_communality=unlist(ifelse(latents %in% c("out","both"),
                                                 list(NDAout_min_communality),
                                                 list(NULL))),
                    NDAout_com_communalities=unlist(ifelse(latents %in% c("out","both"),
                                                      list(NDAout_com_communalities),
                                                      list(NULL))),
                    NDAout_min_R=unlist(ifelse(latents %in% c("out","both"),
                                                      list(NDAout_min_R),
                                                      list(NULL))),
                    fits = fits,
                    optimized=optimized,
                    NSGA=unlist(ifelse(optimized==TRUE,
                                       list(NSGA),
                                       list(NULL))),
                    extra_vars.X=extra_vars.X,
                    extra_vars.Y=extra_vars.Y,
                    dircon_X=unlist(ifelse((extra_vars.X==TRUE)&&latents %in% c("in","both"),
                                     list(dircon_X),
                                     list(NULL))),
                    dircon_Y=unlist(ifelse((extra_vars.Y==TRUE)&&latents %in% c("out","both"),
                                         list(dircon_Y),
                                         list(NULL))),
                    fn=fn)
    print.ndrlm(object)
  }
}
