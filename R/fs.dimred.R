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

fs.dimred<-function(fn,DF,min_comm=0.25,com_comm=0.25){
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  s<-deparse(fn$Call)
  p<-fn
  v<-as.character(fn$Call)
  if (length(v)<2){stop(
    "Callback must be at least two elements!",
    call. = FALSE
  )}
  s<-gsub(v[2],"DF",s,fixed=TRUE) # replace dataset name to "DF"
  if ("principal" %in% as.character(fn$Call)) {
    s<-paste("psych::",s,sep = "") # works with psych functions
  }else{
    if ("fa" %in% as.character(fn$Call)) {
      s<-paste("psych::",s,sep = "") # works with psych functions
    }else{
      if ("ndr" %in% as.character(fn$Call)) {
        s<-paste("nda::",s,sep = "") # works with nda functions
      }else{stop(
        "Feature selection only works with principal, fa, and ndr functions!",
        call. = FALSE
      )}
    }
  }
  dropped_low<-NULL
  loop=TRUE
  while(loop==TRUE){ # Drop low communality values
    p<-eval(str2lang(s))
    if (is.null(p$communality)==TRUE){loop=FALSE}else{
      if (min(p$communality)>=min_comm){loop=FALSE}else{
        i<-which.min(p$communality)
        if (is.null(p$scores)==TRUE){
          DF<-DF[-i,-i] # there is no score value => correlation matrix is
          #investigated
        }else{
          if (is.null(dropped_low)==TRUE){
            dropped_low<-eval(str2lang(paste("as.",class(DF[1]),"(DF[,i])",sep="")))
            names(dropped_low)[1]<-names(DF)[i]
          }else{
            dropped_low<-cbind(dropped_low,DF[i])
          }
          DF<-DF[,-i]
        }
      }
    }
    if (ncol(DF)<3){
      loop=FALSE
    }
  }
  dropped_com<-NULL
  repeat{
    p<-eval(str2lang(s))
    if (is.null(p$communality)==TRUE){
      break
    }else{
      if (is.null(p$loadings)==TRUE){
        break
        }else{
        if (ncol(p$loadings)<2){
          loop=FALSE
          }else{
          l<-abs(p$loadings)
          c<-matrix(0,ncol=1,nrow=nrow(l))
          for (i in 1:nrow(l)){
            r<-l[i,]
            m1<-max(r) # highest loading value
            m2<-max(r[-which.max(r)]) # 2nd highest loading value
            if ((m1<2*m2)&(m1<(m2+com_comm))){
              c[i]<-1
            }
          }
          if (sum(c)<1){
            break
          }
        }
        sel<-setdiff(as.vector(c*1:nrow(as.matrix(p$communality))),0)
        i<-sel[which.min(p$communality[sel])]
        if (is.null(p$scores)==TRUE){
          DF<-DF[-i,-i] # there is no score value => correlation matrix is
          #investigated
        }else{
          if (is.null(dropped_com)==TRUE){
            dropped_com<-eval(str2lang(paste("as.",class(DF)[1],"(DF[,i])",sep="")))
            names(dropped_com)[1]<-names(DF)[i]
          }else{
            dropped_com<-cbind(dropped_com,DF[i])
          }
          DF<-DF[,-i]
        }
      }
    }
    if (ncol(DF)<3){
      break
    }
  }
  p$dropped_low<-dropped_low
  p$dropped_com<-dropped_com
  p$retained_DF<-DF
  return(p)
}
