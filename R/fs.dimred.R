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

fs.dimred<-function(fn,data,min_comm=0.25){
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  s<-deparse(fn$Call)
  p<-NULL
  v<-as.character(fn$Call)
  if (length(v)>1){
    s<-gsub(v[2],"data",s) # replace dataset name to "data"
    s<-paste("psych::",s,sep = "") # works wit psych functions
  }else{
    stop(
      "Callback must be at least two elements!",
      call. = FALSE
    )
  }
  if ("principal" %in% as.character(fn$Call)) {
    p<-eval(str2lang(s))
  }else{
    if ("fa" %in% as.character(fn$Call)) {
      p<-eval(str2lang(s))
    }
  }
  return(p)
}
