#-----------------------------------------------------------------------------#
#                                                                             #
#  NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)                  #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona       #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: May 2022                                                     #
#-----------------------------------------------------------------------------#

#' @export

biplot.nda <- function(x, main=NULL,...){
  oldw <- getOption("warn")
  options(warn = -1)
  if (class(x)=="nda"){
    par(mfrow=c(x$factors,x$factors))
    op <- par(pty = "s")
    if(!is.null(main))
      op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    for (i in c(1:x$factors)){
      for (j in c(1:x$factors)){
        if (i==j){
          hist(x$scores[,i],col="cyan",prob=TRUE,main = paste("NDA",i,sep=""),xlab="",ylab="")
          lines(density(x$scores[,i]),col="red",lwd=2)
        }else{
          biplot(x$scores[,c(i,j)],x$loadings[,c(i,j)],xlab="",ylab="")
        }
      }
    }
    if(!is.null(main))
      mtext(main, line = -1.2, outer = TRUE)
  }else{
    biplot(x,main,...)
  }
  options(warn = oldw)
}


