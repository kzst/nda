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
biplot.nda <- function(x, main=NULL,...){
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop(
      "Package \"graphics\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if ("nda" %in% class(x)){
    if (is.null(x$scores)){
      stop("Biplot requires component scores. You need to run ndr from the raw data",
           call. = FALSE)
    }else{
      oldpar<-graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(oldpar))
      graphics::par(mfrow=c(x$factors,x$factors))
      op <- graphics::par(mar = rep(2.0,4))
      if(!is.null(main))
        op <- c(op, graphics::par(mar = graphics::par("mar")+c(0,0,1,0)))
      for (i in c(1:x$factors)){
        for (j in c(1:x$factors)){
          if (i==j){
            graphics::hist(x$scores[,i],col="cyan",prob=TRUE,
                           main = paste("NDA",i,sep=""),xlab="",ylab="")
            graphics::lines(stats::density(x$scores[,i]),col="red",lwd=2)
          }else{
            stats::biplot(x$scores[,c(i,j)],x$loadings[,c(i,j)],xlab="",ylab="")
          }
        }
      }
      if(!is.null(main))
        graphics::mtext(main, line = -1.2, outer = TRUE)
    }
  }else{
    stats::biplot(x,main,...)
  }
}
