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

###### BIPLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ###

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
  if (methods::is(x,"nda")){
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

#DATA GENERATION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#

data_gen<-function(n,m,nfactors=2,lambda=1){
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  M<-NA
  if (n>=1)
  {
    if (m>=1)
    {
      M<-matrix(0,nrow=n,ncol=m)
      if (nfactors>=1)
      {
        L<-replicate(nfactors,matrix(1,ceiling(n/nfactors),
                                     ceiling(m/nfactors)),simplify=FALSE)
        M<-Matrix::bdiag(L)
        M<-as.matrix(M[1:n,1:m])
        N<-matrix(stats::runif(n*m),n,m)
        M<-M-N*M/exp(lambda)
      }
      else
      {
        warning("nfactors must be equal to or greater than 1!")
      }
    }
    else
    {
      warning("m must be equal to or greater than 1!")
    }
  }
  else
  {
    warning("n must be equal to or greater than 1!")
  }
  return(as.data.frame(M))
}

######## MATRIX-BASED DISTANCE CORRELATION ########

dCor<-function(x,y=NULL){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(x)|is.matrix(x)){
    if (min(dim(x))>1){
      if (!is.null(y)){
        warning("x is a matrix or a data.frame with at least two columns, y is neglected.")
      }
      y<-NULL
    }
  }
  if (is.null(y)){
    if (is.data.frame(x)|is.matrix(x)){
      dC<-matrix(0,nrow=ncol(x),ncol=ncol(x))
      for (i in c(1:ncol(x))){
        for (j in c(1:ncol(x))){
          dC[i,j]<-energy::dcor(x[,i],x[,j])
        }
      }
      rownames(dC)<-colnames(x)
      colnames(dC)<-colnames(x)
      dCor<-dC
      dCor
    }else{
      dCor<-NULL
      stop("Error: x must be a matrix or a dataframe!")
    }
  }else{
    dCor<-energy::dcor(x,y)
    dCor
  }
}


######## MATRIX-BASED DISTANCE COVARIANCE ########

dCov<-function(x,y=NULL){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.null(y)){
    if (is.data.frame(x)|is.matrix(x)){
      dC<-matrix(0,nrow=ncol(x),ncol=ncol(x))
      for (i in c(1:ncol(x))){
        for (j in c(1:ncol(x))){
          dC[i,j]<-energy::dcov(x[,i],x[,j])
        }
      }
      rownames(dC)<-colnames(x)
      colnames(dC)<-colnames(x)
      dCov<-dC
      dCov
    }else{
      stop("Error: x must be a matrix or a dataframe!")
      dCov<-NULL
    }
  }else{
    dCov<-energy::dcov(x,y)
    dCov
  }
}

######## MATRIX-BASED DISTANCE PARTIAL CORRELATION ########

pdCor<-function(x){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package \"MASS\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))

  # sample number
  n <- dim(x)[1]

  # given variables' number
  gp <- dim(x)[2]-2

  # covariance matrix
  cvx <- dCov(x)

  # inverse covariance matrix
  if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
    icvx <- MASS::ginv(cvx)
  }else
    icvx <- Rfast::spdinv(cvx)

  rownames(icvx)<-rownames(cvx)
  colnames(icvx)<-colnames(cvx)
  # partial correlation
  pcor <- -stats::cov2cor(icvx)
  diag(pcor) <- 1
  pdCor<-pcor
  pdCor
}



######## MATRIX-BASED DISTANCE SEMI-PARTIAL CORRELATION ########

spdCor<-function(x){
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package \"MASS\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x))
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))

  # sample number
  n <- dim(x)[1]

  # given variables' number
  gp <- dim(x)[2]-2

  # covariance matrix
  cvx <- dCov(x)

  # inverse covariance matrix
  if(det(cvx) < .Machine$double.eps){
    warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
    icvx <- MASS::ginv(cvx)
  }else
    icvx <- Rfast::spdinv(cvx)

  rownames(icvx)<-rownames(cvx)
  colnames(icvx)<-colnames(cvx)

  # semi-partial correlation
  spcor <- -stats::cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
  diag(spcor) <- 1
  spdCor<-spcor
  spdCor
}


#### GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (GNDA) ####

ndr<-function(r,covar=FALSE,cor_method=1,cor_type=1,min_R=0,min_comm=2,Gamma=1,
              null_model_type=4,mod_mode=6,min_evalue=0,
              min_communality=0,com_communalities=0,use_rotation=FALSE,
              rotation="oblimin",weight=NULL,seed=NULL){

  cl<-match.call()
  if (!requireNamespace("energy", quietly = TRUE)) {
    stop(
      "Package \"energy\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ppcor", quietly = TRUE)) {
    stop(
      "Package \"ppcor\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("leidenAlg", quietly = TRUE)) {
    stop(
      "Package \"leidenAlg\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!is.numeric(as.matrix(r))) {
    stop(
      "The data should be numeric matrix or data.frame!",
      call. = FALSE
    )
  }
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  if (is.null(weight)){
    weight=rep(1,ncol(r))
  }
  r<-t(t(r)*weight)
  weight[is.na(weight)]<-0
  if (is.na(min_R)) {min_R<-0}
  if (is.na(min_evalue)) {min_evalue<-0}
  if (is.na(min_communality)) {min_communality<-0}
  if (is.na(com_communalities)) {com_communalities<-0}
  DATA<-r
  X<-r

  # Prepare correlation matrix

  if (covar==FALSE){
    if (cor_type==1){ # Bivariate correlations
      COR=switch(
        cor_method,
        "1"=stats::cor(X),
        "2"=stats::cor(X,method="spearman"),
        "3"=stats::cor(X,method="kendall"),
        "4"=dCor(X)
      )
    }else{
      if (cor_type==2){ # Partial correlations
        COR=switch(
          cor_method,
          "1"=ppcor::pcor(X)$estimate,
          "2"=ppcor::pcor(X,method="spearman")$estimate,
          "3"=ppcor::pcor(X,method="kendall")$estimate,
          "4"=pdCor(X)
        )
      }else{ # Semi-partial correlations
        COR=switch(
          cor_method,
          "1"=ppcor::spcor(X)$estimate,
          "2"=ppcor::spcor(X,method="spearman")$estimate,
          "3"=ppcor::spcor(X,method="kendall")$estimate,
          "4"=spdCor(X)
        )
      }
    }
  }else{
    COR<-X
  }
  COR[is.na(COR)]<-0
  issymm<-isSymmetric(as.matrix(COR))

  R<-COR^2
  R<-as.data.frame(R)
  colnames(R)<-colnames(r)
  rownames(R)<-colnames(r)
  remove(COR)

  R<-R-diag(nrow(R))

  R[R<min_R]<-0

  ## Calculate null modell

  kin<-colSums(R)
  kout<-rowSums(R)
  l=sum(R)
  N<-(kout %*% t(kin))/l

  # Calculate modularity

  coords<-matrix(1,nrow(R),1)

  Gamma<-1
  null_modell_type<-4

  MTX=switch(
    null_modell_type,
    "1"=R-N*Gamma,
    "2"=R-matrix(mean(R[R>0])*Gamma,nrow(R),ncol(R)),
    "3"=R-matrix(min_R*Gamma,nrow(R),ncol(R)),
    "4"=R
  )
  MTX[MTX<0]<-0
  cor_method<-1 # Non-linear correlation only used for the correlation graph
  if (issymm==TRUE) {
    modular=switch(
      mod_mode,
      "1"=igraph::cluster_louvain(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                      mode = "undirected", weighted = TRUE, diag = FALSE)),
      "2"=igraph::cluster_fast_greedy(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                          mode = "undirected", weighted = TRUE, diag = FALSE)),
      "3"=igraph::cluster_leading_eigen(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                            mode = "undirected", weighted = TRUE, diag = FALSE)),
      "4"=igraph::cluster_infomap(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                      mode = "undirected", weighted = TRUE, diag = FALSE)),
      "5"=igraph::cluster_walktrap(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                       mode = "undirected", weighted = TRUE, diag = FALSE)),
      "6"=if (inherits(try(leidenAlg::leiden.community(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "undirected", weighted = TRUE, diag = FALSE)),silent = TRUE),"try-error"))
      {igraph::cluster_leiden(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "undirected", weighted = TRUE, diag = FALSE),objective_function = "modularity")}
      else{leidenAlg::leiden.community(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "undirected", weighted = TRUE, diag = FALSE))}
    )
  }else{
    modular=switch(
      mod_mode,
      "1"=igraph::cluster_louvain(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                      mode = "max", weighted = TRUE, diag = FALSE)),
      "2"=igraph::cluster_fast_greedy(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                          mode = "max", weighted = TRUE, diag = FALSE)),
      "3"=igraph::cluster_leading_eigen(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                            mode = "max", weighted = TRUE, diag = FALSE)),
      "4"=igraph::cluster_infomap(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                      mode = "directed", weighted = TRUE, diag = FALSE)),
      "5"=igraph::cluster_walktrap(igraph::graph_from_adjacency_matrix(as.matrix(MTX),
                                                                       mode = "directed", weighted = TRUE, diag = FALSE)),
      "6"=if (inherits(try(leidenAlg::leiden.community(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "directed", weighted = TRUE, diag = FALSE)),silent = TRUE),"try-error"))
      {igraph::cluster_leiden(igraph::as.undirected(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "directed", weighted = TRUE, diag = FALSE)),objective_function = "modularity")}
      else{leidenAlg::leiden.community(igraph::graph_from_adjacency_matrix(as.matrix(MTX),mode = "directed", weighted = TRUE, diag = FALSE))}
    )
  }

  S<-as.numeric(modular$membership)

  for (i in 1: max(S)){
    if (nrow(as.matrix(coords[S==i]))<min_comm){
      coords[S==i]<-0
    }
  }

  S[coords==0]<-0

  # Estimate latent variables

  M<-sort(unique(S))
  if (min(M)>0)
  {
    M2=min(M):(length(M))
  }else{
    M2=min(M):(length(M)-1)
  }
  S<-M2[match(S,M)]
  M<-M2
  if (M[1]==0){
    M<-M[-1]
  }

  if (covar==FALSE){
    r<-X;
    is.na(r)<-sapply(r, is.infinite)
    r[is.na(r)]<-0
  }
  # Feature selection (1) - Drop peripheric items
  Coords<-c(1:nrow(as.matrix(S)))
  L<-matrix(0,nrow(DATA),nrow(as.matrix(M))) # Factor scores

  EVCs<-list()
  DATAs<-list()
  for (i in 1:nrow(as.matrix(M))){
    Coordsi<-Coords[(S==M[i])&(coords==1)]
    if (issymm==TRUE) {
      EVC<-as.matrix(igraph::eigen_centrality(igraph::graph_from_adjacency_matrix(
        as.matrix(R[Coordsi,Coordsi]), mode = "undirected",
        weighted = TRUE, diag = FALSE))$vector)
    }else{
      EVC<-as.matrix(igraph::eigen_centrality(igraph::graph_from_adjacency_matrix(
        as.matrix(R[Coordsi,Coordsi]), mode = "directed",
        weighted = TRUE, diag = FALSE))$vector)
    }
    if ((nrow(as.matrix(EVC[EVC>min_evalue]))>2)&(nrow(EVC)>2)){
      L[,i]<-if (inherits(try(as.matrix(rowSums(r[,
                                                  Coordsi[EVC>min_evalue]] * EVC[EVC>min_evalue])),
                              silent = TRUE),"try-error")) {as.matrix(rowSums(r[,
                                                                                Coordsi[EVC>min_evalue]] %*% EVC[EVC>min_evalue]))}
      else{as.matrix(rowSums(r[,Coordsi[EVC>min_evalue]] * EVC[EVC>min_evalue]))}
      coords[Coordsi[EVC<=min_evalue]]<-0
      coords[Coordsi[EVC<=min_evalue]]<-0
      S[Coordsi[EVC<=min_evalue]]<-0
    }else{
      L[,i]<-if (inherits(try(as.matrix(rowSums(r[,Coordsi] * EVC)),silent = TRUE),"try-error"))
      {as.matrix(rowSums(r[,Coordsi] %*% EVC))}else{as.matrix(rowSums(r[,Coordsi] * EVC))}
    }
    EVCs[[i]]=EVC[EVC>min_evalue]
    DATAs[[i]]=r[,S==M[i]];
  }
  if (ncol(L)>1 && use_rotation==TRUE){
    L<-psych::principal(L,nfactors = dim(L)[2],
                        rotate = rotation)$scores
  }else{
    L<-scale(L)
  }

  C=switch(
    cor_method,
    "1"=stats::cor(L),
    "2"=stats::cor(L,method="spearman"),
    "3"=stats::cor(L,method="kendall"),
    "4"=dCor(L)
  )
  CoordsS<-Coords[S!=0]
  CoordsC<-c(1:nrow(as.matrix(CoordsS)))
  if (covar==FALSE){
    LOADING=switch(
      cor_method,
      "1"=stats::cor(r[,S>0],L),
      "2"=stats::cor(r[,S>0],L,method="spearman"),
      "3"=stats::cor(r[,S>0],L,method="kendall"),
      "4"=dCor(r[,S>0],L)
    )
  }else{
    LOADING<-matrix(0,length(S),nrow(as.matrix(M))) # Factor scores
    for (i in 1:nrow(as.matrix(M))){
      LOADING[Coords[S==i],i]<-EVCs[[i]]
    }
    LOADING<-as.matrix(LOADING[Coords[S!=0],])
    rownames(LOADING)<-names(as.data.frame(r))[S>0]
  }
  COMMUNALITY<-t(apply(LOADING^2,1,max))

  # Feature selection (2) - Drop items with low communality

  COMMUNALITY<-t(apply(LOADING^2,1,max))
  COMMUNALITY[is.na(COMMUNALITY)]<-0
  max_it<-100
  it<-1
  while ((min(COMMUNALITY)<min_communality)&&(it<max_it)){
    it<-it+1
    COMMUNALITY<-t(apply(LOADING^2,1,max))
    COMMUNALITY[is.na(COMMUNALITY)]<-0
    CoordsS<-Coords[S!=0]
    CoordsC<-c(1:nrow(as.matrix(CoordsS)))
    s<-S[S!=0]
    coordsS<-coords[S!=0]
    for (i in 1:nrow(as.matrix(M))){
      Coordsi<-Coords[(S==M[i])&(coords==1)]
      CoordsiC<-CoordsC[(s==M[i])&(coordsS==1)]
      COM<-COMMUNALITY[CoordsiC]
      com_min<-min(COM)
      if (sum(COM>min_communality)>=2){
        S[Coordsi[COM<=min_communality]]<-0
        coords[Coordsi[COM<=min_communality]]<-0
        EVC<-EVCs[[i]]
        EVC<-EVC[COM>min_communality]
        EVCs[[i]]<-EVC

        L[,i]<-if (inherits(try(as.matrix(rowSums(r[,Coordsi[COM>min_communality]] * EVC)),silent = TRUE),"try-error"))
        {as.matrix(rowSums(r[,Coordsi[COM>min_communality]] %*% EVC))}else{
          as.matrix(rowSums(r[,Coordsi[COM>min_communality]] * EVC))
        }
      }else{
        EVC<-EVCs[[i]]
        L[,i]<-if (inherits(try(as.matrix(rowSums(r[,Coordsi] * EVC)),silent = TRUE),"try-error"))
        {as.matrix(rowSums(r[,Coordsi] %*% EVC))}else{as.matrix(rowSums(r[,Coordsi] * EVC))}
      }
    }
    if (ncol(L)>1 && use_rotation==TRUE){
      L<-psych::principal(L,nfactors = dim(L)[2],
                          rotate = rotation)$scores
    }else{
      L<-scale(L)
    }
    C=switch(
      cor_method,
      "1"=stats::cor(L),
      "2"=stats::cor(L,method="spearman"),
      "3"=stats::cor(L,method="kendall"),
      "4"=dCor(L)
    )
    if (covar==FALSE){
      LOADING=switch(
        cor_method,
        "1"=stats::cor(r[,S>0],L),
        "2"=stats::cor(r[,S>0],L,method="spearman"),
        "3"=stats::cor(r[,S>0],L,method="kendall"),
        "4"=dCor(r[,S>0],L)
      )
    }else{
      LOADING<-matrix(0,length(S),nrow(as.matrix(M))) # Factor scores
      for (i in 1:nrow(as.matrix(M))){
        LOADING[Coords[S==i],i]<-EVCs[[i]]
      }
      LOADING<-as.matrix(LOADING[Coords[S!=0],])
      rownames(LOADING)<-names(as.data.frame(r))[S>0]
    }
    COMMUNALITY<-t(apply(LOADING^2,1,max))
  }

  # Feature selection (3) - Drop items with high common communalities

  l<-FALSE
  while(l==FALSE){
    l<-TRUE
    CCs<-matrix(0,nrow(as.matrix(LOADING)),1)
    if (ncol(LOADING)>1){
      CoordsC=Coords[S!=0]
      L2<-LOADING^2
      nL2<-nrow(L2)
      for (I in 1:nL2){
        CJ<-max(L2[I,])
        CJ2<-max(L2[I,L2[I,]!=CJ]) #2nd maximal value;
        if ((CJ>=CJ2+com_communalities)||(CJ>2*CJ2)){

        }else{
          CCs[I]<-1
        }
      }
    }
    if (sum(CCs)>0){
      Coords_real<-CoordsC[CCs==1]
      COM<-COMMUNALITY[CCs==1]
      com<-sort(COM,index.return=TRUE)
      O_COM<-com[[1]]
      P_COM<-com[[2]]
      remove(com)
      Coords_real=Coords_real[P_COM]
      l<-TRUE
      i<-1
      if (nrow(as.matrix(S[S==S[Coords_real[i]]]))>2){
        l<-FALSE
        S[Coords_real]<-0
        coords[Coords_real]<-0
      }
      i<-i+1
    }
    for (i in 1:nrow(as.matrix(M))){
      Coordsi=Coords[(S==M[i])&(coords==1)]
      if (issymm==TRUE) {
        EVC<-as.matrix(igraph::eigen_centrality(igraph::graph_from_adjacency_matrix(
          as.matrix(R[Coordsi,Coordsi]), mode = "undirected",
          weighted = TRUE, diag = FALSE))$vector)
      }else{
        EVC<-as.matrix(igraph::eigen_centrality(igraph::graph_from_adjacency_matrix(
          as.matrix(R[Coordsi,Coordsi]), mode = "directed",
          weighted = TRUE, diag = FALSE))$vector)
      }
      EVCs[[i]]<-EVC
      result<-NA
      try(result <- as.matrix(rowSums(r[,Coordsi] %*% EVC)),silent=TRUE)
      if (is.null(nrow(is.nan(result)))){
        try(result <- as.matrix(rowSums(r[,Coordsi] * EVC)),silent=TRUE)
      }
      L[,i]<-result
    }
    centers<-colMeans(L)
    if (ncol(L)>1 && use_rotation==TRUE){
      L<-psych::principal(L,nfactors = dim(L)[2],
                          rotate = rotation)$scores
    }else{
      L<-scale(L)
    }
    C=switch(
      cor_method,
      "1"=stats::cor(L),
      "2"=stats::cor(L,method="spearman"),
      "3"=stats::cor(L,method="kendall"),
      "4"=dCor(L)
    )
    if (covar==FALSE){
      LOADING=switch(
        cor_method,
        "1"=stats::cor(r[,S>0],L),
        "2"=stats::cor(r[,S>0],L,method="spearman"),
        "3"=stats::cor(r[,S>0],L,method="kendall"),
        "4"=dCor(r[,S>0],L)
      )
    }else{
      LOADING<-matrix(0,length(S),nrow(as.matrix(M))) # Factor scores
      for (i in 1:nrow(as.matrix(M))){
        LOADING[Coords[S==i],i]<-EVCs[[i]]
      }
      LOADING<-as.matrix(LOADING[Coords[S!=0],])
      rownames(LOADING)<-names(as.data.frame(r))[S>0]
    }
    COMMUNALITY<-t(apply(LOADING^2,1,max))
  }

  P<-list()
  P$communality<-COMMUNALITY
  P$loadings<-LOADING
  colnames(P$loadings)<-paste("NDA",1:nrow(as.matrix(M)),sep = "")
  P$uniqueness<-1-COMMUNALITY
  P$factors<-nrow(as.matrix(M))
  if (covar==FALSE){
    P$scores<-L
    rownames(P$scores)<-rownames(DATA)
    colnames(P$scores)<-paste("NDA",1:nrow(as.matrix(M)),sep = "")
  }
  P$n.obs<-nrow(DATA)
  P$R<-R
  P$EVCs<-EVCs
  P$center<-centers
  P$membership<-S
  P$weight<-weight
  P$use_rotation<-use_rotation
  P$rotation<-rotation
  P$fn<-"NDA"
  P$seed<-seed
  P$Call<-cl
  class(P) <- c("nda","list")
  return(P)
}




###### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ######

plot.nda <- function(x,cuts=0.3,interactive=TRUE,edgescale=1.0,labeldist=-1.5,
                     show_weights=FALSE,...){
  if (methods::is(x,"nda")){
    R2<-G<-nodes<-edges<-NULL
    R2<-x$R
    R2[R2<cuts]<-0
    if (isSymmetric(as.matrix(R2))){
      G=igraph::graph.adjacency(as.matrix(R2), mode = "undirected",
                                weighted = TRUE, diag = FALSE)
    }else{
      G=igraph::graph.adjacency(as.matrix(R2), mode = "directed",
                                weighted = TRUE, diag = FALSE)
    }
    nodes<-as.data.frame(igraph::V(G)$name)
    nodes$label<-rownames(x$R)
    nodes$size<-igraph::evcent(G)$vector*10+5
    nodes$color<-grDevices::hsv(x$membership/max(x$membership),
                                alpha=0.4)
    nodes[x$membership==0,"color"]<-"#000000"
    colnames(nodes)<-c("id","title","size","color")
    edges<-as.data.frame(igraph::as_edgelist(G))
    edges <- data.frame(
      from=edges$V1,
      to=edges$V2,
      arrows=ifelse(igraph::is.directed(G),c("middle"),""),
      smooth=c(FALSE),
      label=unlist(ifelse(show_weights==TRUE,list(paste(round(igraph::E(G)$weight,2))),list(""))),
      width=(igraph::E(G)$weight)*edgescale,
      color="#5080b1"
    )

    nw <-
      visNetwork::visIgraphLayout(
        visNetwork::visNodes(
          visNetwork::visInteraction(
            visNetwork::visOptions(
              visNetwork::visEdges(
                visNetwork::visNetwork(
                  nodes, edges, height = "1000px", width = "100%"),
                font = list(size = 6)),
              highlightNearest = TRUE, selectedBy = "label"),
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE,
            hideEdgesOnDrag = FALSE),physics=FALSE, size=16,
          borderWidth = 1,
          font=list(face="calibri")),layout = "layout_nicely",
        physics = TRUE, type="full"
      )

    if (interactive==FALSE){
      g <- igraph::graph_from_data_frame(edges,vertices=nodes,
                                         directed = igraph::is.directed(G))
      igraph::E(g)$weight<-igraph::E(G)$weight
      igraph::E(g)$size<-igraph::E(G)$weight
      igraph::plot.igraph(g, vertex.label.dist = labeldist,vertex.size=nodes$size,edge.width=(igraph::E(g)$size*5+1)*edgescale,edge.arrow.size=0.2)
      return(invisible(g))
    }else{
      nw
    }
  }
}


#SUMMARY FUNCTION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#

summary.nda <- function(object,  digits =  getOption("digits"), ...) {
  if (methods::is(object,"nda")){
    communality <- object$communality
    loadings <- object$loadings
    uniqueness <- object$uniqueness
    factors <- object$factors
    scores <- object$scores
    n.obs <- object$n.obs
    factors <- object$factors
    if (!is.null(scores)){
      results<-list(cummunality = communality, loadings = loadings,
                    uniqueness = uniqueness,
                    factors = factors,
                    scores = scores,
                    n.obs = n.obs)
    }else{
      results<-list(cummunality = communality, loadings = loadings,
                    uniqueness = uniqueness,
                    factors = factors,
                    n.obs = n.obs)
    }
    return(results)
    print.nda(object)
  }
}


# PRINT FUNCTION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#

print.nda <- function(x,  digits =  getOption("digits"), ...) {
  if (methods::is(x,"nda")){
    communality <- x$communality
    loadings <- x$loadings
    uniqueness <- x$uniqueness
    factors <- x$factors
    scores <- x$scores
    n.obs <- x$n.obs
    factors <- x$factors
    cat("\nPrint of the NDA:\n")
    cat("\nNumber of latent variables: ",factors)
    cat("\nNumber of observations: ",n.obs)
    cat("\nCommunalities:\n")
    print(communality,digits = digits, ...)
    cat("\nFactor loadings:\n")
    print(loadings,digits = digits, ...)
    if (!is.null(scores)){
      cat("\nFactor scores:\n")
      print(scores,digits = digits, ...)
      cat("\n\nCorrelation matrix of factor scores:\n")
      print(stats::cor(scores),digits = digits, ...)
    }
  }
}


######### Feature selection for KMO #######

fs.KMO<-function(data,min_MSA=0.5,cor.mtx=FALSE){
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.data.frame(data)|is.matrix(data)){
    if (ncol(data)>=2){
      x<-data
      loop=TRUE
      while(loop==TRUE){
        kmo<-psych::KMO(x)
        if (min(kmo$MSAi)>min_MSA){loop=FALSE}else{
          i<-which.min(kmo$MSAi)
          if (cor.mtx==TRUE){
            x<-x[-i,-i]
          }else{
            x<-x[,-i]
          }
        }
        if (ncol(x)<2){
          loop=FALSE
        }
      }
      return(x)
    }else{
      stop("Error: data must contain at least 2 columns!")
      step.KMO<-NULL
      return(step.KMO)
    }
  }else{
    stop("Error: data must be a matrix or a dataframe!")
    step.KMO<-NULL
    return(step.KMO)
  }
}


######### Feature selection for PCA/FA/NDA #######

fs.dimred<-function(fn,DF,min_comm=0.25,com_comm=0.25){
  if (!requireNamespace("psych", quietly = TRUE)) {
    stop(
      "Package \"psych\" must be installed to use this function.",
      call. = FALSE
    )
  }
  DF<-as.data.frame(DF)
  s<-deparse1(fn$Call)
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
          break
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

######### Normalize entire data, row, or column #######

normalize <- function(x,type="all")
{
  results<-NULL
  if ((is.data.frame(x))|(is.matrix(x))|
      (is.array(x))){
    results<-((x - min(x)) / (max(x) - min(x)))
    if ("row" %in% type){
      for (i in 1:nrow(x)){
        results[i,]<-((x[i,] - min(x[i,])) / (max(x[i,]) - min(x[i,])))
      }
    }else{
      if ("col" %in% type){
        for (i in 1:ncol(x)){
          results[,i]<-((x[,i] - min(x[,i])) / (max(x[,i]) - min(x[,i])))
        }
      }
    }
  }else{
    stop(
      "Only matrix, array or data.frame can be used in this function!",
      call. = FALSE
    )
  }
  return(results)
}

### GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (GNDR) ##

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

## PRINT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
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

## SUMMARY FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##
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


### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ####
plot.ndrlm <- function(x,sig=0.05,interactive=FALSE,...){
  if (methods::is(x,"ndrlm")){
    latents<-x$latents
    extra_vars.X<-x$extra_vars.X
    extra_vars.Y<-x$extra_vars.Y

    X<-x$X
    Y<-x$Y

    nY<-ncol(x$Y)
    nSin<-0
    if (latents %in% c("in","both")){
      nSin<-ncol(x$NDAin$scores)
      membership.X<-x$NDAin$membership
      loadings.X<-x$NDAin$loadings
    }
    nSout<-0
    if (latents %in% c("out","both")){
      nSout<-ncol(x$NDAout$scores)
      membership.Y<-x$NDAout$membership
      loadings.Y<-x$NDAout$loadings
    }
    nX<-ncol(x$X)

    node_ID<-1:(nY+nSout+nSin+nX)
    node_label<-c(colnames(x$Y),
                  unlist(ifelse(latents %in% c("out","both"),
                                list(paste("NDAout",1:x$NDAout$factors,sep="")),
                                list(NULL))),
                  unlist(ifelse(latents %in% c("in","both"),
                                list(paste("NDAin",1:x$NDAin$factors,sep="")),
                                list(NULL))),colnames(x$X))

    node_shape<-c(rep("rectangle",nY),
                  unlist(ifelse(latents %in% c("out","both"),
                                list(rep("circle",nSout)),
                                list(NULL))),
                  unlist(ifelse(latents %in% c("in","both"),
                                list(rep("circle",nSin)),
                                list(NULL))),
                  rep("rectangle",nX))

    node_color<-c(unlist(ifelse(latents %in% c("out","both"),
                                list(x$NDAout$membership),
                                list(rep(0,nY)))),
                  unlist(ifelse(latents %in% c("out","both"),
                                list(1:nSout),
                                list(NULL))),
                  unlist(ifelse(latents %in% c("in","both"),
                                list(1:nSin),
                                list(NULL))),
                  unlist(ifelse(latents %in% c("in","both"),
                                list(x$NDAin$membership),
                                list(rep(0,nX)))))
    nodes<-data.frame(id=node_ID,label=node_label,shape=node_shape,
                      color=node_color)
    edges <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(edges) <- c('from', 'to', 'weight' , 'color' , 'lty' , 'dashes')

    dep<-Y
    if (latents %in% c("out","both")){
      if (extra_vars.Y==TRUE){
        dep<-cbind(x$NDAout$scores,x$Y[,x$NDAout$membership==0])
        dep<-as.data.frame(dep)
        colnames(dep)<-c(paste("NDAout",1:x$NDAout$factors,sep=""),
                         colnames(x$Y)[x$NDAout$membership==0])
      }else{
        dep<-x$NDAout$scores
        colnames(dep)<-paste("NDAout",1:x$NDAout$factors,sep="")
      }
    }
    indep<-X
    if (latents %in% c("in","both")){
      if (extra_vars.X==TRUE){
        indep<-cbind(x$NDAin$scores,x$X[,x$NDAin$membership==0])
        indep<-as.data.frame(indep)
        colnames(indep)<-c(paste("NDAin",1:x$NDAin$factors,sep=""),
                           colnames(x$X)[x$NDAin$membership==0])



      }else{
        indep<-x$NDAin$scores
        colnames(indep)<-paste("NDAin",1:x$NDAin$factors,sep="")
      }
    }

    k<-1
    for (i in 1:length(x$fits)){
      coefs<-as.vector(lm.beta::lm.beta(x$fits[[i]])$standardized.coefficients)[-1]
      pvalues<-summary(x$fits[[i]])$coefficients[-1,4]
      indepvars<-colnames(x$fits[[i]]$model)[-1]
      depvar<-colnames(x$fits[[i]]$model)[1]
      for (j in 1:length(coefs)){
        if (pvalues[j]<sig){
          edges[k,"to"]<-node_ID[node_label %in% depvar]
          edges[k,"from"]<-node_ID[node_label %in% indepvars[j]]
          edges[k,"weight"]<-coefs[j]
          edges[k,"color"]<-"black"
          edges[k,"lty"]<-"solid"
          edges[k,"dashes"]<-FALSE
          k<-k+1
        }
      }
    }

    if (latents %in% c("in","both")){
      membership.X<-x$NDAin$membership
      for (i in 1:nSin){
        for (j in 1:length(membership.X)){
          if (membership.X[j]==i){
            edges[k,"from"]<-node_ID[node_label %in% colnames(x$X)[j]]
            edges[k,"to"]<-node_ID[node_label %in% paste("NDAin",i,sep="")]
            edges[k,"weight"]<-loadings.X[colnames(x$X)[j],i]
            edges[k,"color"]<-"grey"
            edges[k,"lty"]<-"dashed"
            edges[k,"dashes"]<-TRUE
            k<-k+1
          }
        }
      }
    }

    if (latents %in% c("out","both")){
      membership.Y<-x$NDAout$membership
      for (i in 1:nSout){
        for (j in 1:length(membership.Y)){
          if (membership.Y[j]==i){
            edges[k,"from"]<-node_ID[node_label %in% colnames(x$Y)[j]]
            edges[k,"to"]<-node_ID[node_label %in% paste("NDAout",i,sep="")]
            edges[k,"weight"]<-loadings.Y[colnames(x$Y)[j],i]
            edges[k,"color"]<-"grey"
            edges[k,"lty"]<-"dashed"
            edges[k,"dashes"]<-TRUE
            k<-k+1
          }
        }
      }
    }


    space<-100
    cust_layout<-matrix(0,ncol=2,nrow=nY+nSin+nSout+nX)
    cust_layout[1:nY,1]<-30
    if (latents %in% c("out","both")){
      cust_layout[sort(membership.Y,index.return=TRUE)$ix,2]<-((1:nY)-mean(1:nY))*space
    }else{
      cust_layout[1:nY,2]<-((1:nY)-mean(1:nY))*space
    }

    if (latents %in% c("out","both")){
      cust_layout[(nY+1):(nY+nSout),1]<-20
      cust_layout[(nY+1):(nY+nSout),2]<-((1:nSout)-mean(1:nSout))*space
    }

    if (latents %in% c("in","both")){
      cust_layout[(nY+nSout+1):(nY+nSin+nSout),1]<-10
      cust_layout[(nY+nSout+1):(nY+nSin+nSout),2]<-((1:nSin)-mean(1:nSin))*space
    }

    cust_layout[(nY+nSin+nSout+1):(nY+nSin+nSout+nX),1]<-0
    if (latents %in% c("in","both")){
      cust_layout[sort(membership.X,index.return=TRUE)$ix+nY+nSin+nSout,2]<-((1:nX)-mean(1:nX))*space
    }else{
      cust_layout[(nY+nSin+nSout+1):(nY+nSin+nSout+nX),2]<-((1:nX)-mean(1:nX))*space
    }
    if (latents %in% c("in","both")){
      for (i in 1:nSin){
        cust_layout[(nY+nSout+i),2]<-mean(cust_layout[cust_layout[,1]==0,2]*(membership.X==i))
      }
    }
    if (latents %in% c("out","both")){
      for (i in 1:nSout){
        cust_layout[(nY+i),2]<-mean(cust_layout[cust_layout[,1]==30,2]*(membership.Y==i))
      }
    }

    G<-igraph::graph_from_data_frame(edges,
                                     directed=TRUE,
                                     vertices=nodes)

    if (interactive==TRUE){
      edges$arrows<-ifelse(igraph::is.directed(G),c("to"),"")
      edges$width<-(abs(igraph::E(G)$weight))
      nodes$color<-grDevices::hsv((node_color+1)/max(node_color+1),
                                  alpha=0.4)
      nodes$shape<-gsub("rectangle","box",nodes$shape)
      nodes$shape<-gsub("circle","ellipse",nodes$shape)
      edges$label<-as.vector(paste(round(edges$weight,2),sep=""))
      nw <-
        visNetwork::visIgraphLayout(
          visNetwork::visNodes(
            visNetwork::visInteraction(
              visNetwork::visOptions(
                visNetwork::visEdges(
                  visNetwork::visNetwork(
                    nodes, edges, height = "1000px", width = "100%"),
                  font = list(size = 6),color="#555555",
                  label=edges$label),
                highlightNearest = TRUE, selectedBy = "label"),
              dragNodes = TRUE,
              dragView = TRUE,
              zoomView = TRUE,
              hideEdgesOnDrag = FALSE),physics=FALSE, size=16,
            borderWidth = 1,
            shape=nodes$shape,
            font=list(face="calibri")),layout="layout.norm",
          layoutMatrix = cust_layout,
          physics = FALSE, type="full"
        )
      nw
      return(nw)

    }else{
      igraph::V(G)$color<-grDevices::hsv((node_color+1)/max(node_color+1),
                                         alpha=0.4)
      igraph::plot.igraph(G,layout=cust_layout,edge.width=abs(igraph::E(G)$weight)*2,
                          edge.label=round(igraph::E(G)$weight,2),vertex.size=30)
      return(invisible(G))
    }
  }
}

# FITTINGS FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ##

fitted.ndrlm <- function(object,  ...) {
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
    FITTED<-as.data.frame(matrix(0,nrow=nrow(dep),ncol=ncol(dep)))
    colnames(FITTED)<-colnames(dep)
    rownames(FITTED)<-rownames(dep)
    for (i in 1:length(fits)){
      FITTED[,i]<-stats::fitted(fits[[i]])
    }
    FITTED<-FITTED[,1:length(fits)]
    return(FITTED)
  }
}

### PREDICT SCORES NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ##

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

### PREDICT SCORES NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ##

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

