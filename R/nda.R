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

###### BIPLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) #####

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

# DATA GENERATION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) #

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
      stop("Error: x must be a matrix or a dataframe!")
      dCor<-NULL
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
    icvx <- solve(cvx)

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
    icvx <- solve(cvx)

  # semi-partial correlation
  spcor <- -stats::cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
  diag(spcor) <- 1
  spdCor<-spcor
  spdCor
}


########### NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ###########

ndr<-function(r,covar=FALSE,cor_method=1,cor_type=1,min_R=0,min_comm=2,Gamma=1,
              null_modell_type=4,mod_mode=6,min_evalue=0,
              min_communality=0,com_communalities=0,use_rotation=FALSE){

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
  if (!requireNamespace("leidenAlg", quietly = TRUE)) {
    stop(
      "Package \"leidenAlg\" must be installed to use this function.",
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
  if (issymm==FALSE){
    if (mod_mode<4){
      stop(
        "If correlation/simmilarity matrix is non-symmetric only InfoMap/Walktrap/Leiden modularities can be used.",
        call. = FALSE
      )
    }
  }
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
      "1"=igraph::cluster_louvain(igraph::graph.adjacency(as.matrix(MTX),
                                                          mode = "undirected", weighted = TRUE, diag = FALSE)),
      "2"=igraph::cluster_fast_greedy(igraph::graph.adjacency(as.matrix(MTX),
                                                              mode = "undirected", weighted = TRUE, diag = FALSE)),
      "3"=igraph::cluster_leading_eigen(igraph::graph.adjacency(as.matrix(MTX),
                                                                mode = "undirected", weighted = TRUE, diag = FALSE)),
      "4"=igraph::cluster_infomap(igraph::graph.adjacency(as.matrix(MTX),
                                                          mode = "undirected", weighted = TRUE, diag = FALSE)),
      "5"=igraph::cluster_walktrap(igraph::graph.adjacency(as.matrix(MTX),
                                                           mode = "undirected", weighted = TRUE, diag = FALSE)),
      "6"=leidenAlg::leiden.community(igraph::graph.adjacency(as.matrix(MTX),
                                                              mode = "undirected", weighted = TRUE, diag = FALSE))
    )
  }else{
    modular=switch(
      mod_mode,
      "1"=igraph::cluster_louvain(igraph::graph.adjacency(MTX,
                                                          mode = "directed", weighted = TRUE, diag = FALSE)),
      "2"=igraph::cluster_fast_greedy(igraph::graph.adjacency(as.matrix(MTX),
                                                              mode = "directed", weighted = TRUE, diag = FALSE)),
      "3"=igraph::cluster_leading_eigen(igraph::graph.adjacency(as.matrix(MTX),
                                                                mode = "directed", weighted = TRUE, diag = FALSE)),
      "4"=igraph::cluster_infomap(igraph::graph.adjacency(as.matrix(MTX),
                                                          mode = "directed", weighted = TRUE, diag = FALSE)),
      "5"=igraph::cluster_walktrap(igraph::graph.adjacency(as.matrix(MTX),
                                                           mode = "directed", weighted = TRUE, diag = FALSE)),
      "6"=leidenAlg::leiden.community(igraph::graph.adjacency(as.matrix(MTX),
                                                              mode = "directed", weighted = TRUE, diag = FALSE))
    )
  }

  S<-as.numeric(modular$membership)

  # igraph::sizes(modular)

  for (i in 1: max(S)){
    if (nrow(as.matrix(coords[S==i]))<min_comm){
      coords[S==i]<-0
    }
  }

  S[coords==0]<-0

  # Estimate latent variables

  M<-sort(unique(S))
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
      EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(
        as.matrix(R[Coordsi,Coordsi]), mode = "undirected",
        weighted = TRUE, diag = FALSE))$vector)
    }else{
      EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(
        as.matrix(R[Coordsi,Coordsi]), mode = "directed",
        weighted = TRUE, diag = FALSE))$vector)
    }
    if ((nrow(as.matrix(EVC[EVC>min_evalue]))>2)&(nrow(EVC)>2)){
      L[,i]<-as.matrix(rowSums(r[,
                                 Coordsi[EVC>min_evalue]] * EVC[EVC>min_evalue]))
      coords[Coordsi[EVC<=min_evalue]]<-0
      coords[Coordsi[EVC<=min_evalue]]<-0
      S[Coordsi[EVC<=min_evalue]]<-0
    }else{
      L[,i]<-as.matrix(rowSums(r[,Coordsi] * EVC))
    }
    EVCs[[i]]=EVC[EVC>min_evalue]
    DATAs[[i]]=r[,S==M[i]];
  }
  if (ncol(L)>1 && use_rotation==TRUE){
    L<-psych::principal(L,nfactors = dim(L)[2])$scores
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
    LOADING<-LOADING[-Coords[S==0],]
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
        L[,i]<-as.matrix(rowSums(r[,Coordsi[COM>min_communality]] * EVC))
      }else{
        EVC<-EVCs[[i]]
        L[,i]<-as.matrix(rowSums(r[,Coordsi] * EVC))
      }
    }
    if (ncol(L)>1 && use_rotation==TRUE){
      L<-psych::principal(L,nfactors = dim(L)[2])$scores
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
      LOADING<-LOADING[-Coords[S==0],]
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
        EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(
          as.matrix(R[Coordsi,Coordsi]), mode = "undirected",
          weighted = TRUE, diag = FALSE))$vector)
      }else{
        EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(
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
    if (ncol(L)>1 && use_rotation==TRUE){
      L<-psych::principal(L,nfactors = dim(L)[2])$scores
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
      LOADING<-LOADING[-Coords[S==0],]
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
  P$membership<-S
  P$fn<-"NDA"
  P$Call<-cl
  class(P) <- "nda"
  return(P)
}

####### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) ######

plot.nda <- function(x,cuts=0.3,...){
  if ("nda" %in% class(x)){
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
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
      stop(
        "Package \"visNetwork\" must be installed to use this function.",
        call. = FALSE
      )
    }
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
    nodes$color<-grDevices::hsv(x$membership/max(x$membership))
    nodes[x$membership==0,"color"]<-"#000000"
    colnames(nodes)<-c("id","title","color")
    edges<-as.data.frame(igraph::as_edgelist(G))
    edges <- data.frame(
      from=edges$V1,
      to=edges$V2,
      smooth=c(FALSE),
      width=igraph::E(G)$weight,
      color="#5080b1"
    )

    nw <-
      visNetwork::visIgraphLayout(
        visNetwork::visNodes(
          visNetwork::visInteraction(
            visNetwork::visOptions(
              visNetwork::visNetwork(
                nodes, edges, height = "1000px", width = "100%"),
              highlightNearest = TRUE, selectedBy = "label"),
            dragNodes = TRUE,
            dragView = TRUE,
            zoomView = TRUE,
            hideEdgesOnDrag = FALSE),physics=FALSE, size=16,
          borderWidth = 1,
          font=list(face="calibri")),layout = "layout_nicely",
        physics = TRUE, type="full"
      )
    nw
  }else{
    plot(x,...)
  }
}

# SUMMARY FUNCTION FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA)#

summary.nda <- function(object,  digits =  getOption("digits"), ...) {
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop(
      "Package \"stats\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if ("nda" %in% class(object)){
    communality <- object$communality
    loadings <- object$loadings
    uniqueness <- object$uniqueness
    factors <- object$factors
    scores <- object$scores
    n.obs <- object$n.obs
    factors <- object$factors
    cat("\nSummary of the NDA:\n")
    cat("\nNumber of latent variables: ",factors)
    cat("\nNumber of observations: ",n.obs)
    cat("\nCommunalities:\n")
    print(communality,digits = digits, ...)
    cat("\nFactor loadings:\n")
    print(loadings,digits = digits, ...)
    cat("\n\nCorrelaction matrix of factor scores:\n")
    print(stats::cor(scores),digits = digits, ...)
  }else{
    summary(object,...)
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
