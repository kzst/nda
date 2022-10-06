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

ndr<-function(data,cor_method=1,min_R=0,min_comm=2,Gamma=1,null_modell_type=4,mod_mode=6,min_evalue=0,min_communality=0,com_communalities=0,use_rotation=FALSE){


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
  DATA<-data
  X<-data

  # Prepare correlation matrix

  COR=switch(
    cor_method,
    "1"=stats::cor(X),
    "2"=stats::cor(X,method="spearman"),
    "3"=stats::cor(X,method="kendall"),
    "4"=dCor(X)
  )
  COR[is.na(COR)]<-0
  R<-COR^2
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
  modular=switch(
    mod_mode,
    "1"=igraph::cluster_louvain(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE)),
    "2"=igraph::cluster_fast_greedy(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE)),
    "3"=igraph::cluster_leading_eigen(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE)),
    "4"=igraph::cluster_infomap(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE)),
    "5"=igraph::cluster_walktrap(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE)),
    "6"=leidenAlg::leiden.community(igraph::graph.adjacency(MTX, mode = "undirected", weighted = TRUE, diag = FALSE))
  )

  S<-as.numeric(modular$membership)

  igraph::sizes(modular)

  for (i in 1: max(S)){
    if (nrow(as.matrix(coords[S==i]))<min_comm){
      coords[S==i]<-0
    }else{
      #image(rotateImage(R[S==i,S==i],angle=270), axes = FALSE, col = brewer.pal(9, "Blues"))
    }
  }

  S[coords==0]<-0

  # Estimate latent variables

  M<-sort(unique(S))
  if (M[1]==0){
    M<-M[-1]
  }

  data<-X;
  is.na(data)<-sapply(data, is.infinite)
  data[is.na(data)]<-0

  # Feature selection (1) - Drop peripheric items

  Coords<-c(1:nrow(as.matrix(S)))
  L<-matrix(0,nrow(DATA),nrow(as.matrix(M)))

  EVCs<-list()
  DATAs<-list()

  for (i in 1:nrow(as.matrix(M))){
    Coordsi<-Coords[(S==M[i])&(coords==1)]
    EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(R[Coordsi,Coordsi], mode = "undirected", weighted = TRUE, diag = FALSE))$vector)
    if ((nrow(as.matrix(EVC[EVC>min_evalue]))>2)&(nrow(EVC)>2)){
      L[,i]<-as.matrix(rowSums(data[,Coordsi[EVC>min_evalue]] * EVC[EVC>min_evalue]))
      coords[Coordsi[EVC<=min_evalue]]<-0
      coords[Coordsi[EVC<=min_evalue]]<-0
      S[Coordsi[EVC<=min_evalue]]<-0
    }else{
      L[,i]<-as.matrix(rowSums(data[,Coordsi] * EVC))
    }
    EVCs[[i]]=EVC[EVC>min_evalue]
    DATAs[[i]]=data[,S==M[i]];
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
  LOADING=switch(
    cor_method,
    "1"=stats::cor(data[,S>0],L),
    "2"=stats::cor(data[,S>0],L,method="spearman"),
    "3"=stats::cor(data[,S>0],L,method="kendall"),
    "4"=dCor(data[,S>0],L)
  )
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
        L[,i]<-as.matrix(rowSums(data[,Coordsi[COM>min_communality]] * EVC))
      }else{
        EVC<-EVCs[[i]]
        L[,i]<-as.matrix(rowSums(data[,Coordsi] * EVC))
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
    LOADING=switch(
      cor_method,
      "1"=stats::cor(data[,S>0],L),
      "2"=stats::cor(data[,S>0],L,method="spearman"),
      "3"=stats::cor(data[,S>0],L,method="kendall"),
      "4"=dCor(data[,S>0],L)
    )
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
      EVC<-as.matrix(igraph::eigen_centrality(igraph::graph.adjacency(R[Coordsi,Coordsi], mode = "undirected", weighted = TRUE, diag = FALSE))$vector)
      EVCs[[i]]<-EVC
      result<-NA
      try(result <- as.matrix(rowSums(data[,Coordsi] %*% EVC)),silent=TRUE)
      if (is.null(nrow(is.nan(result)))){
        try(result <- as.matrix(rowSums(data[,Coordsi] * EVC)),silent=TRUE)
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
    LOADING=switch(
      cor_method,
      "1"=stats::cor(data[,S>0],L),
      "2"=stats::cor(data[,S>0],L,method="spearman"),
      "3"=stats::cor(data[,S>0],L,method="kendall"),
      "4"=dCor(data[,S>0],L)
    )
    COMMUNALITY<-t(apply(LOADING^2,1,max))
  }

  P<-list()
  P$communality<-COMMUNALITY
  P$loadings<-LOADING
  colnames(P$loadings)<-paste("NDA",1:nrow(as.matrix(M)),sep = "")
  P$uniqueness<-1-COMMUNALITY
  P$factors<-nrow(as.matrix(M))
  P$scores<-L
  rownames(P$scores)<-rownames(DATA)
  colnames(P$scores)<-paste("NDA",1:nrow(as.matrix(M)),sep = "")
  P$n.obs<-nrow(DATA)
  P$R<-R
  P$membership<-S
  P$fn<-"NDA"
  class(P) <- "nda"
  return(P)
}

