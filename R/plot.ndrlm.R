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
### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND REGRESSION (NDRLM) ####
#' @export
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
