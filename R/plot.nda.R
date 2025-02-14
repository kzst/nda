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
###### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) #####
#' @export
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
