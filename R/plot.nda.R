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
# Last modified: August 2026                                                #
#-----------------------------------------------------------------------------#


#' Plot an NDA association network
#' @export
plot.nda <- function(x, cuts = 0.3, interactive = TRUE, edgescale = 1,
                     labeldist = -1.5, show_weights = FALSE,
                     shapes = "circle", edge.color.links = TRUE,
                     edge.colors = c("darkred", "grey", "darkblue"),
                     edge.label.size = 10, edge.label.color = "#333333",
                     vertex.label.size = 12, vertex.label.color = "black",
                     layout = "nicely", igraph_args = list(),
                     visnetwork_args = list(), ...) {
  if (!inherits(x, "nda")) stop("x must inherit from 'nda'.", call. = FALSE)
  similarity <- as.matrix(x$R)
  values <- similarity[similarity > 0]
  threshold <- cuts
  if (isTRUE(x$covar) && length(values)) {
    threshold <- min(values) + cuts * (max(values) - min(values))
  }
  if (length(threshold) == 1L) {
    similarity[similarity < threshold] <- 0
  } else if (length(threshold) == 2L) {
    similarity[similarity < min(threshold) | similarity > max(threshold)] <- 0
  } else stop("cuts must have length one or two.", call. = FALSE)
  directed <- !isSymmetric(similarity)
  graph <- igraph::graph_from_adjacency_matrix(similarity,
    mode = if (directed) "directed" else "undirected", weighted = TRUE,
    diag = FALSE)
  membership <- x$membership
  palette <- grDevices::hcl.colors(max(1L, max(membership)), "Dark 3")
  vertex_color <- ifelse(membership > 0, palette[pmax(1L, membership)], "grey80")
  shape <- rep(shapes, length.out = igraph::vcount(graph))
  edge_pairs <- igraph::as_edgelist(graph, names = FALSE)
  association_values <- if (nrow(edge_pairs)) x$correlation[edge_pairs] else numeric()
  scale_colors <- function(z) {
    pal <- grDevices::colorRampPalette(edge.colors)(101)
    if (!length(z)) return(character())
    if (isTRUE(x$covar)) {
      low <- min(z, na.rm = TRUE); middle <- stats::median(z, na.rm = TRUE)
      high <- max(z, na.rm = TRUE)
      index <- ifelse(z <= middle,
        1 + 49 * (z - low) / max(middle - low, .Machine$double.eps),
        51 + 50 * (z - middle) / max(high - middle, .Machine$double.eps))
    } else index <- 1 + 100 * (pmax(-1, pmin(1, z)) + 1) / 2
    pal[pmax(1L, pmin(101L, round(index)))]
  }
  edge_color <- if (edge.color.links) scale_colors(association_values) else "grey55"
  edge_label <- if (show_weights) format(round(association_values, 3), trim = TRUE) else ""
  igraph::V(graph)$color <- vertex_color
  igraph::V(graph)$shape <- shape
  igraph::V(graph)$label.color <- vertex.label.color
  igraph::V(graph)$label.cex <- vertex.label.size / 12
  igraph::E(graph)$color <- edge_color
  igraph::E(graph)$label <- edge_label
  igraph::E(graph)$label.color <- edge.label.color
  igraph::E(graph)$label.cex <- edge.label.size / 12
  igraph::E(graph)$width <- 1 + edgescale * 4 * igraph::E(graph)$weight

  if (interactive) {
    nodes <- data.frame(
      id = seq_len(igraph::vcount(graph)),
      label = igraph::V(graph)$name,
      color = vertex_color,
      shape = gsub("rectangle", "box", gsub("circle", "ellipse", shape)),
      font.size = vertex.label.size,
      stringsAsFactors = FALSE)
    edges <- data.frame(
      from = edge_pairs[, 1L], to = edge_pairs[, 2L],
      width = 1 + edgescale * 4 * igraph::E(graph)$weight,
      color = edge_color, label = edge_label,
      arrows = if (directed) "to" else "",
      font.size = edge.label.size, font.color = edge.label.color,
      stringsAsFactors = FALSE)
    widget <- do.call(visNetwork::visNetwork,
      c(list(nodes = nodes, edges = edges, height = "800px", width = "100%"),
        visnetwork_args$network %||% list()))
    widget <- do.call(visNetwork::visNodes,
      c(list(graph = widget), visnetwork_args$nodes %||% list()))
    widget <- do.call(visNetwork::visEdges,
      c(list(graph = widget), visnetwork_args$edges %||% list()))
    widget <- do.call(visNetwork::visOptions,
      c(list(graph = widget, highlightNearest = TRUE, selectedBy = "label"),
        visnetwork_args$options %||% list()))
    widget <- do.call(visNetwork::visIgraphLayout,
      c(list(graph = widget, layout = paste0("layout_", layout), physics = TRUE),
        visnetwork_args$layout %||% list()))
  } else widget <- NULL
  out <- list(graph = graph, widget = widget, interactive = interactive,
              layout = layout, labeldist = labeldist,
              igraph_args = c(igraph_args, list(...)))
  class(out) <- c("nda_plot", "list")
  out
}

#' @export
print.nda_plot <- function(x, ...) {
  if (x$interactive) {
    print(x$widget)
  } else {
    layout_value <- if (is.character(x$layout)) {
      layout_name <- switch(x$layout, nicely = "layout_nicely", circle = "layout_in_circle",
                            sphere = "layout_on_sphere", grid = "layout_on_grid",
                            paste0("layout_with_", x$layout))
      fun <- get(layout_name, envir = asNamespace("igraph"), inherits = FALSE)
      fun(x$graph)
    } else x$layout
    args <- c(list(x = x$graph, layout = layout_value,
                   vertex.label.dist = x$labeldist), x$igraph_args, list(...))
    do.call(igraph::plot.igraph, args)
  }
  invisible(x)
}
