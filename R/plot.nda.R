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
plot.nda <- function(x, cuts = 0.3, interactive = FALSE, edgescale = 1,
                     labeldist = -1.5, show_weights = FALSE,
                     shapes = NULL, shape = "automatic",
                     static.engine = c("grid", "igraph"),
                     edge.color.links = TRUE, edge.colors = NULL,
                     positive.color = "#009E73", negative.color = "#C62828",
                     neutral.color = "grey65", positive.lty = 1,
                     negative.lty = 1, neutral.lty = 1,
                     edge.width.range = c(1, 6),
                     edge.label.size = 10, edge.label.color = "#333333",
                     edge.label.background = "#FFFFFFE6",
                     edge.label.offset = 0.022,
                     vertex.label.size = 12,
                     vertex.label.color = "automatic",
                     node.border.color = "#25313C",
                     layout = "nicely", igraph_args = list(),
                     visnetwork_args = list(), ...) {
  if (!inherits(x, "nda")) stop("x must inherit from 'nda'.", call. = FALSE)
  static.engine <- match.arg(static.engine)
  similarity <- x$R %||% x$correlation
  if (is.null(similarity)) {
    stop("x does not contain an association matrix in $R or $correlation.",
         call. = FALSE)
  }
  similarity <- as.matrix(similarity)
  if (!is.numeric(similarity) || nrow(similarity) != ncol(similarity)) {
    stop("The association matrix must be numeric and square.", call. = FALSE)
  }
  if (is.null(colnames(similarity))) {
    colnames(similarity) <- rownames(similarity) <-
      paste0("V", seq_len(ncol(similarity)))
  } else if (is.null(rownames(similarity))) {
    rownames(similarity) <- colnames(similarity)
  } else if (setequal(rownames(similarity), colnames(similarity))) {
    similarity <- similarity[colnames(similarity), colnames(similarity), drop = FALSE]
  }
  raw_similarity <- similarity
  similarity[!is.finite(similarity)] <- 0
  values <- similarity[is.finite(similarity) & similarity > 0]
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
  if (is.null(membership)) membership <- rep(0L, ncol(similarity))
  if (!is.null(names(membership)) &&
      all(colnames(similarity) %in% names(membership))) {
    membership <- membership[colnames(similarity)]
  }
  membership <- suppressWarnings(as.integer(membership))
  if (length(membership) != ncol(similarity)) {
    stop("x$membership must contain one value per association-matrix variable.",
         call. = FALSE)
  }
  membership[is.na(membership) | membership < 0L] <- 0L
  groups <- if (length(membership)) max(membership, 0L) else 0L
  palette <- grDevices::hcl.colors(max(1L, groups), "Set 2")
  vertex_color <- ifelse(membership > 0L,
                         palette[pmax(1L, membership)], "grey88")
  if (!is.null(shapes)) shape <- shapes
  shape <- .nda_plot_shapes(shape, igraph::vcount(graph), "rectangle")
  edge_pairs <- igraph::as_edgelist(graph, names = FALSE)
  association <- x$correlation
  if (is.null(association)) association <- raw_similarity
  association <- as.matrix(association)
  if (!all(dim(association) == dim(raw_similarity))) association <- raw_similarity
  if (!is.null(rownames(association)) && !is.null(colnames(association)) &&
      setequal(rownames(association), colnames(raw_similarity)) &&
      setequal(colnames(association), colnames(raw_similarity))) {
    association <- association[colnames(raw_similarity), colnames(raw_similarity),
                               drop = FALSE]
  }
  association_values <- if (nrow(edge_pairs)) association[edge_pairs] else numeric()
  graph_weights <- if (nrow(edge_pairs)) igraph::E(graph)$weight else numeric()
  association_values[!is.finite(association_values)] <-
    graph_weights[!is.finite(association_values)]
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
  if (!is.null(edge.colors)) {
    if (!is.character(edge.colors) || length(edge.colors) < 2L) {
      stop("edge.colors must contain at least two colours.", call. = FALSE)
    }
    edge_color <- if (edge.color.links) scale_colors(association_values) else neutral.color
  } else {
    edge_color <- if (edge.color.links) {
      ifelse(association_values < 0, negative.color,
             ifelse(association_values > 0, positive.color, neutral.color))
    } else rep(neutral.color, length(association_values))
  }
  edge_lty <- ifelse(association_values < 0, negative.lty,
                     ifelse(association_values > 0, positive.lty, neutral.lty))
  edge_width <- .nda_edge_widths(association_values, edge.width.range)
  if (!is.null(edgescale) && is.numeric(edgescale) && length(edgescale) == 1L &&
      is.finite(edgescale) && edgescale >= 0) edge_width <- edge_width * edgescale
  edge_label <- if (show_weights) {
    format(round(association_values, 3), trim = TRUE)
  } else rep("", length(association_values))
  node_label_color <- if (identical(vertex.label.color, "automatic")) {
    .nda_contrast_color(vertex_color)
  } else rep(vertex.label.color, length.out = length(vertex_color))
  igraph::V(graph)$color <- vertex_color
  igraph::V(graph)$shape <- shape
  igraph::V(graph)$label.color <- node_label_color
  igraph::V(graph)$label.cex <- vertex.label.size / 12
  igraph::E(graph)$color <- edge_color
  igraph::E(graph)$lty <- edge_lty
  igraph::E(graph)$label <- edge_label
  igraph::E(graph)$label.color <- edge.label.color
  igraph::E(graph)$label.cex <- edge.label.size / 12
  igraph::E(graph)$width <- edge_width

  nodes <- data.frame(
    id = seq_len(igraph::vcount(graph)),
    label = igraph::V(graph)$name,
    color = vertex_color,
    shape = shape,
    label_color = node_label_color,
    stringsAsFactors = FALSE
  )
  render_edges <- data.frame(
    from = if (nrow(edge_pairs)) edge_pairs[, 1L] else integer(),
    to = if (nrow(edge_pairs)) edge_pairs[, 2L] else integer(),
    width = edge_width,
    color = edge_color,
    lty = edge_lty,
    label = edge_label,
    stringsAsFactors = FALSE
  )
  coordinates <- .nda_igraph_layout(graph, layout)
  static_grob <- .nda_network_grob(nodes, render_edges, coordinates,
    directed = directed, vertex.label.size = vertex.label.size,
    edge.label.size = edge.label.size,
    edge.label.color = edge.label.color,
    edge.label.background = edge.label.background,
    edge.label.offset = edge.label.offset,
    node.border.color = node.border.color)

  if (interactive) {
    vis_nodes <- nodes
    vis_nodes$shape <- gsub("rectangle", "box",
      gsub("circle", "ellipse", vis_nodes$shape))
    vis_nodes$font.size <- vertex.label.size
    vis_nodes$font.color <- vis_nodes$label_color
    vis_nodes$borderWidth <- 1.5
    vis_nodes$title <- ifelse(membership > 0L,
      paste0(vis_nodes$label, " - community ", membership),
      paste0(vis_nodes$label, " - excluded"))
    edges <- data.frame(
      from = edge_pairs[, 1L], to = edge_pairs[, 2L],
      width = edge_width,
      color = edge_color, label = edge_label,
      dashes = edge_lty != 1,
      arrows = if (directed) "to" else "",
      font.size = edge.label.size, font.color = edge.label.color,
      stringsAsFactors = FALSE)
    edges$font.background <- edge.label.background
    widget <- do.call(visNetwork::visNetwork,
      .nda_merge_args(list(nodes = vis_nodes, edges = edges,
        height = "800px", width = "100%"),
        visnetwork_args$network %||% list()))
    widget <- do.call(visNetwork::visNodes,
      .nda_merge_args(list(graph = widget, margin = 12,
        shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.16)",
          size = 5, x = 1, y = 2)),
        visnetwork_args$nodes %||% list()))
    widget <- do.call(visNetwork::visEdges,
      .nda_merge_args(list(graph = widget,
        font = list(background = edge.label.background, align = "middle",
          strokeWidth = 0),
        smooth = list(enabled = TRUE, type = "dynamic", roundness = 0.12)),
        visnetwork_args$edges %||% list()))
    widget <- do.call(visNetwork::visOptions,
      .nda_merge_args(list(graph = widget, highlightNearest = TRUE,
        selectedBy = "label"), visnetwork_args$options %||% list()))
    widget <- do.call(visNetwork::visInteraction,
      .nda_merge_args(list(graph = widget, dragNodes = TRUE,
        dragView = TRUE, zoomView = TRUE, hideEdgesOnDrag = FALSE),
        visnetwork_args$interaction %||% list()))
    vis_layout <- if (is.character(layout)) switch(layout,
      nicely = "layout_nicely", circle = "layout_in_circle",
      sphere = "layout_on_sphere", grid = "layout_on_grid",
      paste0("layout_with_", layout)) else "layout_nicely"
    widget <- do.call(visNetwork::visIgraphLayout,
      .nda_merge_args(list(graph = widget, layout = vis_layout,
        physics = TRUE), visnetwork_args$layout %||% list()))
  } else widget <- NULL
  out <- list(graph = graph, nodes = nodes, edges = render_edges,
              widget = widget, interactive = interactive,
              static.engine = static.engine, static_grob = static_grob,
              layout = layout, coordinates = coordinates,
              labeldist = labeldist,
              igraph_args = c(igraph_args, list(...)))
  class(out) <- c("nda_plot", "list")
  out
}

#' @export
print.nda_plot <- function(x, ...) {
  if (x$interactive) {
    print(x$widget)
  } else if (identical(x$static.engine, "grid")) {
    .nda_draw_network_grob(x$static_grob)
  } else {
    args <- c(list(x = x$graph, layout = x$coordinates,
                   vertex.label.dist = x$labeldist), x$igraph_args, list(...))
    do.call(igraph::plot.igraph, args)
  }
  invisible(x)
}
