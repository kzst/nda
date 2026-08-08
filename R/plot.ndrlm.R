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


#' Plot an NDRLM path network
#' @export
plot.ndrlm <- function(x, sig = 0.05, interactive = FALSE,
                       show_edge_labels = TRUE, edge.label.size = 9,
                       vertex.label.size = 11, edgescale = 2,
                       igraph_args = list(), visnetwork_args = list(), ...) {
  if (!inherits(x, "ndrlm")) stop("x must inherit from 'ndrlm'.", call. = FALSE)
  node_rows <- list()
  add_nodes <- function(id, label, role, column) {
    data.frame(id = id, label = label, role = role, column = column,
               stringsAsFactors = FALSE)
  }
  node_rows[[length(node_rows) + 1L]] <- add_nodes(
    paste0("X:", colnames(x$X)), colnames(x$X), "indicator_in", 0)
  if (x$latents %in% c("in", "both")) {
    names_in <- paste0("NDAin", seq_len(x$NDAin$factors))
    node_rows[[length(node_rows) + 1L]] <- add_nodes(
      paste0("Lin:", names_in), names_in, "latent_in", 1)
  }
  if (x$latents %in% c("out", "both")) {
    names_out <- paste0("NDAout", seq_len(x$NDAout$factors))
    node_rows[[length(node_rows) + 1L]] <- add_nodes(
      paste0("Lout:", names_out), names_out, "latent_out", 2)
  }
  node_rows[[length(node_rows) + 1L]] <- add_nodes(
    paste0("Y:", colnames(x$Y)), colnames(x$Y), "indicator_out", 3)
  nodes <- do.call(rbind, node_rows)
  role_colors <- c(indicator_in = "#9ecae1", latent_in = "#3182bd",
                   latent_out = "#e6550d", indicator_out = "#fdae6b")
  nodes$color <- unname(role_colors[nodes$role])
  nodes$shape <- ifelse(grepl("latent", nodes$role), "ellipse", "box")
  nodes$y <- 0
  for (column in unique(nodes$column)) {
    index <- which(nodes$column == column)
    nodes$y[index] <- seq(1, -1, length.out = length(index))
  }
  edges <- data.frame(from = character(), to = character(), weight = numeric(),
                      color = character(), dashes = logical(),
                      stringsAsFactors = FALSE)
  add_edge <- function(from, to, weight, color, dashes) {
    data.frame(from = from, to = to, weight = as.numeric(weight),
               color = color, dashes = dashes, stringsAsFactors = FALSE)
  }
  if (x$latents %in% c("in", "both")) {
    for (j in which(x$NDAin$membership > 0L)) {
      g <- x$NDAin$membership[j]
      edges <- rbind(edges, add_edge(paste0("X:", colnames(x$X)[j]),
        paste0("Lin:NDAin", g), x$NDAin$loadings[j, g], "grey55", TRUE))
    }
  }
  if (x$latents %in% c("out", "both")) {
    for (j in which(x$NDAout$membership > 0L)) {
      g <- x$NDAout$membership[j]
      edges <- rbind(edges, add_edge(paste0("Y:", colnames(x$Y)[j]),
        paste0("Lout:NDAout", g), x$NDAout$loadings[j, g], "grey55", TRUE))
    }
  }
  independent_id <- function(name) {
    if (startsWith(name, "NDAin")) paste0("Lin:", name) else paste0("X:", name)
  }
  dependent_id <- function(name) {
    if (startsWith(name, "NDAout")) paste0("Lout:", name) else paste0("Y:", name)
  }
  for (i in seq_along(x$fits)) {
    fit <- x$fits[[i]]
    coef <- fit$coefficients
    p_value <- fit$p_values
    for (j in seq_along(coef)) {
      if (!is.finite(coef[j]) || coef[j] == 0) next
      if (is.finite(p_value[j]) && p_value[j] >= sig) next
      edges <- rbind(edges, add_edge(independent_id(names(coef)[j]),
        dependent_id(names(x$fits)[i]), coef[j],
        if (coef[j] < 0) "darkred" else "darkblue", FALSE))
    }
  }
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  igraph::V(graph)$color <- nodes$color
  igraph::V(graph)$shape <- ifelse(nodes$shape == "box", "rectangle", "circle")
  igraph::V(graph)$label <- nodes$label
  igraph::V(graph)$label.cex <- vertex.label.size / 12
  igraph::E(graph)$color <- edges$color
  igraph::E(graph)$lty <- ifelse(edges$dashes, 2, 1)
  igraph::E(graph)$width <- 1 + edgescale * abs(edges$weight)
  igraph::E(graph)$label <- if (show_edge_labels) round(edges$weight, 3) else ""
  igraph::E(graph)$label.cex <- edge.label.size / 12
  layout_matrix <- as.matrix(nodes[, c("column", "y")])

  if (interactive) {
    vis_nodes <- data.frame(id = nodes$id, label = nodes$label,
      color = nodes$color, shape = nodes$shape, x = nodes$column * 260,
      y = nodes$y * 500, fixed = FALSE, font.size = vertex.label.size,
      stringsAsFactors = FALSE)
    vis_edges <- edges
    vis_edges$width <- 1 + edgescale * abs(vis_edges$weight)
    vis_edges$arrows <- "to"
    vis_edges$label <- if (show_edge_labels) {
      format(round(vis_edges$weight, 3), trim = TRUE)
    } else ""
    vis_edges$font.size <- edge.label.size
    widget <- do.call(visNetwork::visNetwork,
      c(list(nodes = vis_nodes, edges = vis_edges, height = "850px", width = "100%"),
        visnetwork_args$network %||% list()))
    widget <- do.call(visNetwork::visOptions,
      c(list(graph = widget, highlightNearest = TRUE, selectedBy = "label"),
        visnetwork_args$options %||% list()))
    widget <- do.call(visNetwork::visPhysics,
      c(list(graph = widget, enabled = FALSE), visnetwork_args$physics %||% list()))
  } else widget <- NULL
  out <- list(graph = graph, layout = layout_matrix, interactive = interactive,
              widget = widget, igraph_args = c(igraph_args, list(...)))
  class(out) <- c("ndrlm_plot", "list")
  out
}

#' @export
print.ndrlm_plot <- function(x, ...) {
  if (x$interactive) print(x$widget) else {
    do.call(igraph::plot.igraph,
      c(list(x = x$graph, layout = x$layout, edge.arrow.size = 0.35,
             vertex.size = 26), x$igraph_args, list(...)))
  }
  invisible(x)
}
