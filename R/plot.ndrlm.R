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
                       vertex.label.size = 11, edgescale = 1,
                       shape = "automatic",
                       positive.color = "#009E73", negative.color = "#C62828",
                       loading.color = "grey60", positive.lty = 1,
                       negative.lty = 1, loading.lty = 1,
                       edge.width.range = c(1, 6),
                       igraph_args = list(), visnetwork_args = list(), ...) {
  if (!inherits(x, "ndrlm")) stop("x must inherit from 'ndrlm'.", call. = FALSE)
  if (is.null(x$X) || is.null(x$Y) || is.null(x$fits)) {
    stop("x must contain $X, $Y, and $fits.", call. = FALSE)
  }
  input_names <- colnames(x$X) %||% paste0("X", seq_len(ncol(x$X)))
  output_names <- colnames(x$Y) %||% paste0("Y", seq_len(ncol(x$Y)))
  latents <- x$latents %||% if (!is.null(x$NDAin) && !is.null(x$NDAout)) {
    "both"
  } else if (!is.null(x$NDAin)) {
    "in"
  } else if (!is.null(x$NDAout)) {
    "out"
  } else {
    "none"
  }
  latents <- tolower(as.character(latents)[1L])
  if (!latents %in% c("in", "out", "both", "none")) latents <- "none"
  factor_count <- function(model) {
    if (is.null(model)) return(0L)
    value <- model$factors
    if (is.null(value) && !is.null(model$scores)) value <- ncol(model$scores)
    if (is.null(value) && length(model$membership)) {
      value <- max(suppressWarnings(as.integer(model$membership)), 0L,
                   na.rm = TRUE)
    }
    value <- suppressWarnings(as.integer(value)[1L])
    if (!length(value) || is.na(value) || value < 0L) 0L else value
  }
  factors_in <- if (latents %in% c("in", "both")) factor_count(x$NDAin) else 0L
  factors_out <- if (latents %in% c("out", "both")) factor_count(x$NDAout) else 0L
  node_rows <- list()
  add_nodes <- function(id, label, role, column) {
    data.frame(id = id, label = label, role = role, column = column,
               stringsAsFactors = FALSE)
  }
  node_rows[[length(node_rows) + 1L]] <- add_nodes(
    paste0("X:", input_names), input_names, "indicator_in", 0)
  if (factors_in > 0L) {
    names_in <- paste0("NDAin", seq_len(factors_in))
    node_rows[[length(node_rows) + 1L]] <- add_nodes(
      paste0("Lin:", names_in), names_in, "latent_in", 1)
  }
  if (factors_out > 0L) {
    names_out <- paste0("NDAout", seq_len(factors_out))
    node_rows[[length(node_rows) + 1L]] <- add_nodes(
      paste0("Lout:", names_out), names_out, "latent_out", 2)
  }
  node_rows[[length(node_rows) + 1L]] <- add_nodes(
    paste0("Y:", output_names), output_names, "indicator_out", 3)
  nodes <- do.call(rbind, node_rows)
  role_colors <- c(indicator_in = "#9ecae1", latent_in = "#3182bd",
                   latent_out = "#CC79A7", indicator_out = "#E69F00")
  nodes$color <- unname(role_colors[nodes$role])
  automatic_shape <- ifelse(grepl("latent", nodes$role), "circle", "rectangle")
  nodes$shape <- .nda_plot_shapes(shape, nrow(nodes), automatic_shape)
  nodes$y <- 0
  for (column in unique(nodes$column)) {
    index <- which(nodes$column == column)
    nodes$y[index] <- seq(1, -1, length.out = length(index))
  }
  edges <- data.frame(from = character(), to = character(), weight = numeric(),
                      color = character(), lty = numeric(), role = character(),
                      stringsAsFactors = FALSE)
  add_edge <- function(from, to, weight, color, lty, role) {
    data.frame(from = from, to = to, weight = as.numeric(weight),
               color = color, lty = lty, role = role,
               stringsAsFactors = FALSE)
  }
  normalized_membership <- function(model, variable_names, factors) {
    value <- model$membership %||% rep(0L, length(variable_names))
    if (!is.null(names(value)) && all(variable_names %in% names(value))) {
      value <- value[variable_names]
    }
    value <- suppressWarnings(as.integer(value))
    if (length(value) != length(variable_names)) return(rep(0L, length(variable_names)))
    value[is.na(value) | value < 1L | value > factors] <- 0L
    value
  }
  normalized_loadings <- function(model, variable_names, factors) {
    value <- as.matrix(model$loadings)
    if (!length(value) || nrow(value) != length(variable_names) ||
        ncol(value) < factors) {
      return(matrix(NA_real_, nrow = length(variable_names), ncol = factors,
                    dimnames = list(variable_names, NULL)))
    }
    if (!is.null(rownames(value)) && all(variable_names %in% rownames(value))) {
      value <- value[variable_names, , drop = FALSE]
    }
    value
  }
  membership_in <- if (factors_in > 0L) {
    normalized_membership(x$NDAin, input_names, factors_in)
  } else rep(0L, length(input_names))
  membership_out <- if (factors_out > 0L) {
    normalized_membership(x$NDAout, output_names, factors_out)
  } else rep(0L, length(output_names))
  total_factors <- factors_in + factors_out
  community_palette <- if (total_factors > 0L) {
    grDevices::hcl.colors(total_factors, "Dark 3")
  } else character()
  if (factors_in > 0L) {
    for (j in which(membership_in > 0L)) {
      nodes$color[nodes$id == paste0("X:", input_names[j])] <-
        community_palette[membership_in[j]]
    }
    for (g in seq_len(factors_in)) {
      nodes$color[nodes$id == paste0("Lin:NDAin", g)] <- community_palette[g]
    }
  }
  if (factors_out > 0L) {
    offset <- factors_in
    for (j in which(membership_out > 0L)) {
      nodes$color[nodes$id == paste0("Y:", output_names[j])] <-
        community_palette[offset + membership_out[j]]
    }
    for (g in seq_len(factors_out)) {
      nodes$color[nodes$id == paste0("Lout:NDAout", g)] <-
        community_palette[offset + g]
    }
  }
  if (factors_in > 0L) {
    loadings_in <- normalized_loadings(x$NDAin, input_names, factors_in)
    for (j in which(membership_in > 0L)) {
      g <- membership_in[j]
      weight <- loadings_in[j, g]
      if (!is.finite(weight)) next
      edges <- rbind(edges, add_edge(paste0("X:", input_names[j]),
        paste0("Lin:NDAin", g), weight, loading.color, loading.lty,
        "loading"))
    }
  }
  if (factors_out > 0L) {
    loadings_out <- normalized_loadings(x$NDAout, output_names, factors_out)
    for (j in which(membership_out > 0L)) {
      g <- membership_out[j]
      weight <- loadings_out[j, g]
      if (!is.finite(weight)) next
      edges <- rbind(edges, add_edge(paste0("Y:", output_names[j]),
        paste0("Lout:NDAout", g), weight, loading.color, loading.lty,
        "loading"))
    }
  }
  independent_id <- function(name) {
    if (length(name) && !is.na(name) && startsWith(name, "NDAin")) {
      paste0("Lin:", name)
    } else paste0("X:", name)
  }
  dependent_id <- function(name) {
    if (length(name) && !is.na(name) && startsWith(name, "NDAout")) {
      paste0("Lout:", name)
    } else paste0("Y:", name)
  }
  for (i in seq_along(x$fits)) {
    fit <- x$fits[[i]]
    coef <- .ndrlm_fit_coefficients(fit)
    p_value <- .ndrlm_fit_p_values(fit)
    if (length(p_value) != length(coef)) p_value <- rep(NA_real_, length(coef))
    dependent_name <- names(x$fits)[i]
    if (is.null(dependent_name) || is.na(dependent_name) ||
        !nzchar(dependent_name)) {
      dependent_name <- .ndrlm_fit_response(fit) %||% output_names[i]
    }
    for (j in seq_along(coef)) {
      if (!is.finite(coef[j]) || coef[j] == 0) next
      if (is.finite(p_value[j]) && p_value[j] >= sig) next
      from <- independent_id(names(coef)[j])
      to <- dependent_id(dependent_name)
      if (!from %in% nodes$id || !to %in% nodes$id) next
      edges <- rbind(edges, add_edge(from, to, coef[j],
        if (coef[j] < 0) negative.color else positive.color,
        if (coef[j] < 0) negative.lty else positive.lty, "path"))
    }
  }
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
  igraph::V(graph)$color <- nodes$color
  igraph_shape <- gsub("box", "rectangle", gsub("ellipse", "circle", nodes$shape))
  igraph::V(graph)$shape <- igraph_shape
  igraph::V(graph)$label <- nodes$label
  igraph::V(graph)$label.cex <- vertex.label.size / 12
  igraph::E(graph)$color <- edges$color
  igraph::E(graph)$lty <- edges$lty
  edge_width <- numeric(nrow(edges))
  for (edge_role in unique(edges$role)) {
    index <- which(edges$role == edge_role)
    edge_width[index] <- .nda_edge_widths(edges$weight[index], edge.width.range)
  }
  if (!is.null(edgescale) && is.numeric(edgescale) && length(edgescale) == 1L &&
      is.finite(edgescale) && edgescale >= 0) edge_width <- edge_width * edgescale
  igraph::E(graph)$width <- edge_width
  igraph::E(graph)$label <- if (show_edge_labels) round(edges$weight, 3) else ""
  igraph::E(graph)$label.cex <- edge.label.size / 12
  layout_matrix <- as.matrix(nodes[, c("column", "y")])

  if (interactive) {
    vis_shape <- gsub("rectangle", "box", gsub("circle", "ellipse", nodes$shape))
    vis_nodes <- data.frame(id = nodes$id, label = nodes$label,
      color = nodes$color, shape = vis_shape, x = nodes$column * 260,
      y = nodes$y * 500, fixed = FALSE, font.size = vertex.label.size,
      stringsAsFactors = FALSE)
    vis_edges <- edges
    vis_edges$width <- edge_width
    vis_edges$dashes <- vis_edges$lty != 1
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
