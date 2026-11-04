# Internal publication-quality network drawing helpers.

.nda_centered_positions <- function(n, span = 0.82) {
  n <- suppressWarnings(as.integer(n)[1L])
  if (!length(n) || is.na(n) || n <= 0L) return(numeric())
  if (n == 1L) return(0)
  seq(span, -span, length.out = n)
}

.nda_contrast_color <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  luminance <- 0.2126 * rgb[1L, ] + 0.7152 * rgb[2L, ] +
    0.0722 * rgb[3L, ]
  ifelse(luminance < 0.48, "white", "#17202A")
}

.nda_lighten_color <- function(fill, amount = 0, alpha = 1) {
  if (!is.numeric(amount) || length(amount) != 1L || !is.finite(amount) ||
      amount < 0 || amount > 1) {
    stop("amount must be a number between zero and one.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a number between zero and one.", call. = FALSE)
  }
  rgba <- grDevices::col2rgb(fill, alpha = TRUE) / 255
  rgba[seq_len(3L), ] <- rgba[seq_len(3L), , drop = FALSE] +
    (1 - rgba[seq_len(3L), , drop = FALSE]) * amount
  grDevices::rgb(rgba[1L, ], rgba[2L, ], rgba[3L, ],
                 alpha = rgba[4L, ] * alpha)
}

.nda_merge_args <- function(defaults, overrides = NULL) {
  if (is.null(overrides) || !length(overrides)) return(defaults)
  if (is.null(names(overrides)) || any(!nzchar(names(overrides)))) {
    stop("Forwarded plotting arguments must be named.", call. = FALSE)
  }
  utils::modifyList(defaults, overrides, keep.null = TRUE)
}

.nda_igraph_layout <- function(graph, layout) {
  node_count <- igraph::vcount(graph)
  if (node_count == 1L) return(matrix(c(0.5, 0.5), nrow = 1L))
  coordinates <- if (is.character(layout)) {
    layout_name <- switch(layout,
      nicely = "layout_nicely",
      circle = "layout_in_circle",
      sphere = "layout_on_sphere",
      grid = "layout_on_grid",
      paste0("layout_with_", layout))
    fun <- get(layout_name, envir = asNamespace("igraph"), inherits = FALSE)
    fun(graph)
  } else layout
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) != node_count || ncol(coordinates) < 2L) {
    stop("layout must provide two coordinates per graph node.", call. = FALSE)
  }
  igraph::norm_coords(coordinates[, 1:2, drop = FALSE],
    xmin = 0.08, xmax = 0.92, ymin = 0.08, ymax = 0.92)
}

.nda_network_grob <- function(nodes, edges, coordinates, directed = FALSE,
                              vertex.label.size = 11,
                              edge.label.size = 9,
                              edge.label.color = "#333333",
                              edge.label.background = "#FFFFFFE6",
                              edge.label.offset = 0.022,
                              node.border.color = "#25313C") {
  coordinates <- as.matrix(coordinates)
  if (nrow(coordinates) != nrow(nodes) || ncol(coordinates) < 2L) {
    stop("coordinates must contain two positions per node.", call. = FALSE)
  }
  coordinates[, 1L] <- pmin(0.94, pmax(0.06, coordinates[, 1L]))
  coordinates[, 2L] <- pmin(0.94, pmax(0.06, coordinates[, 2L]))
  node_index <- stats::setNames(seq_len(nrow(nodes)), as.character(nodes$id))
  grobs <- list()
  add_grob <- function(grob) grobs[[length(grobs) + 1L]] <<- grob

  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      from_index <- unname(node_index[as.character(edges$from[i])])
      to_index <- unname(node_index[as.character(edges$to[i])])
      if (!length(from_index) || !length(to_index) ||
          is.na(from_index) || is.na(to_index)) next
      x0 <- coordinates[from_index, 1L]
      y0 <- coordinates[from_index, 2L]
      x1 <- coordinates[to_index, 1L]
      y1 <- coordinates[to_index, 2L]
      dx <- x1 - x0
      dy <- y1 - y0
      distance <- sqrt(dx^2 + dy^2)
      if (!is.finite(distance) || distance <= .Machine$double.eps) next
      ux <- dx / distance
      uy <- dy / distance
      start_gap <- 0.045
      end_gap <- if (directed) 0.064 else 0.045
      add_grob(grid::segmentsGrob(
        x0 = grid::unit(x0 + ux * start_gap, "npc"),
        y0 = grid::unit(y0 + uy * start_gap, "npc"),
        x1 = grid::unit(x1 - ux * end_gap, "npc"),
        y1 = grid::unit(y1 - uy * end_gap, "npc"),
        arrow = if (directed) grid::arrow(type = "closed",
          length = grid::unit(2.4, "mm")) else NULL,
        gp = grid::gpar(col = edges$color[i],
          lwd = max(0.5, edges$width[i]), lty = edges$lty[i],
          lineend = "round", linejoin = "round")
      ))
      label <- as.character(edges$label[i])
      if (length(label) && !is.na(label) && nzchar(label)) {
        direction <- if (i %% 2L) 1 else -1
        label_x <- (x0 + x1) / 2 - uy * edge.label.offset * direction
        label_y <- (y0 + y1) / 2 + ux * edge.label.offset * direction
        label_grob <- grid::textGrob(label,
          x = grid::unit(label_x, "npc"), y = grid::unit(label_y, "npc"),
          gp = grid::gpar(col = edge.label.color,
            fontsize = edge.label.size, fontfamily = "sans"))
        add_grob(grid::roundrectGrob(
          x = grid::unit(label_x, "npc"), y = grid::unit(label_y, "npc"),
          width = grid::grobWidth(label_grob) + grid::unit(2.8, "mm"),
          height = grid::grobHeight(label_grob) + grid::unit(1.6, "mm"),
          r = grid::unit(1.2, "mm"),
          gp = grid::gpar(fill = edge.label.background, col = NA)))
        add_grob(label_grob)
      }
    }
  }

  for (i in seq_len(nrow(nodes))) {
    x <- coordinates[i, 1L]
    y <- coordinates[i, 2L]
    label <- as.character(nodes$label[i])
    shape <- tolower(as.character(nodes$shape[i]))
    label_color <- nodes$label_color[i]
    if (!length(label_color) || is.na(label_color)) {
      label_color <- .nda_contrast_color(nodes$color[i])
    }
    label_grob <- grid::textGrob(label,
      x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
      gp = grid::gpar(col = label_color, fontsize = vertex.label.size,
        fontfamily = "sans", fontface = "plain"))
    if (shape %in% c("circle", "ellipse")) {
      radius <- max(0.052, min(0.082, 0.046 + 0.004 * nchar(label)))
      add_grob(grid::circleGrob(
        x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
        r = grid::unit(radius, "snpc"),
        gp = grid::gpar(fill = nodes$color[i], col = node.border.color,
          lwd = 1.2)))
    } else {
      square <- shape %in% c("square", "squared")
      width <- if (square) grid::unit(0.09, "snpc") else
        grid::grobWidth(label_grob) + grid::unit(7, "mm")
      height <- if (square) grid::unit(0.09, "snpc") else
        grid::grobHeight(label_grob) + grid::unit(4.5, "mm")
      add_grob(grid::roundrectGrob(
        x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
        width = width, height = height,
        r = grid::unit(1.8, "mm"),
        gp = grid::gpar(fill = nodes$color[i], col = node.border.color,
          lwd = 1.2)))
    }
    add_grob(label_grob)
  }
  grid::gTree(children = do.call(grid::gList, grobs))
}

.nda_draw_network_grob <- function(grob) {
  grid::grid.newpage()
  grid::grid.draw(grob)
  invisible(grob)
}


