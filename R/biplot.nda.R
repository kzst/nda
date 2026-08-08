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


#' Biplot NDA scores and loadings
#' @export
biplot.nda <- function(x, main = NULL, dims = c(1, 2), view_3D = FALSE,
                       triples = NULL, interactive = FALSE, angle = 45,
                       point_cex = 0.7, label_cex = 0.9, ...) {
  if (!inherits(x, "nda")) stop("x must inherit from 'nda'.", call. = FALSE)
  if (is.null(x$scores)) {
    stop("A biplot requires a model fitted to row-level data.", call. = FALSE)
  }
  if (view_3D && x$factors < 3L) stop("At least three factors are required for a 3D view.", call. = FALSE)
  if (is.null(triples) && view_3D) {
    triples <- utils::combn(seq_len(x$factors), 3L, simplify = FALSE)
  }
  one_dimensional <- x$factors == 1L && !view_3D
  if (!one_dimensional && !view_3D &&
      (length(dims) != 2L || any(!dims %in% seq_len(x$factors)))) {
    stop("dims must select two fitted factors.", call. = FALSE)
  }
  out <- list(scores = x$scores, loadings = x$loadings, main = main,
              dims = dims, view_3D = view_3D, triples = triples,
              interactive = interactive, angle = angle,
              point_cex = point_cex, label_cex = label_cex,
              one_dimensional = one_dimensional, dots = list(...))
  class(out) <- c("nda_biplot", "list")
  out
}

.nda_draw_biplot_2d <- function(scores, loadings, dims, main, point_cex,
                                label_cex, dots) {
  s <- scores[, dims, drop = FALSE]
  l <- loadings[, dims, drop = FALSE]
  limits <- apply(s, 2L, range, finite = TRUE)
  args <- c(list(x = s[, 1L], y = s[, 2L], pch = 16, col = grDevices::adjustcolor("grey35", .45),
                 cex = point_cex, xlab = colnames(s)[1L], ylab = colnames(s)[2L],
                 main = main), dots)
  do.call(graphics::plot, args)
  factor <- min(diff(limits[, 1L]), diff(limits[, 2L])) /
    max(2 * max(abs(l)), .Machine$double.eps)
  graphics::arrows(0, 0, l[, 1L] * factor, l[, 2L] * factor,
                   length = 0.08, col = "darkred")
  graphics::text(l[, 1L] * factor, l[, 2L] * factor,
                 labels = rownames(l), cex = label_cex, pos = 3)
}

.nda_project_3d <- function(z, angle) {
  theta <- angle * pi / 180
  cbind(z[, 1L] + cos(theta) * z[, 3L],
        z[, 2L] + sin(theta) * z[, 3L])
}

#' @export
print.nda_biplot <- function(x, ...) {
  if (x$one_dimensional) {
    graphics::hist(x$scores[, 1L], probability = TRUE, col = "lightblue",
                   main = x$main %||% colnames(x$scores)[1L],
                   xlab = colnames(x$scores)[1L])
    if (nrow(x$scores) > 1L) graphics::lines(stats::density(x$scores[, 1L]),
                                              col = "darkred", lwd = 2)
    return(invisible(x))
  }
  if (!x$view_3D) {
    .nda_draw_biplot_2d(x$scores, x$loadings, x$dims, x$main,
                        x$point_cex, x$label_cex, c(x$dots, list(...)))
    return(invisible(x))
  }
  for (triple in x$triples) {
    s <- x$scores[, triple, drop = FALSE]
    l <- x$loadings[, triple, drop = FALSE]
    title <- x$main %||% paste(colnames(s), collapse = " / ")
    if (x$interactive) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Package 'plotly' is required for interactive 3D biplots.", call. = FALSE)
      }
      widget <- plotly::plot_ly(x = s[, 1L], y = s[, 2L], z = s[, 3L],
        type = "scatter3d", mode = "markers", marker = list(size = 3),
        name = "observations")
      widget <- plotly::add_trace(widget, x = l[, 1L], y = l[, 2L], z = l[, 3L],
        type = "scatter3d", mode = "markers+text", text = rownames(l),
        textposition = "top center", marker = list(color = "darkred", size = 4),
        name = "loadings")
      print(plotly::layout(widget, title = title))
    } else {
      projected <- .nda_project_3d(s, x$angle)
      projected_loadings <- .nda_project_3d(l, x$angle)
      graphics::plot(projected, pch = 16,
        col = grDevices::adjustcolor("grey35", .45), cex = x$point_cex,
        xlab = paste0(colnames(s)[1L], " + projected ", colnames(s)[3L]),
        ylab = paste0(colnames(s)[2L], " + projected ", colnames(s)[3L]),
        main = title)
      factor <- max(abs(projected)) / max(2 * max(abs(projected_loadings)),
                                          .Machine$double.eps)
      graphics::arrows(0, 0, projected_loadings[, 1L] * factor,
                       projected_loadings[, 2L] * factor, length = 0.08,
                       col = "darkred")
      graphics::text(projected_loadings[, 1L] * factor,
                     projected_loadings[, 2L] * factor,
                     labels = rownames(l), cex = x$label_cex, pos = 3)
    }
  }
  invisible(x)
}
