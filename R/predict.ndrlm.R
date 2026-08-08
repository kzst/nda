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


#' Predict from an NDRLM model
#' @export
predict.ndrlm <- function(object, newdata = NULL, se.fit = FALSE,
                          interval = c("none", "confidence", "prediction"),
                          level = 0.95, ...) {
  if (!inherits(object, "ndrlm")) stop("object must inherit from 'ndrlm'.", call. = FALSE)
  interval <- match.arg(interval)
  if (is.null(newdata)) return(stats::fitted(object))
  data <- as.data.frame(as.matrix(.nda_numeric_matrix(newdata, "newdata")))
  missing_x <- setdiff(colnames(object$X), colnames(data))
  if (length(missing_x)) {
    stop("newdata is missing: ", paste(missing_x, collapse = ", "), call. = FALSE)
  }
  new_x <- data[, colnames(object$X), drop = FALSE]
  if (object$latents %in% c("in", "both")) {
    indep <- as.data.frame(stats::predict(object$NDAin, new_x))
    colnames(indep) <- paste0("NDAin", seq_len(object$NDAin$factors))
    if (object$extra_vars.X) indep <- cbind(indep, new_x[, object$dircon_X, drop = FALSE])
  } else indep <- new_x
  indep <- indep[, object$indep_names, drop = FALSE]
  out <- lapply(object$fits, .ndrlm_predict_one, newdata = indep,
                se.fit = se.fit, interval = interval, level = level, ...)
  names(out) <- names(object$fits)
  if (!se.fit && interval == "none" && all(vapply(out, is.atomic, logical(1)))) {
    matrix_out <- do.call(cbind, out)
    dimnames(matrix_out) <- list(rownames(newdata), names(out))
    return(as.data.frame(matrix_out))
  }
  out
}
