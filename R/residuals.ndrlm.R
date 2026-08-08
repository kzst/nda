#' Residuals from an NDRLM model
#' @export
residuals.ndrlm <- function(object, ...) {
  if (!inherits(object, "ndrlm")) stop("object must inherit from 'ndrlm'.", call. = FALSE)
  out <- vapply(object$fits, function(fit) fit$residuals, numeric(nrow(object$dep)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  dimnames(out) <- list(rownames(object$dep), names(object$fits))
  as.data.frame(out)
}
