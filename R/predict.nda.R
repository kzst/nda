#' Predict NDA latent-variable scores
#' @export
predict.nda <- function(object, newdata, ...) {
  if (!inherits(object, "nda")) stop("object must inherit from 'nda'.", call. = FALSE)
  if (missing(newdata) || is.null(newdata)) return(object$scores)
  if (is.null(object$scores)) {
    stop("Scores cannot be predicted from a model fitted with covar = TRUE.", call. = FALSE)
  }
  x <- as.matrix(.nda_numeric_matrix(newdata, "newdata"))
  required <- names(object$weight)
  missing_names <- setdiff(required, colnames(x))
  if (length(missing_names)) {
    stop("newdata is missing: ", paste(missing_names, collapse = ", "), call. = FALSE)
  }
  x <- x[, required, drop = FALSE]
  scaled <- .nda_scale(x, object$feature_center, object$feature_scale,
                       object$standardized)
  x <- sweep(scaled$x, 2L, object$weight, "*")
  raw <- matrix(0, nrow(x), object$factors)
  for (g in seq_len(object$factors)) {
    index <- match(names(object$EVCs[[g]]), colnames(x))
    raw[, g] <- as.numeric(x[, index, drop = FALSE] %*% object$EVCs[[g]])
  }
  scores <- sweep(sweep(raw, 2L, object$center, "-"), 2L, object$scale, "/")
  scores <- scores %*% object$rotation_matrix
  dimnames(scores) <- list(rownames(newdata), paste0("NDA", seq_len(object$factors)))
  scores
}
