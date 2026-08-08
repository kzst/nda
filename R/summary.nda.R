#' Summarize an NDA model
#' @export
summary.nda <- function(object, digits = getOption("digits"), ...) {
  if (!inherits(object, "nda")) stop("object must inherit from 'nda'.", call. = FALSE)
  out <- list(
    Call = object$Call,
    factors = object$factors,
    n.obs = object$n.obs,
    modularity = object$modularity,
    membership = object$membership,
    community_sizes = object$stats$community_sizes,
    loadings = object$loadings,
    communality = object$communality,
    uniqueness = object$uniqueness,
    stats = object$stats,
    settings = list(cor_method = object$cor_method, cor_type = object$cor_type,
                    centrality = object$centrality,
                    null_model = object$stats$null_model,
                    standardized = object$standardized,
                    sparsecalc = object$sparsecalc))
  class(out) <- c("summary.nda", "list")
  out
}

#' @export
print.summary.nda <- function(x, digits = getOption("digits"), ...) {
  cat("Generalized network-based dimensionality reduction\n")
  cat("Call:\n")
  print(x$Call)
  cat("\nFactors:", x$factors)
  if (is.finite(x$n.obs)) cat(" | Observations:", x$n.obs)
  cat(" | Modularity:", formatC(x$modularity, digits = digits, format = "fg"), "\n")
  cat("Association:", x$settings$cor_type, x$settings$cor_method,
      "| Centrality:", x$settings$centrality,
      "| Null model:", x$settings$null_model, "\n")
  cat("Community sizes:", paste(x$community_sizes, collapse = ", "), "\n\n")
  cat("Loadings:\n")
  print(round(x$loadings, digits))
  cat("\nCommunalities:\n")
  print(round(x$communality, digits))
  if (length(x$stats$dropped)) {
    cat("\nDropped indicators:", paste(x$stats$dropped, collapse = ", "), "\n")
  }
  invisible(x)
}
