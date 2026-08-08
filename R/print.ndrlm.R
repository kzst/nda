#' Print an NDRLM model
#' @export
print.ndrlm <- function(x, digits = getOption("digits"), ...) {
  print(summary.ndrlm(x), digits = digits, ...)
  invisible(x)
}
