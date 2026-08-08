#' Print an NDA model
#' @export
print.nda <- function(x, digits = getOption("digits"), ...) {
  print(summary.nda(x), digits = digits, ...)
  invisible(x)
}
