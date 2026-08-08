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


#' Min-max normalization
#' @export
normalize <- function(x, type = c("all", "row", "col"), range = c(0, 1),
                      na.rm = TRUE) {
  type <- match.arg(type)
  if (!is.numeric(range) || length(range) != 2L || any(!is.finite(range)) ||
      range[1L] >= range[2L]) {
    stop("range must contain two increasing finite numbers.", call. = FALSE)
  }
  input_df <- is.data.frame(x)
  m <- as.matrix(.nda_numeric_matrix(x))
  rescale <- function(z) {
    limits <- base::range(z, na.rm = na.rm)
    span <- diff(limits)
    if (!is.finite(span) || span == 0) return(rep(mean(range), length(z)))
    range[1L] + (z - limits[1L]) * diff(range) / span
  }
  out <- switch(type,
    all = matrix(rescale(as.vector(m)), nrow(m), ncol(m), dimnames = dimnames(m)),
    row = t(apply(m, 1L, rescale)),
    col = apply(m, 2L, rescale))
  dimnames(out) <- dimnames(m)
  if (input_df) as.data.frame(out) else out
}
