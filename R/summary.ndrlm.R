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


#' Summarize an NDRLM model
#' @export
summary.ndrlm <- function(object, digits = getOption("digits"), ...) {
  if (!inherits(object, "ndrlm")) stop("object must inherit from 'ndrlm'.", call. = FALSE)
  out <- list(Call = object$Call, target = object$target, fval = object$fval,
              fit_stats = object$fit_stats, regression = object$regression,
              optimized = object$optimized, pareto = object$pareto,
              latents = object$latents, independent = object$indep_names,
              dependent = object$dep_names,
              NDAin = object$NDAin %||% NULL, NDAout = object$NDAout %||% NULL,
              fits = object$fits, extra_vars.X = object$extra_vars.X,
              extra_vars.Y = object$extra_vars.Y,
              dircon_X = object$dircon_X, dircon_Y = object$dircon_Y)
  class(out) <- c("summary.ndrlm", "list")
  out
}

#' @export
print.summary.ndrlm <- function(x, digits = getOption("digits"), ...) {
  cat("Network-based dimensionality reduction and regression\n")
  cat("Call:\n")
  print(x$Call)
  cat("\nRegression:", x$regression, "| Latent mode:", x$latents,
      "| Optimized:", if (x$optimized) "yes" else "no", "\n")
  cat("Target:", x$target, "| Objective:",
      paste(formatC(x$fval, digits = digits, format = "fg"), collapse = ", "), "\n")
  cat("Independent variables:", paste(x$independent, collapse = ", "), "\n")
  cat("Dependent variables:", paste(x$dependent, collapse = ", "), "\n\n")
  metric_table <- x$fit_stats
  metric_table$metric <- round(metric_table$metric, digits)
  print(metric_table, row.names = FALSE)
  for (name in names(x$fits)) {
    fit <- x$fits[[name]]
    cat("\nFit for", name, "\n")
    if (fit$method == "lm") {
      print(summary(fit$fit), digits = digits)
    } else {
      table <- data.frame(coefficient = names(fit$coefficients),
                          estimate = as.numeric(fit$coefficients), row.names = NULL)
      table$estimate <- round(table$estimate, digits)
      print(table, row.names = FALSE)
      if (!is.null(fit$lambda)) cat("Selected lambda:", format(fit$lambda), "\n")
    }
  }
  invisible(x)
}
