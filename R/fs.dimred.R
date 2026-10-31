#' Feature selection from dimensionality-reduction loadings
#' @export
fs.dimred <- function(fn, DF, min_comm = 0.25, com_comm = 0.25, ...) {
  if (!is.numeric(min_comm) || length(min_comm) != 1L ||
      !is.finite(min_comm) || min_comm < 0) {
    stop("min_comm must be one non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(com_comm) || length(com_comm) != 1L ||
      !is.finite(com_comm) || com_comm < 0) {
    stop("com_comm must be one non-negative finite number.", call. = FALSE)
  }

  numeric_data <- as.matrix(.nda_numeric_matrix(DF, "DF", min_cols = 2L))
  correlation_input <- nrow(numeric_data) == ncol(numeric_data) &&
    !is.null(rownames(numeric_data)) && isSymmetric(numeric_data) &&
    setequal(rownames(numeric_data), colnames(numeric_data))
  working <- if (correlation_input) {
    numeric_data[colnames(numeric_data), colnames(numeric_data), drop = FALSE]
  } else if (is.data.frame(DF)) {
    as.data.frame(DF)
  } else {
    numeric_data
  }

  model_call <- fn$Call %||% fn$call
  can_refit <- !is.null(model_call) && is.call(model_call)

  extract_loadings <- function(model, variable_names) {
    value <- if (inherits(model, "nda")) {
      model$loadings
    } else if (inherits(model, "prcomp")) {
      model$rotation
    } else if (inherits(model, "princomp")) {
      unclass(model$loadings)
    } else if (!is.null(model$loadings)) {
      unclass(model$loadings)
    } else if (!is.null(model$rotation)) {
      model$rotation
    } else if (!is.null(model$var$coord)) {
      model$var$coord
    } else {
      stop("fn does not expose loadings, a rotation matrix, or variable coordinates.",
           call. = FALSE)
    }
    value <- as.matrix(value)
    if (nrow(value) != length(variable_names)) {
      stop("The number of loading rows does not match the variables in DF.",
           call. = FALSE)
    }
    if (is.null(rownames(value))) rownames(value) <- variable_names
    if (!setequal(rownames(value), variable_names)) {
      stop("The loading row names do not match the variables in DF.",
           call. = FALSE)
    }
    value[variable_names, , drop = FALSE]
  }

  extract_communality <- function(model, loadings) {
    value <- model$communality
    if (is.null(value)) value <- rowSums(loadings^2, na.rm = TRUE)
    value <- as.numeric(value)
    if (length(value) != nrow(loadings)) {
      stop("The communality vector does not match the loading matrix.",
           call. = FALSE)
    }
    names(value) <- rownames(loadings)
    value
  }

  refit <- function(data) {
    call <- model_call
    call_name <- paste(deparse(call[[1L]]), collapse = "")
    function_name <- sub("^.*::+", "", call_name)
    namespace <- if (function_name %in% c("fa", "principal")) {
      "psych"
    } else if (function_name %in% c("prcomp", "princomp")) {
      "stats"
    } else if (function_name == "ndr") {
      "nda"
    } else {
      stop("The fitted model cannot be refitted safely. Supported functions are ",
           "psych::fa, psych::principal, stats::prcomp, stats::princomp, and nda::ndr.",
           call. = FALSE)
    }
    call[[1L]] <- getExportedValue(namespace, function_name)
    argument_names <- names(call)
    data_position <- match(c("r", "x", "data"), argument_names, nomatch = 0L)
    data_position <- data_position[data_position > 0L][1L]
    if (!length(data_position) || is.na(data_position)) data_position <- 2L
    call[[data_position]] <- quote(.fs_dimred_data)
    tryCatch(
      eval(call, envir = list(.fs_dimred_data = data), enclos = parent.frame()),
      error = function(e) stop("Refitting the dimensionality-reduction model failed: ",
                               conditionMessage(e), call. = FALSE)
    )
  }

  drop_variable <- function(data, variable) {
    keep <- setdiff(colnames(data), variable)
    if (correlation_input) {
      data[keep, keep, drop = FALSE]
    } else {
      data[, keep, drop = FALSE]
    }
  }

  current_model <- fn
  dropped_low_names <- character()
  dropped_com_names <- character()

  if (can_refit) {
    repeat {
      loadings <- extract_loadings(current_model, colnames(working))
      communality <- extract_communality(current_model, loadings)
      candidates <- which(is.finite(communality) & communality < min_comm)
      if (!length(candidates) || ncol(working) <= 2L) break
      index <- candidates[which.min(communality[candidates])]
      variable <- names(communality)[index]
      dropped_low_names <- c(dropped_low_names, variable)
      working <- drop_variable(working, variable)
      current_model <- refit(working)
    }

    repeat {
      loadings <- extract_loadings(current_model, colnames(working))
      communality <- extract_communality(current_model, loadings)
      if (ncol(loadings) < 2L || ncol(working) <= 2L) break
      absolute_loadings <- abs(loadings)
      ordered <- t(apply(absolute_loadings, 1L, sort, decreasing = TRUE))
      candidates <- which(is.finite(ordered[, 1L]) &
                          is.finite(ordered[, 2L]) &
                          ordered[, 1L] < 2 * ordered[, 2L] &
                          ordered[, 1L] < ordered[, 2L] + com_comm)
      if (!length(candidates)) break
      index <- candidates[which.min(communality[candidates])]
      variable <- rownames(loadings)[index]
      dropped_com_names <- c(dropped_com_names, variable)
      working <- drop_variable(working, variable)
      current_model <- refit(working)
    }
  } else {
    loadings <- extract_loadings(current_model, colnames(working))
    communality <- extract_communality(current_model, loadings)
    dropped_low_names <- names(communality)[
      is.finite(communality) & communality < min_comm]
    candidates <- character()
    if (ncol(loadings) > 1L) {
      ordered <- t(apply(abs(loadings), 1L, sort, decreasing = TRUE))
      candidates <- rownames(loadings)[
        is.finite(ordered[, 1L]) & is.finite(ordered[, 2L]) &
        ordered[, 1L] < 2 * ordered[, 2L] &
        ordered[, 1L] < ordered[, 2L] + com_comm]
    }
    dropped_com_names <- setdiff(candidates, dropped_low_names)
    for (variable in c(dropped_low_names, dropped_com_names)) {
      if (ncol(working) <= 2L) break
      working <- drop_variable(working, variable)
    }
  }

  dropped_data <- function(variable_names) {
    if (correlation_input || !length(variable_names)) return(NULL)
    out <- DF[, variable_names, drop = FALSE]
    if (is.data.frame(DF)) as.data.frame(out) else out
  }

  current_model$dropped_low <- dropped_data(dropped_low_names)
  current_model$dropped_com <- dropped_data(dropped_com_names)
  current_model$retained_DF <- working
  current_model$remain_DF <- working
  current_model$model <- fn
  current_model
}
