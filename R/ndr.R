# Generalized network-based dimensionality reduction.

.nda_centrality <- function(adjacency, method = "eigenvector", directed = FALSE) {
  n <- nrow(adjacency)
  if (n == 1L) return(1)
  graph <- igraph::graph_from_adjacency_matrix(adjacency,
    mode = if (directed) "directed" else "undirected", weighted = TRUE,
    diag = FALSE)
  if (igraph::ecount(graph) == 0L) return(rep(1 / n, n))
  value <- switch(method,
    eigenvector = igraph::eigen_centrality(graph, directed = directed,
                                           weights = igraph::E(graph)$weight)$vector,
    authority = igraph::authority_score(graph, weights = igraph::E(graph)$weight)$vector,
    hub = igraph::hub_score(graph, weights = igraph::E(graph)$weight)$vector,
    pagerank = igraph::page_rank(graph, directed = directed,
                                 weights = igraph::E(graph)$weight)$vector)
  value <- abs(as.numeric(value))
  value[!is.finite(value)] <- 0
  if (sum(value) <= .Machine$double.eps) value[] <- 1
  value / sum(value)
}

.nda_detect_communities <- function(adjacency, mod_mode, directed, Gamma) {
  graph <- igraph::graph_from_adjacency_matrix(adjacency,
    mode = if (directed) "directed" else "undirected", weighted = TRUE,
    diag = FALSE)
  if (igraph::ecount(graph) == 0L) return(seq_len(igraph::vcount(graph)))
  method <- .nda_choice(mod_mode,
    c("louvain", "fast_greedy", "leading_eigen", "infomap", "walktrap", "leiden"),
    c(`1` = "louvain", `2` = "fast_greedy", `3` = "leading_eigen",
      `4` = "infomap", `5` = "walktrap", `6` = "leiden"), "mod_mode")
  undirected <- if (directed) igraph::as.undirected(graph, mode = "collapse",
                                                    edge.attr.comb = "sum") else graph
  clustering <- switch(method,
    louvain = tryCatch(igraph::cluster_louvain(undirected,
      weights = igraph::E(undirected)$weight, resolution = Gamma),
      error = function(e) igraph::cluster_louvain(undirected,
        weights = igraph::E(undirected)$weight)),
    fast_greedy = igraph::cluster_fast_greedy(undirected,
      weights = igraph::E(undirected)$weight),
    leading_eigen = igraph::cluster_leading_eigen(undirected,
      weights = igraph::E(undirected)$weight),
    infomap = igraph::cluster_infomap(graph, e.weights = igraph::E(graph)$weight),
    walktrap = igraph::cluster_walktrap(undirected,
      weights = igraph::E(undirected)$weight),
    leiden = tryCatch(igraph::cluster_leiden(undirected,
      objective_function = "modularity", weights = igraph::E(undirected)$weight,
      resolution = Gamma), error = function(e) igraph::cluster_louvain(undirected,
        weights = igraph::E(undirected)$weight)))
  as.integer(igraph::membership(clustering))
}

.nda_renumber <- function(membership) {
  positive <- sort(unique(membership[membership > 0]))
  out <- integer(length(membership))
  for (i in seq_along(positive)) out[membership == positive[i]] <- i
  out
}

.nda_latent_state <- function(x_score, x_loading, association, similarity,
                              membership, centrality, directed, covar,
                              use_rotation, rotation) {
  groups <- sort(unique(membership[membership > 0]))
  k <- length(groups)
  p <- length(membership)
  weights <- vector("list", k)
  names(weights) <- paste0("NDA", seq_len(k))
  raw_scores <- if (!covar) matrix(0, nrow(x_score), k) else NULL
  loadings <- matrix(0, p, k,
    dimnames = list(rownames(association), paste0("NDA", seq_len(k))))
  for (g in seq_along(groups)) {
    index <- which(membership == groups[g])
    w <- .nda_centrality(similarity[index, index, drop = FALSE], centrality, directed)
    names(w) <- rownames(association)[index]
    weights[[g]] <- w
    if (!covar) {
      raw_scores[, g] <- as.numeric(x_score[, index, drop = FALSE] %*% w)
    } else {
      latent_var <- as.numeric(crossprod(w, association[index, index, drop = FALSE] %*% w))
      if (!is.finite(latent_var) || latent_var <= 0) latent_var <- 1
      loadings[, g] <- as.numeric(association[, index, drop = FALSE] %*% w) /
        sqrt(pmax(diag(association), .Machine$double.eps) * latent_var)
    }
  }
  rotation_matrix <- diag(k)
  latent_center <- rep(0, k)
  latent_scale <- rep(1, k)
  if (!covar) {
    scaled <- .nda_scale(raw_scores, enabled = TRUE)
    scores <- scaled$x
    latent_center <- scaled$center
    latent_scale <- scaled$scale
    if (k > 1L && use_rotation) {
      rotated <- psych::principal(scores, nfactors = k, rotate = rotation,
                                  scores = TRUE)$scores
      rotation_matrix <- tryCatch(qr.solve(scores, rotated),
                                  error = function(e) diag(k))
      scores <- scores %*% rotation_matrix
    }
    loadings <- stats::cor(x_loading, scores, use = "pairwise.complete.obs")
    loadings[!is.finite(loadings)] <- 0
    colnames(loadings) <- paste0("NDA", seq_len(k))
  } else if (k > 1L && use_rotation) {
    rotated <- tryCatch(stats::varimax(loadings)$loadings,
                        error = function(e) loadings)
    rotation_matrix <- tryCatch(qr.solve(loadings, as.matrix(rotated)),
                                error = function(e) diag(k))
    loadings <- as.matrix(rotated)
  } else {
    scores <- NULL
  }
  communality <- apply(loadings^2, 1L, max, na.rm = TRUE)
  communality[!is.finite(communality)] <- 0
  list(scores = if (covar) NULL else scores, loadings = loadings,
       communality = communality, uniqueness = pmax(0, 1 - communality),
       EVCs = weights, latent_center = latent_center,
       latent_scale = latent_scale, rotation_matrix = rotation_matrix)
}

.nda_indicator_statistics <- function(loadings, n, adjust = "BH") {
  if (is.null(n) || !is.finite(n) || n <= 2L) return(NULL)
  value <- pmin(abs(loadings), 1 - sqrt(.Machine$double.eps))
  t_value <- value * sqrt((n - 2) / pmax(1 - value^2, .Machine$double.eps))
  p <- 2 * stats::pt(-t_value, df = n - 2)
  p[] <- stats::p.adjust(as.vector(p), method = adjust)
  dimnames(p) <- dimnames(loadings)
  p
}

.nda_bootstrap_loadings <- function(x_score, x_loading, membership, evcs,
                                    latent_center, latent_scale, rotation_matrix,
                                    n_boot, conf, parallel, cores) {
  n <- nrow(x_score)
  p <- ncol(x_loading)
  k <- length(evcs)
  one <- function(i) {
    rows <- sample.int(n, n, replace = TRUE)
    scores <- matrix(0, n, k)
    for (g in seq_len(k)) {
      index <- match(names(evcs[[g]]), colnames(x_score))
      scores[, g] <- as.numeric(x_score[rows, index, drop = FALSE] %*% evcs[[g]])
    }
    scores <- sweep(sweep(scores, 2L, latent_center, "-"), 2L, latent_scale, "/")
    scores <- scores %*% rotation_matrix
    out <- stats::cor(x_loading[rows, , drop = FALSE], scores,
                      use = "pairwise.complete.obs")
    out[!is.finite(out)] <- 0
    out
  }
  jobs <- as.list(seq_len(n_boot))
  estimates <- if (parallel && cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(jobs, one, mc.cores = min(cores, n_boot))
  } else lapply(jobs, one)
  array_est <- array(unlist(estimates), dim = c(p, k, n_boot))
  alpha <- (1 - conf) / 2
  list(lower = apply(array_est, c(1, 2), stats::quantile, probs = alpha,
                     na.rm = TRUE),
       upper = apply(array_est, c(1, 2), stats::quantile, probs = 1 - alpha,
                     na.rm = TRUE),
       replicates = n_boot, confidence = conf)
}

#' Generalized network-based dimensionality reduction and analysis
#' @export
ndr <- function(r, covar = FALSE, cor_method = 1, cor_type = 1, min_R = 0,
                min_comm = 2, Gamma = 1, null_model_type = 4, mod_mode = 6,
                min_evalue = 0, min_communality = 0,
                com_communalities = 0, use_rotation = FALSE,
                rotation = "oblimin", weight = NULL, seed = NULL,
                parallel = FALSE, cores = 1L, sparsecalc = FALSE,
                membership = NULL,
                centrality = c("eigenvector", "authority", "hub", "pagerank"),
                standardized = TRUE, P_mat = NULL, test = FALSE,
                p_adjust = "BH", alpha = NULL, bootstrap = FALSE,
                n_boot = 200L, conf = 0.95) {
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  centrality <- match.arg(centrality)
  cor_method <- .nda_cor_method(cor_method)
  cor_type <- .nda_cor_type(cor_type)
  null_model_type <- as.integer(null_model_type)
  if (!null_model_type %in% 1:4) stop("null_model_type must be 1, 2, 3, or 4.", call. = FALSE)
  if (!is.numeric(Gamma) || length(Gamma) != 1L || !is.finite(Gamma) || Gamma <= 0) {
    stop("Gamma must be a positive finite number.", call. = FALSE)
  }
  if (!p_adjust %in% stats::p.adjust.methods) {
    stop("p_adjust must be one of stats::p.adjust.methods.", call. = FALSE)
  }
  if (!is.null(alpha) && (!is.numeric(alpha) || length(alpha) != 1L ||
                          !is.finite(alpha) || alpha <= 0 || alpha >= 1)) {
    stop("alpha must be a number between zero and one.", call. = FALSE)
  }
  if (!is.numeric(cores) || length(cores) != 1L || cores < 1) {
    stop("cores must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1) {
    stop("conf must be a number between zero and one.", call. = FALSE)
  }
  if (!is.numeric(min_comm) || min_comm < 1) stop("min_comm must be positive.", call. = FALSE)
  input <- .nda_numeric_matrix(r)
  p <- ncol(input)
  variable_names <- colnames(input)
  if (is.null(weight)) weight <- rep(1, p)
  if (!is.numeric(weight) || length(weight) != p || any(!is.finite(weight))) {
    stop("weight must contain one finite number per variable.", call. = FALSE)
  }
  names(weight) <- variable_names

  if (covar) {
    association <- as.matrix(input)
    if (nrow(association) != ncol(association)) {
      stop("With covar = TRUE, r must be a square association matrix.", call. = FALSE)
    }
    if (is.null(rownames(association))) rownames(association) <- variable_names
    if (!setequal(rownames(association), colnames(association))) {
      stop("Association-matrix row and column names must match.", call. = FALSE)
    }
    association <- association[colnames(association), colnames(association), drop = FALSE]
    association_p <- NULL
    x_score <- x_loading <- NULL
    feature_center <- rep(0, p)
    feature_scale <- rep(1, p)
    n_obs <- NA_integer_
  } else {
    dense <- as.matrix(input)
    scaled <- .nda_scale(dense, enabled = standardized)
    x_loading <- scaled$x
    x_score <- sweep(x_loading, 2L, weight, "*")
    feature_center <- scaled$center
    feature_scale <- scaled$scale
    n_obs <- nrow(dense)
    assoc_result <- .nda_association(x_score, cor_method, cor_type, test,
                                     p_adjust, alpha, parallel, cores)
    association <- if (test) assoc_result$r else assoc_result
    association_p <- if (test) assoc_result$p else NULL
  }
  association[!is.finite(association)] <- 0
  diag(association) <- 1
  similarity <- association^2
  diag(similarity) <- 0
  similarity[similarity < min_R] <- 0
  directed <- !isSymmetric(association, tol = sqrt(.Machine$double.eps))

  total <- sum(similarity)
  expected <- if (total > 0) outer(rowSums(similarity), colSums(similarity)) / total else
    matrix(0, p, p)
  if (!is.null(P_mat)) {
    P_mat <- as.matrix(P_mat)
    if (!all(dim(P_mat) == c(p, p)) || any(!is.finite(P_mat))) {
      stop("P_mat must be a finite matrix with the same dimensions as the association matrix.",
           call. = FALSE)
    }
    expected <- P_mat
    network <- pmax(similarity - Gamma * expected, 0)
    null_label <- "manual"
  } else {
    positive <- similarity[similarity > 0]
    network <- switch(as.character(null_model_type),
      `1` = pmax(similarity - Gamma * expected, 0),
      `2` = pmax(similarity - Gamma * if (length(positive)) mean(positive) else 0, 0),
      `3` = pmax(similarity - Gamma * min_R, 0),
      `4` = similarity)
    null_label <- c("configuration", "mean", "threshold", "none")[null_model_type]
  }
  diag(network) <- 0
  dimnames(network) <- list(variable_names, variable_names)

  if (is.null(membership)) {
    membership <- .nda_detect_communities(if (sparsecalc) Matrix::Matrix(network, sparse = TRUE) else network,
                                          mod_mode, directed, Gamma)
  } else {
    if (length(membership) != p || anyNA(membership)) {
      stop("membership must contain one non-missing value per variable.", call. = FALSE)
    }
    excluded <- membership == 0
    membership <- as.integer(factor(membership, levels = unique(membership)))
    membership[excluded] <- 0L
  }
  sizes <- table(membership)
  small <- as.integer(names(sizes)[sizes < min_comm])
  if (p > 1L && length(small)) membership[membership %in% small] <- 0L
  if (!any(membership > 0L)) {
    largest <- as.integer(names(which.max(sizes)))
    membership[membership == largest] <- 1L
  }
  membership <- .nda_renumber(membership)

  # Centrality-based periphery filtering is performed once, retaining a valid
  # community even when only one indicator remains.
  for (g in sort(unique(membership[membership > 0]))) {
    index <- which(membership == g)
    values <- .nda_centrality(similarity[index, index, drop = FALSE], centrality, directed)
    keep <- values > min_evalue
    minimum <- if (length(index) == 1L) 1L else min(as.integer(min_comm), length(index))
    if (sum(keep) >= minimum) membership[index[!keep]] <- 0L
  }
  membership <- .nda_renumber(membership)

  fit_state <- function(current) .nda_latent_state(
    x_score, x_loading, association, similarity, current, centrality, directed,
    covar, use_rotation, rotation)
  state <- fit_state(membership)

  # Iterative feature filtering only recomputes latent scores when a feature is
  # actually removed. This avoids the repeated graph reconstruction in v0.2.x.
  for (iteration in seq_len(max(1L, p))) {
    candidates <- which(membership > 0 & state$communality < min_communality)
    if (!length(candidates)) break
    removed <- integer()
    for (g in sort(unique(membership[candidates]))) {
      group <- which(membership == g)
      group_candidates <- group[order(state$communality[group])]
      can_remove <- max(0L, length(group) - min(as.integer(min_comm), length(group)))
      if (can_remove > 0L) {
        removed <- c(removed, utils::head(group_candidates, can_remove))
      }
    }
    if (!length(removed)) break
    membership[removed] <- 0L
    membership <- .nda_renumber(membership)
    state <- fit_state(membership)
  }

  if (com_communalities > 0 && ncol(state$loadings) > 1L) {
    for (iteration in seq_len(max(1L, p))) {
      squared <- state$loadings^2
      ordered <- t(apply(squared, 1L, sort, decreasing = TRUE))
      ambiguous <- which(membership > 0 & ordered[, 1L] < ordered[, 2L] + com_communalities &
                           ordered[, 1L] <= 2 * ordered[, 2L])
      if (!length(ambiguous)) break
      candidate <- ambiguous[which.min(state$communality[ambiguous])]
      group <- which(membership == membership[candidate])
      if (length(group) <= min(as.integer(min_comm), length(group))) break
      membership[candidate] <- 0L
      membership <- .nda_renumber(membership)
      state <- fit_state(membership)
    }
  }

  groups <- sort(unique(membership[membership > 0]))
  factor_names <- paste0("NDA", seq_along(groups))
  colnames(state$loadings) <- factor_names
  if (!is.null(state$scores)) {
    colnames(state$scores) <- factor_names
    rownames(state$scores) <- rownames(input)
  }
  names(state$EVCs) <- factor_names
  reconstruction <- if (!is.null(state$scores)) {
    state$scores %*% t(state$loadings)
  } else {
    state$loadings %*% t(state$loadings)
  }
  indicator_p <- .nda_indicator_statistics(state$loadings, n_obs, p_adjust)
  bootstrap_stats <- NULL
  if (bootstrap) {
    if (covar) {
      warning("Bootstrap statistics require row-level data and are unavailable with covar = TRUE.",
              call. = FALSE)
    } else {
      if (n_boot < 2L) stop("n_boot must be at least two.", call. = FALSE)
      bootstrap_stats <- .nda_bootstrap_loadings(
        x_score, x_loading, membership, state$EVCs, state$latent_center,
        state$latent_scale, state$rotation_matrix, as.integer(n_boot), conf,
        parallel, as.integer(cores))
      dimnames(bootstrap_stats$lower) <- dimnames(state$loadings)
      dimnames(bootstrap_stats$upper) <- dimnames(state$loadings)
    }
  }

  graph <- igraph::graph_from_adjacency_matrix(network,
    mode = if (directed) "directed" else "undirected", weighted = TRUE,
    diag = FALSE)
  active <- membership > 0L
  active_graph <- igraph::induced_subgraph(graph, which(active))
  modularity_value <- if (sum(active) < 2L || igraph::ecount(active_graph) == 0L) {
    0
  } else tryCatch(igraph::modularity(active_graph, membership[active],
    weights = igraph::E(active_graph)$weight, resolution = Gamma),
    error = function(e) tryCatch(igraph::modularity(active_graph, membership[active],
      weights = igraph::E(active_graph)$weight), error = function(e2) NA_real_))

  result <- list(
    communality = state$communality,
    loadings = state$loadings,
    uniqueness = state$uniqueness,
    factors = length(groups), scores = state$scores,
    reconstruction = reconstruction, n.obs = n_obs,
    R = if (sparsecalc) Matrix::Matrix(similarity, sparse = TRUE) else similarity,
    correlation = association, correlation_p = association_p,
    EVCs = state$EVCs, center = state$latent_center,
    scale = state$latent_scale, feature_center = feature_center,
    feature_scale = feature_scale, membership = membership, weight = weight,
    centrality = centrality, standardized = standardized, covar = covar,
    sparsecalc = sparsecalc, use_rotation = use_rotation, rotation = rotation,
    rotation_matrix = state$rotation_matrix, modularity = modularity_value,
    expected = expected, null_model_type = null_model_type,
    stats = list(indicator_p = indicator_p, bootstrap = bootstrap_stats,
                 community_sizes = table(membership[membership > 0]),
                 retained = variable_names[membership > 0],
                 dropped = variable_names[membership == 0],
                 null_model = null_label, directed = directed),
    cor_method = cor_method, cor_type = cor_type, min_R = min_R,
    min_comm = min_comm, Gamma = Gamma, mod_mode = mod_mode,
    fn = "NDA", seed = seed, Call = call)
  class(result) <- c("nda", "list")
  result
}
