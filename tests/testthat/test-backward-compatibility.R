test_that("plot.nda accepts legacy nda objects without correlation", {
  association <- matrix(c(
    1.0, 0.8, 0.1,
    0.8, 1.0, 0.6,
    0.1, 0.6, 1.0
  ), 3, 3, byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  legacy <- structure(list(
    R = association,
    covar = FALSE,
    membership = c(1, 1, 0),
    loadings = matrix(c(0.8, 0.7, 0), 3, 1,
                      dimnames = list(c("a", "b", "c"), "NDA1")),
    factors = 1L
  ), class = c("nda", "list"))

  result <- plot(legacy, interactive = FALSE)
  expect_s3_class(result, "nda_plot")
  expect_equal(igraph::vcount(result$graph), 3L)
  expect_gt(igraph::ecount(result$graph), 0L)
  expect_length(unique(igraph::V(result$graph)$shape), 1L)
  expect_true(all(igraph::E(result$graph)$lty == 1))
  expect_true(all(igraph::E(result$graph)$width >= 1 &
                  igraph::E(result$graph)$width <= 6))
})

test_that("plot.nda separates community by colour and supports explicit shapes", {
  similarity <- matrix(c(0, .81, .49, .81, 0, .36, .49, .36, 0), 3, 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  association <- matrix(c(1, .9, -.7, .9, 1, -.6, -.7, -.6, 1), 3, 3,
    dimnames = dimnames(similarity))
  fit <- structure(list(R = similarity, correlation = association,
    covar = FALSE, membership = c(1, 1, 2)), class = c("nda", "list"))

  automatic <- plot(fit, cuts = 0, interactive = FALSE)
  expect_length(unique(igraph::V(automatic$graph)$shape), 1L)
  expect_true(all(igraph::V(automatic$graph)$shape == "rectangle"))
  expect_true(all(igraph::E(automatic$graph)$lty == 1))
  expect_true("#009E73" %in% igraph::E(automatic$graph)$color)
  expect_true("#C62828" %in% igraph::E(automatic$graph)$color)

  custom <- plot(fit, cuts = 0, interactive = FALSE,
                 shape = c("circle", "circle", "square"))
  expect_equal(igraph::V(custom$graph)$shape,
               c("circle", "circle", "square"))
})

test_that("plot.ndrlm accepts legacy lm fits", {
  data <- data.frame(
    x1 = seq_len(20),
    x2 = rep(c(0, 1), 10)
  )
  data$y <- 2 + 0.5 * data$x1 - data$x2 +
    rep(c(-0.15, 0.05, 0.10, -0.05), 5L)
  fit <- stats::lm(y ~ x1 + x2, data = data)
  legacy <- structure(list(
    X = data[, c("x1", "x2")],
    Y = data["y"],
    latents = "none",
    fits = list(y = fit)
  ), class = c("ndrlm", "list"))

  result <- plot(legacy, sig = 1, interactive = FALSE)
  expect_s3_class(result, "ndrlm_plot")
  expect_equal(igraph::vcount(result$graph), 3L)
  expect_equal(igraph::ecount(result$graph), 2L)
  expect_true(all(igraph::E(result$graph)$lty == 1))
  expect_true(all(igraph::E(result$graph)$width >= 1 &
                  igraph::E(result$graph)$width <= 6))
})

test_that("optimized ndrlm records an NSGA-II evaluation trace", {
  set.seed(42)
  X <- data.frame(
    x1 = stats::rnorm(24), x2 = stats::rnorm(24),
    x3 = stats::rnorm(24), x4 = stats::rnorm(24)
  )
  Y <- data.frame(y = X$x1 - X$x3 + stats::rnorm(24, sd = 0.2))
  model <- ndrlm(Y, X, latents = "in", optimize = TRUE,
    membership.X = c(1, 1, 2, 2), popsize = 4, generations = 1,
    lower.bounds.x = rep(0.5, 4), upper.bounds.x = rep(1.5, 4),
    lower.bounds.latentx = rep(0, 4),
    upper.bounds.latentx = rep(1e-8, 4),
    seed = 42)

  expect_s3_class(model, "ndrlm")
  expect_true(is.data.frame(model$NSGA_trace))
  expect_true(all(c("restart", "evaluation", "generation", "target_value",
                    "cost", "best_target") %in% names(model$NSGA_trace)))
  expect_gt(nrow(model$NSGA_trace), 0L)

  graph <- plot(model, sig = 1, interactive = FALSE)$graph
  node_color <- stats::setNames(igraph::V(graph)$color, igraph::V(graph)$name)
  expect_identical(node_color[["X:x1"]], node_color[["Lin:NDAin1"]])
  expect_identical(node_color[["X:x2"]], node_color[["Lin:NDAin1"]])
  expect_identical(node_color[["X:x3"]], node_color[["Lin:NDAin2"]])
  expect_equal(unname(igraph::V(graph)$shape[igraph::V(graph)$name == "X:x1"]),
               "rectangle")
  expect_equal(unname(igraph::V(graph)$shape[
    igraph::V(graph)$name == "Lin:NDAin1"]), "circle")
})

test_that("ndrlm exposes lm fits while retaining wrapped-fit compatibility", {
  data <- mtcars[seq_len(20), c("mpg", "disp", "hp", "wt")]
  model <- ndrlm(data["mpg"], data[c("disp", "hp", "wt")],
                 latents = "none")

  expect_s3_class(model$fits[[1L]], "lm")
  expect_s3_class(summary(model$fits[[1L]]), "summary.lm")
  expect_silent(stats::AIC(model$fits[[1L]]))
  expect_silent(stats::BIC(model$fits[[1L]]))
  if (requireNamespace("lm.beta", quietly = TRUE)) {
    expect_silent(summary(lm.beta::lm.beta(model$fits[[1L]])))
  }
  expect_equal(nrow(fitted(model)), nrow(data))
  expect_equal(nrow(residuals(model)), nrow(data))
  expect_equal(nrow(predict(model, data[c("disp", "hp", "wt")])), nrow(data))

  raw_fit <- model$fits[[1L]]
  metadata <- attr(raw_fit, "ndrlm_metadata", exact = TRUE)
  wrapped_fit <- c(list(method = "lm", fit = raw_fit), metadata)
  wrapped_model <- model
  wrapped_model$fits[[1L]] <- wrapped_fit

  expect_equal(fitted(wrapped_model), fitted(model))
  expect_equal(residuals(wrapped_model), residuals(model))
  expect_equal(predict(wrapped_model, data[c("disp", "hp", "wt")]),
               predict(model, data[c("disp", "hp", "wt")]))
  expect_output(print(summary(wrapped_model)),
                "Network-based dimensionality reduction and regression")
})

test_that("ndrlm keeps the new non-lm regression engines", {
  data <- mtcars[seq_len(20), c("mpg", "disp", "hp", "wt")]
  Y <- data["mpg"]
  X <- data[c("disp", "hp", "wt")]

  X_loess <- X["wt"]
  loess_model <- ndrlm(Y, X_loess, latents = "none",
                       regression = "loess", loess_span = 1)
  expect_identical(loess_model$fits[[1L]]$method, "loess")
  expect_false(inherits(loess_model$fits[[1L]], "lm"))
  expect_equal(nrow(predict(loess_model, X_loess)), nrow(data))
  expect_equal(nrow(fitted(loess_model)), nrow(data))
  expect_equal(nrow(residuals(loess_model)), nrow(data))

  skip_if_not_installed("glmnet")
  for (engine in c("ridge", "elastic_net")) {
    model <- ndrlm(Y, X, latents = "none", regression = engine,
                   lambda = 0.1)
    expect_identical(model$fits[[1L]]$method, engine)
    expect_false(inherits(model$fits[[1L]], "lm"))
    expect_equal(nrow(predict(model, X)), nrow(data))
    expect_equal(nrow(fitted(model)), nrow(data))
    expect_equal(nrow(residuals(model)), nrow(data))
    expect_output(print(summary(model)),
                  "Network-based dimensionality reduction and regression")
  }
})

test_that("fs.dimred refits after each removal", {
  set.seed(42)
  n <- 600
  f1 <- stats::rnorm(n)
  f2 <- stats::rnorm(n)
  data <- data.frame(
    a = 0.90 * f1 + stats::rnorm(n, sd = 0.20),
    b = 0.80 * f1 + stats::rnorm(n, sd = 0.25),
    c = 0.90 * f2 + stats::rnorm(n, sd = 0.20),
    d = 0.80 * f2 + stats::rnorm(n, sd = 0.25),
    cross = 0.65 * f1 + 0.65 * f2 + stats::rnorm(n, sd = 0.25),
    weak = stats::rnorm(n)
  )
  model <- psych::fa(data, nfactors = 2, rotate = "varimax", fm = "ml")
  result <- fs.dimred(model, data)

  expect_true(inherits(result, class(model)[1L]))
  expect_true(ncol(result$retained_DF) < ncol(data))
  removed <- c(names(result$dropped_low), names(result$dropped_com))
  expect_length(intersect(removed, names(result$retained_DF)), 0L)
  expect_equal(sort(c(removed, names(result$retained_DF))), sort(names(data)))
})
