test_that("static NDA plots use the readable visual grammar", {
  association <- matrix(c(
    1, .8, -.6,
    .8, 1, -.4,
    -.6, -.4, 1
  ), 3, 3, byrow = TRUE,
  dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
  fit <- structure(list(
    R = association^2,
    correlation = association,
    covar = FALSE,
    membership = c(a = 1L, b = 1L, c = 2L)
  ), class = c("nda", "list"))

  result <- plot(fit, cuts = 0, interactive = FALSE, show_weights = TRUE)

  expect_s3_class(result, "nda_plot")
  expect_s3_class(result$static_grob, "gTree")
  expect_true(all(result$nodes$shape == "rectangle"))
  expect_true(all(result$edges$lty == 1))
  expect_true("#009E73" %in% result$edges$color)
  expect_true("#C62828" %in% result$edges$color)
})

test_that("NDRLM reproduces the layered NDR_PAPER path structure", {
  n <- 24L
  input <- data.frame(
    x1 = seq_len(n), x2 = seq_len(n) / 2,
    x3 = sin(seq_len(n)), x4 = cos(seq_len(n))
  )
  output <- data.frame(
    y1 = seq_len(n) / 3, y2 = seq_len(n) / 4,
    y3 = sin(seq_len(n) / 2), y4 = cos(seq_len(n) / 2)
  )
  latent_data <- data.frame(
    NDAin1 = scale(rowMeans(input))[, 1L]
  )
  latent_data$NDAout1 <- 0.8 * latent_data$NDAin1 +
    rep(c(-0.03, 0.01, 0.02), length.out = n)
  path_fit <- stats::lm(NDAout1 ~ NDAin1, data = latent_data)
  nda_in <- list(
    factors = 1L,
    membership = stats::setNames(rep(1L, 4L), names(input)),
    loadings = matrix(c(.90, .82, .76, .70), ncol = 1L,
      dimnames = list(names(input), "NDAin1"))
  )
  nda_out <- list(
    factors = 1L,
    membership = stats::setNames(rep(1L, 4L), names(output)),
    loadings = matrix(c(.88, .80, .74, .68), ncol = 1L,
      dimnames = list(names(output), "NDAout1"))
  )
  model <- structure(list(
    X = input, Y = output, latents = "both",
    NDAin = nda_in, NDAout = nda_out,
    fits = list(NDAout1 = path_fit)
  ), class = c("ndrlm", "list"))

  result <- plot(model, sig = 1, interactive = FALSE,
                 show_edge_labels = TRUE)
  node_y <- stats::setNames(result$layout[, 2L], result$nodes$id)
  node_x <- stats::setNames(result$layout[, 1L], result$nodes$id)
  path_edges <- result$edges[result$edges$role == "path", ]
  expected_standardized <- unname(stats::coef(path_fit)["NDAin1"]) *
    stats::sd(latent_data$NDAin1) / stats::sd(latent_data$NDAout1)

  expect_s3_class(result, "ndrlm_plot")
  expect_s3_class(result$static_grob, "gTree")
  expect_equal(unname(node_y[c("Lin:NDAin1", "Lout:NDAout1")]),
               c(0.5, 0.5))
  expect_true(all(result$nodes$shape[
    grepl("indicator", result$nodes$role)] == "rectangle"))
  expect_true(all(result$nodes$shape[
    grepl("latent", result$nodes$role)] == "circle"))
  expect_true(all(result$edges$lty[result$edges$role == "loading"] == 2))
  expect_true(all(result$edges$lty[result$edges$role == "path"] == 1))
  expect_equal(path_edges$weight, expected_standardized, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(
    path_edges$weight, unname(stats::coef(path_fit)["NDAin1"]))))
  expect_equal(unname(node_x[c("X:x1", "Lin:NDAin1",
                               "Lout:NDAout1", "Y:y1")]),
               c(0.20, 0.40, 0.60, 0.80))
  expect_true(all(substr(result$nodes$color, 8L, 9L) == "E6"))
  raw_result <- plot(model, sig = 1, interactive = FALSE,
                     standardized.paths = FALSE)
  raw_path <- raw_result$edges[raw_result$edges$role == "path", ]
  expect_equal(raw_path$weight,
               unname(stats::coef(path_fit)["NDAin1"]),
               tolerance = 1e-10)

  direct <- withVisible(plot(model, sig = 1, interactive = FALSE))
  assigned <- withVisible(saved <- plot(model, sig = 1,
                                        interactive = FALSE))
  expect_true(direct$visible)
  expect_false(assigned$visible)
  expect_s3_class(saved, "ndrlm_plot")

  legacy <- model
  legacy$fits <- unname(legacy$fits)
  legacy$dep_names <- NULL
  legacy$dep <- NULL
  legacy_result <- plot(legacy, sig = 1, interactive = FALSE)
  legacy_paths <- legacy_result$edges[legacy_result$edges$role == "path", ]
  expect_equal(legacy_paths$from, "Lin:NDAin1")
  expect_equal(legacy_paths$to, "Lout:NDAout1")
  expect_equal(legacy_paths$weight, expected_standardized,
               tolerance = 1e-10)
  expect_false(any(legacy_paths$to %in% paste0("Y:", names(output))))
})

test_that("NDRLM standardizes paths and closes gaps when a layer is absent", {
  x <- data.frame(x1 = 1:20, x2 = (1:20)^2)
  y <- data.frame(y1 = 5 + 3 * x$x1 - 0.04 * x$x2 +
                    rep(c(-1, 0, 1, 0), 5))
  fit <- ndrlm(y, x, latents = "none", optimize = FALSE,
               regression = "lm")
  result <- plot(fit, sig = 1, interactive = FALSE)
  expected <- stats::coef(fit$fits[[1]])[-1L] *
    vapply(x, stats::sd, numeric(1L)) / stats::sd(y$y1)
  paths <- result$edges[result$edges$role == "path", ]
  observed <- stats::setNames(paths$weight,
                              sub("^X:", "", paths$from))
  expect_equal(unname(observed[names(expected)]), unname(expected),
               tolerance = 1e-10)
  active_x <- sort(unique(result$layout[, 1L]))
  expect_equal(active_x, c(0.36, 0.64))
})

test_that("ridge and elastic-net plots use stored standardized coefficients", {
  skip_if_not_installed("glmnet")
  x <- data.frame(
    x1 = seq(-2, 2, length.out = 40),
    x2 = seq(-2, 2, length.out = 40)^2
  )
  y <- data.frame(y1 = 10 + 4 * x$x1 - 0.5 * x$x2)
  for (engine in c("ridge", "elastic_net")) {
    fit <- ndrlm(y, x, latents = "none", optimize = FALSE,
                 regression = engine, alpha = 0.5, lambda = 0.05)
    result <- plot(fit, sig = 1, interactive = FALSE)
    paths <- result$edges[result$edges$role == "path", ]
    observed <- stats::setNames(paths$weight,
                                sub("^X:", "", paths$from))
    expected <- fit$fits[[1]]$standardized_coefficients
    retained <- names(expected)[is.finite(expected) & expected != 0]
    expect_equal(unname(observed[retained]), unname(expected[retained]),
                 tolerance = 1e-10)
  }
})


