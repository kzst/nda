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
})

test_that("plot.ndrlm accepts legacy lm fits", {
  data <- data.frame(
    x1 = seq_len(20),
    x2 = rep(c(0, 1), 10)
  )
  data$y <- 2 + 0.5 * data$x1 - data$x2
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
