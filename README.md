# Generalized network-based dimensionality reduction and analysis

`nda` constructs interpretable latent variables from communities in
squared-association networks. It supports unsupervised dimensionality
reduction with `ndr()` and supervised modelling with `ndrlm()`.

Version 0.3.0 adds faster and safer association calculations, sparse and
parallel workflows, fixed memberships, inference and bootstrap diagnostics,
new centrality definitions, prediction methods, configurable graphics, and
ordinary least squares, ridge, elastic-net, and LOESS response engines.

## Installation

Install the development version from GitHub:

```r
remotes::install_github("kzst/nda")
library(nda)
```

## Main function 1: `ndr()`

`ndr()`:

1. estimates or accepts an association matrix;
2. constructs a squared-association network;
3. detects indicator communities;
4. creates one standardized, centrality-weighted latent score per retained
   community.

```r
fit <- ndr(
  swiss,
  cor_method = "spearman",
  cor_type = "full",
  centrality = "pagerank",
  min_comm = 2,
  seed = 42
)

fit
summary(fit)
round(fit$loadings, 3)

# Score new observations with the training transformation.
predict(fit, swiss[1:5, , drop = FALSE])
```

Network and biplot functions return printable objects. Assignment is silent;
print the object to render it.

```r
network <- plot(
  fit,
  interactive = FALSE,
  show_weights = TRUE,
  shape = "automatic"
)
network

# The same visual grammar as an embedded interactive HTML widget.
interactive_network <- plot(fit, interactive = TRUE, show_weights = TRUE)
interactive_network

score_plot <- biplot(fit, view_3D = FALSE, interactive = FALSE)
score_plot
```

Community membership is shown with a light qualitative node-colour palette.
Indicators use one common rectangular shape by default; positive and negative
associations use green and red edges, and absolute association controls edge
width. Shapes, colours, line types, and the minimum/maximum width remain
configurable.

`covar = TRUE` accepts an arbitrary square similarity matrix. For example,
NDR can cluster observations from normalized Euclidean similarities:

```r
z <- scale(swiss)
d <- as.matrix(dist(z))
similarity <- 1 - d / max(d)
cluster_fit <- ndr(similarity, covar = TRUE, seed = 42)
table(cluster_fit$membership)
```

Important options include:

- full, partial, and semipartial associations;
- Pearson, Spearman, Kendall, distance, Gaussian-rank, and biweight methods;
- eigenvector, authority, hub, and PageRank centrality;
- fixed memberships and custom expected-edge matrices;
- FDR-adjusted association screening;
- bootstrap loading intervals;
- sparse networks and multicore distance calculations;
- prediction using stored training transformations.

## Main function 2: `ndrlm()`

`ndrlm()` connects original or NDA-reduced predictors to original or
NDA-reduced responses. Set `latents` to `"none"`, `"in"`, `"out"`, or
`"both"`.

The example below uses the two-community latent input structure from the
published low-sample-size `mtcars` application:

```r
X <- mtcars[, setdiff(names(mtcars), "mpg")]
Y <- mtcars["mpg"]
membership <- c(
  cyl = 1, disp = 2, hp = 1, drat = 2, wt = 2,
  qsec = 1, vs = 1, am = 2, gear = 2, carb = 1
)

model <- ndrlm(
  Y,
  X,
  latents = "in",
  optimize = FALSE,
  regression = "lm",
  membership.X = membership,
  seed = 42
)

summary(model)
predict(model, X[1:5, , drop = FALSE])
head(fitted(model))
head(residuals(model))
```

For `regression = "lm"`, each `model$fits[[i]]` is a genuine `lm` object and
can be used directly with `summary()`, `AIC()`, `BIC()`, and
`lm.beta::lm.beta()`.

```r
inherits(model$fits[[1]], "lm")
summary(model$fits[[1]])
```

Ridge and elastic net use `glmnet` and retain NDRLM engine wrappers:

```r
ridge <- ndrlm(
  Y,
  X,
  latents = "none",
  optimize = FALSE,
  regression = "ridge",
  lambda = 0.1
)
predict(ridge, X[1:5, , drop = FALSE])
```

The package methods also accept legacy direct `lm` fits and early 0.3.0
wrapped `lm` fits. For engine-independent code, use `summary()`, `predict()`,
`fitted()`, `residuals()`, and `plot()` on the enclosing `ndrlm` object.

```r
static_paths <- plot(
  model,
  interactive = FALSE,
  shape = "automatic"
)
static_paths

interactive_paths <- plot(
  model,
  interactive = TRUE,
  show_edge_labels = TRUE
)
interactive_paths
```

Observed variables are rectangles and latent variables are circles. Each
latent variable and its indicators share one light community colour. Positive
paths are green, negative paths are red, and absolute standardized coefficient
or loading magnitude controls width within its edge role. The plot computes
standardized structural coefficients for ordinary linear, ridge, and
elastic-net regression from the complete training observations; use
`standardized.paths = FALSE` only when raw coefficients are required.
Measurement links between
observed and latent variables are dashed; structural paths are solid. Latent
nodes are centred on their indicator groups, absent layers leave no empty
horizontal column, and the default node fills are slightly translucent and
lightened. Edge labels use a light background offset from the path to avoid
overwriting thick lines. For legacy
unnamed fit lists, the path target is recovered from the fitted response
formula before any original-variable fallback is considered; standardized
coefficients are reconstructed from the legacy fitted model.

Static plots use the package's layered grid renderer by default. The legacy
igraph renderer remains available with `static.engine = "igraph"`. Both static
and interactive variants are printable objects, so the following pattern
works with screen and file devices:

```r
p <- plot(model, interactive = FALSE)  # no drawing during assignment
p                                      # draw now

pdf("ndrlm-paths.pdf", width = 10, height = 6)
p
dev.off()
```

Expensive NSGA-II optimization is optional:

```r
optimized <- ndrlm(
  Y,
  X,
  latents = "in",
  regression = "lm",
  target = "adj.r.square",
  popsize = 20,
  generations = 50,
  restarts = 3,
  seed = 42
)

plot(
  optimized$NSGA_trace$evaluation,
  optimized$NSGA_trace$best_target,
  type = "s",
  xlab = "Objective evaluation",
  ylab = "Best target value"
)
```

## Supporting functions

| Purpose | Functions |
|---|---|
| Structured example data | `data_gen()` |
| Global, row, or column scaling | `normalize()` |
| Distance covariance and correlation | `dCov()`, `dCor()` |
| Partial and semipartial distance correlation | `pdCor()`, `spdCor()` |
| KMO-based feature selection | `fs.KMO()` |
| Loading-based feature selection and iterative refitting | `fs.dimred()` |

## Vignettes

- A visual tour of the nda package:
  [source](vignettes/nda-function-tour.Rmd) |
  [HTML preview](vignette-preview/nda-function-tour.html)
- Generalized NDR workflow:
  [source](vignettes/ndr-workflow.Rmd) |
  [HTML preview](vignette-preview/ndr-workflow.html)
- Supervised NDRLM workflow:
  [source](vignettes/ndrlm-workflow.Rmd) |
  [HTML preview](vignette-preview/ndrlm-workflow.html)

Standalone HTML previews are stored in [`vignette-preview/`](vignette-preview/).
They are excluded from the source package build; `R CMD build` generates the
official installed vignettes from the R Markdown sources.

## Reproducibility and performance

- Use `seed` whenever community detection, bootstrap inference, or NSGA-II is
  involved.
- Use small `n_boot` and optimization settings while developing, then increase
  them for final analysis.
- Parallel workers mainly benefit distance-association and bootstrap tasks.
- `sparsecalc = TRUE` is useful after meaningful edge thresholding.
- Report association, community, centrality, filtering, optimization, and
  regression settings with fitted results.

## References

Kosztyan, Z. T., Kurbucz, M. T., and Katona, A. I. (2022).
Network-based dimensionality reduction of high-dimensional, low-sample-size
datasets. *Knowledge-Based Systems*, 109180.
<https://doi.org/10.1016/j.knosys.2022.109180>

Kosztyan, Z. T., Katona, A. I., Kurbucz, M. T., and Lantos, Z. (2024).
Generalized network-based dimensionality analysis.
*Expert Systems with Applications*, 238, 121779.
<https://doi.org/10.1016/j.eswa.2023.121779>

## Authors

Zsolt T. Kosztyan, Marcell T. Kurbucz, Attila I. Katona, and Zahid Khan.


