# Generalized network-based dimensionality reduction and analysis

`nda` builds interpretable latent variables from communities in a squared-association network. Version 0.3.0 adds a faster and safer computation path while retaining the original numeric method codes.

Core capabilities include:

* full, partial, and semipartial Pearson, rank, robust, and distance associations;
* centrality-weighted latent variables with optional fixed communities;
* sparse networks, multicore pair calculations, FDR screening, and bootstrap loading intervals;
* printable static and interactive network plots and 2D/3D biplots;
* supervised NDRLM models using ordinary least squares, ridge, elastic net, or LOESS.

#### Author

* Zsolt T. Kosztyan
* Marcell T. Kurbucz
* Attila I. Katona
* Zahid Khan

#### Contributor

* Zsolt T. Kosztyan

#### Maintainer

* Zsolt T. Kosztyan

## Installation

Install the development version with:


```
remotes::install_github("kzst/nda")
library(nda)
```

## Quick start

```r
fit <- ndr(swiss, cor_method = "spearman", cor_type = "full",
           centrality = "pagerank")
summary(fit)

# Assignment is silent; printing renders the graph.
network_plot <- plot(fit, interactive = FALSE, show_weights = TRUE)
network_plot

X <- as.data.frame(freeny.x)
Y <- data.frame(income = freeny.y)
model <- ndrlm(Y, X, optimize = FALSE, regression = "ridge")
predict(model, X[1:5, ])
```
