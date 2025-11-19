# Multivariate Truncated Normal Sampler

`tnorm` provides exact Hamiltonian Monte Carlo (HMC) sampling for multivariate
Gaussian distributions subject to linear inequality constraints. The core
sampler is written in C++ following [Pakman & Paninski
(2014)](https://arxiv.org/pdf/1208.4118) and exposed to R via Rcpp for
convenient statistical workflows. It started as a ground-up rewrite of the
unmaintained [`tmg`](https://cran.r-project.org/package=tmg) package.

## Installation

```r
install.packages("remotes")
remotes::install_github("weiyaw/tnorm")
```

## Example

Draw 2,000 samples from a bivariate standard normal distribution truncated to the positive quadrant.
The rows of the returned matrix are independent draws after a burn-in period.

```r
library(tnorm)

# Inequalities Fx + g >= 0 enforce x1 >= 0, x2 >= 0
F <- diag(2)
g <- c(0, 0)

samples <- rmvtnorm(
  n = 2000,
  mean = c(0, 0),
  cov = diag(2),
  initial = c(0.2, 0.2),
  F = F,
  g = g,
  burn = 500
)

head(samples)
colMeans(samples)  # ≈ sqrt(2 / pi)
diag(var(samples)) # ≈ 1 - 2 / pi
```

Because the sampler operates in a transformed standard-normal space, the linear constraints are
enforced exactly and the moments of the samples match known analytical results (see
`tests/testthat/test-rmvtnorm.R` for references).


## Reference
Pakman, A. and Paninski, L. (2014) Exact Hamiltonian Monte Carlo for Truncated Multivariate Gaussians, _Journal of Computational and Graphical Statistics_, 23(2), pp. 518–542
