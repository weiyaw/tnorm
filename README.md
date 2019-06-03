# Truncated Gaussian sampler using exact Hamiltonian Monte Carlo

This is an R package that implements Hamiltonian Monte Carlo
([Pakman and Paninski, 2014](https://arxiv.org/abs/1208.4118)) to sample from a
Gaussian distribution truncated with linear constraints. The sampler was written
in `C++` for fast computation and connected to `R` using `Rcpp`.

This package reuses part of the code from the `tmg` package
([link](https://cran.r-project.org/web/packages/tmg/index.html) to CRAN) and
fixes an edge case in `tmg`.

You can painlessly install this package using `devtools`, which can be
downloaded and installed with the following one-liner:
```
install.packages("devtools")
```

Then, the installation of this package is just
```
devtools::install_github("weiyaw/blackbox")
```
