# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rmvtnorm <- function(n, mean, cov, initial, F, g, burn = 10L) {
    .Call('_tnorm_rmvtnorm', PACKAGE = 'tnorm', n, mean, cov, initial, F, g, burn)
}

