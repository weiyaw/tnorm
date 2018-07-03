#' Generate truncated normal samples.
#'
#' Generate truncated normal samples using Hamiltonian Monte Carlo technique.
#'
#' This is a modified version of \code{tmg} package in CRAN. The samples will satisfy
#' the following inequality \deqn{F \times x + g > 0} where \eqn{x} is a column vector
#' of the generated sample.
#'
#' @param n number of samples.
#' @param mean mean vector.
#' @param cov variance covariance matrix.
#' @param initial an initial value for the Markov chain.
#' @param F,g constraint matrix and vector.
#'
#' @return a matrix with each row corresponding to a sample.
#'
#' @export
rmvtnorm <- function(n, mean = c(0, 0), cov = diag(2), initial = c(0, 0),
                     F = NULL, g = NULL, burn = 10L) {
    if (is.null(F) || is.null(g)) {
        F <- matrix(1L, nrow = 1, ncol = length(mean))
        g <- Inf
    }

    if (any(F %*% initial + g < 0)) {
        stop("Initial value violates constraints.")
    }
    .Call('_tnorm_rmvtnorm', PACKAGE = 'tnorm', n, mean, cov, initial, F, g, burn)
}

