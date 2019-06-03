#' Generate truncated normal samples.
#'
#' Generate truncated normal samples using Hamiltonian Monte Carlo technique.
#'
#' This is a modified version of \code{tmg} package in CRAN. The samples will satisfy
#' the following inequality \deqn{F \times x + g >= 0} where \eqn{x} is a column vector
#' of the generated sample.
#'
#' No check for symmetry is performed on covariance. The initial value must
#' strictly satisfy the inequality constraints, although the generated samples
#' only satisfy them weakly.
#'
#' @param n number of samples.
#' @param mean mean vector.
#' @param cov variance covariance matrix.
#' @param initial an initial value for the Markov chain.
#' @param F,g constraint matrix and vector.
#'
#' @return A matrix with each row corresponding to a sample.
#'
#' @examples
#' plot(rmvtnorm(100))
#'
#' plot(rmvtnorm(100, F = diag(2), g = c(1, -1), initial = c(0, 2)))
#'
#' plot(rmvtnorm(100, mean = c(1, 1), cov = matrix(c(1, 0.5, 0.5, 1), 2), F =
#' diag(2), g = c(1, -1), initial = c(0, 2)))
#'
#' @export
rmvtnorm <- function(n, mean = c(0, 0), cov = diag(length(mean)),
                     initial = rep(0, length(mean)), F = NULL, g = NULL, burn = 10L) {
    if (is.null(F) || is.null(g)) {
        F <- matrix(1L, nrow = 1, ncol = length(mean))
        g <- Inf
    }

    if (!is.vector(mean) && !(is.matrix(mean) && min(dim(mean)) == 1)) {
        stop("Mean must be a vector.")
    }
    if (!is.matrix(cov) || ncol(cov) != nrow(cov)) {
        stop("Covariance must be a square matrix.")
    }
    if (ncol(cov) != length(mean)) {
        stop("Inconsistant dimensions of the mean and covariance.")
    }
    if (length(initial) != length(mean)) {
        stop("Inconsistant dimensions of the mean and the initial value.")
    }
    if (any(F %*% initial + g < 0)) {
        stop("Initial value violates constraints.")
    }
    .Call('_tnorm_rmvtnorm', PACKAGE = 'tnorm', n, mean, cov, initial, F, g, burn)
}


