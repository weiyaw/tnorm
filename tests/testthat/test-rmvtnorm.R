context("rmvtnorm / HMC sampler")

# Helper functions ----------------------------------------------------------

truncated_normal_moments <- function(mean, sd, lower = -Inf, upper = Inf) {
  alpha <- (lower - mean) / sd
  beta <- (upper - mean) / sd
  Z <- pnorm(beta) - pnorm(alpha)
  phi_alpha <- dnorm(alpha)
  phi_beta <- dnorm(beta)
  alpha_phi <- if (is.infinite(alpha)) 0 else alpha * phi_alpha
  beta_phi <- if (is.infinite(beta)) 0 else beta * phi_beta
  expected_mean <- mean + (phi_alpha - phi_beta) / Z * sd
  expected_var <- sd^2 * (1 + (alpha_phi - beta_phi) / Z -
    ((phi_alpha - phi_beta) / Z)^2)
  list(mean = expected_mean, var = expected_var)
}

summarise_samples <- function(samples) {
  samples <- as.matrix(samples)
  list(mean = colMeans(samples), cov = var(samples))
}

# One-dimensional scenarios -------------------------------------------------

test_that("1D half-line truncation matches analytical moments", {
  set.seed(101)

  lower <- 1
  mu <- 0
  sigma <- 1

  samples <- rmvtnorm(
    n = 8000,
    mean = mu,
    cov = matrix(sigma^2),
    initial = 2,
    F = matrix(1, nrow = 1),
    g = -lower,
    burn = 800
  )

  moments <- truncated_normal_moments(mu, sigma, lower = lower)
  stats <- summarise_samples(samples)

  expect_equal(stats$mean[1], moments$mean, tolerance = 0.08)
  expect_equal(stats$cov[1, 1], moments$var, tolerance = 0.08)
})

test_that("1D bounded interval matches analytical moments", {
  set.seed(123)

  mu <- 1.2
  sigma <- 0.8
  lower <- 0.25
  upper <- 2

  F <- rbind(1, -1)
  g <- c(-lower, upper)

  samples <- rmvtnorm(
    n = 5000,
    mean = mu,
    cov = matrix(sigma^2),
    initial = 0.5,
    F = F,
    g = g,
    burn = 500
  )

  moments <- truncated_normal_moments(mu, sigma, lower = lower, upper = upper)
  stats <- summarise_samples(samples)

  expect_equal(stats$mean[1], moments$mean, tolerance = 0.05)
  expect_equal(stats$cov[1, 1], moments$var, tolerance = 0.06)
})

# Multivariate scenarios ----------------------------------------------------

test_that("Positive quadrant truncation preserves marginal independence", {
  set.seed(234)

  mean_vec <- c(0, 0)
  cov_mat <- diag(2)

  F <- diag(2)
  g <- c(0, 0)

  samples <- rmvtnorm(
    n = 9000,
    mean = mean_vec,
    cov = cov_mat,
    initial = c(1, 1),
    F = F,
    g = g,
    burn = 900
  )

  stats <- summarise_samples(samples)
  moments <- truncated_normal_moments(0, 1, lower = 0)
  expect_equal(stats$mean[1], moments$mean, tolerance = 0.05)
  expect_equal(stats$mean[2], moments$mean, tolerance = 0.05)
  expect_equal(stats$cov[1, 1], moments$var, tolerance = 0.05)
  expect_equal(stats$cov[2, 2], moments$var, tolerance = 0.05)
  expect_equal(stats$cov[1, 2], 0, tolerance = 0.05)
})

test_that("Independent rectangular constraints match analytical moments", {
  set.seed(456)

  mean_vec <- c(-0.2, 0.5)
  cov_mat <- diag(c(1.6, 0.4))
  lower_bounds <- c(0, -1)
  upper_bounds <- c(Inf, 1)

  F <- rbind(
    c(1, 0),
    c(0, 1),
    c(0, -1)
  )
  g <- c(-lower_bounds[1], -lower_bounds[2], upper_bounds[2])

  samples <- rmvtnorm(
    n = 6000,
    mean = mean_vec,
    cov = cov_mat,
    initial = c(0.2, 0),
    F = F,
    g = g,
    burn = 800
  )

  stats <- summarise_samples(samples)

  x1_mom <- truncated_normal_moments(mean_vec[1], sqrt(cov_mat[1, 1]), lower = lower_bounds[1])
  x2_mom <- truncated_normal_moments(mean_vec[2], sqrt(cov_mat[2, 2]), lower = lower_bounds[2], upper = upper_bounds[2])

  expect_equal(stats$mean[1], x1_mom$mean, tolerance = 0.05)
  expect_equal(stats$cov[1, 1], x1_mom$var, tolerance = 0.06)

  expect_equal(stats$mean[2], x2_mom$mean, tolerance = 0.05)
  expect_equal(stats$cov[2, 2], x2_mom$var, tolerance = 0.06)

  expect_equal(stats$cov[1, 2], 0, tolerance = 0.04)
})
