#include "hmc_sampler.h"
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <RcppEigen.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

namespace {

using Eigen::LLT;
using Eigen::MatrixXd;
using Eigen::VectorXd;

void validateDimensions(const MatrixXd& cov,
                        const VectorXd& mean,
                        const VectorXd& initial,
                        const MatrixXd& F,
                        const VectorXd& g) {
  const int dim = mean.size();

  if (cov.rows() != cov.cols()) {
    Rcpp::stop("Covariance must be a square matrix.");
  }
  if (cov.rows() != dim) {
    Rcpp::stop("Dimensions of mean and covariance do not match.");
  }
  if (initial.size() != dim) {
    Rcpp::stop("Dimensions of mean and initial value do not match.");
  }
  if (F.cols() != dim) {
    Rcpp::stop("Dimensions of mean and columns of F do not match.");
  }
  if (F.rows() != g.size()) {
    Rcpp::stop("Inconsistent dimensions of F and g.");
  }
}

MatrixXd symmetrise(const MatrixXd& cov) {
  return (cov + cov.transpose()) / 2.0;
}

LLT<MatrixXd> safeCholesky(const MatrixXd& cov) {
  LLT<MatrixXd> cov_llt(cov);
  if (cov_llt.info() == Eigen::NumericalIssue) {
    Rcpp::stop("Covariance must be positive definite.");
  }
  return cov_llt;
}

MatrixXd transformConstraintMatrix(const MatrixXd& F, const MatrixXd& L) {
  return F * L;
}

VectorXd transformConstraintOffset(const MatrixXd& F,
                                   const VectorXd& mean,
                                   const VectorXd& g) {
  return (F * mean) + g;
}

VectorXd transformInitialPoint(const MatrixXd& L,
                               const VectorXd& mean,
                               const VectorXd& initial) {
  const VectorXd diff = initial - mean;
  return L.triangularView<Eigen::Lower>().solve(diff);
}

void addConstraints(HmcSampler& sampler,
                    const MatrixXd& constraints,
                    const VectorXd& offsets) {
  const int num_constraints = constraints.rows();
  for (int i = 0; i < num_constraints; ++i) {
    sampler.addLinearConstraint(constraints.row(i), offsets(i));
  }
}

void runBurnIn(HmcSampler& sampler, const int burn) {
  for (int i = 0; i < burn; ++i) {
    sampler.sampleNext();
  }
}

MatrixXd sampleChain(HmcSampler& sampler,
                     const MatrixXd& L,
                     const VectorXd& mean,
                     const int n) {
  const int dim = mean.size();
  MatrixXd draws = MatrixXd::Zero(n, dim);
  for (int i = 0; i < n; ++i) {
    draws.row(i) = (L * sampler.sampleNext() + mean).transpose();
  }
  return draws;
}

} // namespace

// [[Rcpp::export]]
Eigen::MatrixXd rmvtnorm(const int& n,
			 const Eigen::VectorXd& mean,
			 const Eigen::MatrixXd& cov,
			 const Eigen::VectorXd& initial,
			 const Eigen::MatrixXd& F,
			 const Eigen::VectorXd& g,
			 const int& burn = 10) {

  validateDimensions(cov, mean, initial, F, g);
  if (burn < 0) {
    Rcpp::stop("Burn-in must be non-negative.");
  }

  const int dim = mean.size();
  const MatrixXd cov_sym = symmetrise(cov);
  const LLT<MatrixXd> cov_llt = safeCholesky(cov_sym);
  const MatrixXd L = cov_llt.matrixL();

  const MatrixXd transformed_F = transformConstraintMatrix(F, L);
  const VectorXd transformed_g = transformConstraintOffset(F, mean, g);
  const VectorXd initial_std = transformInitialPoint(L, mean, initial);

  HmcSampler sampler(initial_std, dim, transformed_F, transformed_g);

  if (transformed_F.rows() > 0) {
    addConstraints(sampler, transformed_F, transformed_g);
  }

  runBurnIn(sampler, burn);
  return sampleChain(sampler, L, mean, n);
}
