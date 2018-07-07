#include "my_HmcSampler2.h"
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd rmvtnorm(const int& n,
			 const Eigen::VectorXd& mean,
			 const Eigen::MatrixXd& cov,
			 const Eigen::VectorXd& initial,
			 const Eigen::MatrixXd& F,
			 const Eigen::VectorXd& g,
			 const int& burn = 10) {

  using Eigen::VectorXd;
  using Eigen::MatrixXd;

  // idea : generate standard normal and transform it back with desired mean and variance

  // number of dimensions
  const int dim = mean.size();
  
  if (cov.rows() != cov.cols()) {
    Rcpp::stop("Covariance must be a square matrix.");
  }

  if (cov.rows() != dim) {
    Rcpp::stop("Dimensions of mean and covariance not matching.");
  }

  if (initial.size() != dim) {
    Rcpp::stop("Dimensions of mean and initial value not matching.");
  }
  
  if (F.cols() != dim) {
    Rcpp::stop("Dimension of mean and the columns of F not matching.");
  }
  
  if (F.rows() != g.size()) {
    Rcpp::stop("Inconsistent dimensions of F and g.");
  }
  
  // number of linear constraints
  const int numlin = F.rows();

  if (initial.size() != dim) {
    Rcpp::stop("Dimensions of mean and initial value not matching.");
  }

  // covariance is "symmetricalised"
  MatrixXd cov_sym {(cov.transpose() + cov) / 2};
  Eigen::LLT<MatrixXd> cov_llt {cov_sym.llt()};

  if (cov_llt.info() == Eigen::NumericalIssue) {
    Rcpp::stop("Non positive definitie covariance.");
  }
  
  // Cholesky of covariance
  MatrixXd L {cov_llt.matrixL()};

  // get constraints corresponding to standard normal
  MatrixXd f2 {F * L};
  MatrixXd g2 {(F * mean) + g};

  // initialise hmc sampler
  HmcSampler hmc1(L.inverse() * (initial - mean), dim, f2, g2);

  if (numlin > 0){
    for(int i = 0; i < numlin; i++) {
      hmc1.addLinearConstraint(f2.row(i),g2(i));
    }
  }

  MatrixXd res(n, dim);
  res = res * 0;

  // burn samples
  for (int i = 0; i < burn; i++) {
    hmc1.sampleNext();
  }

  // generate samples and transform them back with appropriate mean and covariance
  for (int i = 0; i < n; i++) {
    res.row(i) = (L * hmc1.sampleNext() + mean).transpose();
  }

  return(res);
}
