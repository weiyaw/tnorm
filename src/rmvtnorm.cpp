#include "my_HmcSampler2.h"
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd rmvtnorm(const int n,
			 const Eigen::VectorXd mean,
			 const Eigen::MatrixXd cov,
			 const Eigen::VectorXd initial, const Eigen::MatrixXd F,
			 const Eigen::VectorXd g,
			 const int burn = 10) {

  using Eigen::VectorXd;
  using Eigen::MatrixXd;

  // idea : generate standard normal and transform it back with desired mean and variance

  // number of dimensions
  const int dim = mean.size();

  // number of linear constraints
  const int numlin = F.rows();

  // Cholesky of covariance
  MatrixXd L {cov.llt().matrixL()};

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
