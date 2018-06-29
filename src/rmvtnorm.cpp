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

  // numlin_ : number of linear constraints

  // f2 = f %*% Ri
  // g2 = as.vector(f %*% Mir + g)
  // R = chol(M)
  // Mir = solve(M, r)
  // Ri = solve(R)

  // idea : generate standard normal and transform it back with desired mean and variance

  // number of dimensions
  const int dim = mean.size();

  HmcSampler hmc1(dim);

  // number of linear constraints
  const int numlin = F.rows();

  // Cholesky of covariance
  MatrixXd L {cov.llt().matrixL()};

  // get constraints corresponding to standard normal
  MatrixXd f2 {F * L};
  MatrixXd g2 {(F * mean) + g};

  if (numlin >0){
    for(int i = 0; i < numlin; i++) {
      hmc1.addLinearConstraint(f2.row(i),g2(i));
    }
  }

  hmc1.setInitialValue(L.inverse() * (initial - mean));

  MatrixXd res(n, dim);
  res = res * 0;

  // burn samples
  for (int i = 0; i < burn; i++) {
    hmc1.sampleNext();
  }

  // generate samples and transform them back with appropriate mean and covariance
  for (int i = 0; i < n; i++) {
    MatrixXd temp = (L * hmc1.sampleNext().transpose() + mean);
    res.row(i) = temp.transpose();
  }

  return(res);
}
