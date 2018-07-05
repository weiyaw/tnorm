/* 
 * File:   HmcSampler.h
 * Author: aripakman
 * Modified by: Kenyon Ng
 * Created on July 4, 2012, 10:44 AM
 */

#ifndef HMCSAMPLER_H
#define	HMCSAMPLER_H

#define _USE_MATH_DEFINES

#include <RcppEigen.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct LinearConstraint{
  VectorXd f;
  double  g;
};

struct QuadraticConstraint{
    MatrixXd A;
    VectorXd B;
    double  C;    
};


class HmcSampler {
 public:
    
  HmcSampler(const VectorXd& init_val, const int& d, const MatrixXd& F, const VectorXd& g);
  void addLinearConstraint(const VectorXd & f, const double & g);
  MatrixXd sampleNext();
    
 private:
  int COUNTER = 0;
  bool TRIGGER;
  MatrixXd F_mat;
  VectorXd g_vec;
    
  int dim;
  VectorXd last_sample;
  static const double EPS;
  std::vector<LinearConstraint> linearConstraints;
    
  void _getNextLinearHitTime(const VectorXd& a, const VectorXd& b,  double& t, int& cn );
  bool _verifyConstraints(const VectorXd& b);
  
  /*     void _updateTrace( VectorXd const & a,  VectorXd const & b, double const & tt, MatrixXd & tracePoints); */
};

#endif	/* HMCSAMPLER_H */

