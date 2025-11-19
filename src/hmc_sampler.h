/* 
 * File:   HmcSampler.h
 * Author: aripakman
 * Modified by: Kenyon Ng
 * Created on July 4, 2012, 10:44 AM
 */

#ifndef HMCSAMPLER_H
#define	HMCSAMPLER_H

#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <RcppEigen.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif

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
  void _sampleVelocity(VectorXd& velocity) const;
  bool _advanceUntilHit(VectorXd& position,
                        VectorXd& velocity,
                        double& time_left,
                        int& hit_constraint);
  bool _reflectVelocity(const LinearConstraint& constraint,
                        const VectorXd& hit_velocity,
                        VectorXd& velocity);
  
  /*     void _updateTrace( VectorXd const & a,  VectorXd const & b, double const & tt, MatrixXd & tracePoints); */
};

#endif	/* HMCSAMPLER_H */
