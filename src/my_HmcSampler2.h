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


class HmcSampler   {
public:
    
    HmcSampler(const int & d);

    void setInitialValue(const VectorXd & initial);
    void addLinearConstraint(const VectorXd & f, const double & g);
    void addQuadraticConstraint(const MatrixXd & A, const VectorXd & B, const double & C);
    MatrixXd sampleNext(bool returnTrace = false);
    
private:
    int dim;
    VectorXd lastSample;    
    static const double min_t; 
    std::vector<LinearConstraint> linearConstraints;
    std::vector<QuadraticConstraint> quadraticConstraints;
    
    void _getNextLinearHitTime(const VectorXd & a, const VectorXd & b,  double & t, int & cn );
    void _getNextQuadraticHitTime(const VectorXd & a, const VectorXd & b, double & t, int & cn, const bool );
    /* double _verifyConstraints(const VectorXd & b); */
    bool _verifyConstraints(const VectorXd & b);
    void _updateTrace( VectorXd const & a,  VectorXd const & b, double const & tt, MatrixXd & tracePoints);
};

#endif	/* HMCSAMPLER_H */

