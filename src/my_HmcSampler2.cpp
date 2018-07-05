/* 
 * File:   HmcSampler.cpp
 * Author: aripakman
 * 
 * Created on July 4, 2012, 10:44 AM
 */

#define _USE_MATH_DEFINES   // for the constant M_PI
#include <cmath>
#include <RcppEigen.h>
#include "my_HmcSampler2.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double HmcSampler::EPS = 0.00000001;

HmcSampler::HmcSampler(const VectorXd& init_val,
		       const int& d,
		       const MatrixXd& F,
		       const VectorXd& g)
  : last_sample(init_val), dim(d), F_mat(F), g_vec(g) {}

MatrixXd HmcSampler::sampleNext() {  
  COUNTER++;
  bool RESAMPLE = false;
  TRIGGER = false;

  double T = M_PI/2;		// sample how much time T to move

  VectorXd a = VectorXd(dim);   // initial velocity 
  int check_interrupt2 = 0;
  
  while(2){
    if (check_interrupt2++ % 50 == 0) {
      Rcpp::checkUserInterrupt();
    }

    VectorXd b = last_sample;
    int hit_wall = -1;		   // the constraint hit by the particle (-1 if not hit)

    if (check_interrupt2 > 5) {
      Rcpp::Rcout << COUNTER << std::endl;
      Rcpp::stop("Too many loops.");
    }
    
    // Sample new initial velocity from normal distributions
    for (int i = 0; i < dim; i++) { 
      a(i) = R::norm_rand();
    }

    double t_left = T;		// t_left is the time left to move 
    int check_interrupt1 = 0;	// check for user interruptions

    while (1){
      // advance the particle for t_left = PI/2 

      if (check_interrupt1++ % 50 == 0) {
	Rcpp::checkUserInterrupt();
      }

      double hit_t = -1; 	// time to hit a wall from the latest b

      if (!linearConstraints.empty()) {
	_getNextLinearHitTime(a, b, hit_t, hit_wall);
      }

      if (hit_t < 0 || t_left < hit_t) {
	// if no wall to be hit (t = -1) or not enough time left to hit the wall (t_left < t) 
	hit_wall = -1;		// no wall is hit
	break;
      } else if (hit_t < EPS) {
	// if getting stuck at a place (i.e. moving will violate constraints.)
	RESAMPLE = true;	// resample a new velocity
	break;
      } else {            
	t_left = t_left - hit_t;
	b = sin(hit_t) * a + cos(hit_t) * b;		    // hit location 
	VectorXd hit_vel = cos(hit_t) * a - sin(hit_t) * b; // hit velocity

	// reflect the velocity and verify that it points in the right direction
	LinearConstraint ql = linearConstraints.at(hit_wall);
	double f2 = ((ql.f).dot((ql.f)));
	double alpha = ((ql.f).dot(hit_vel)) / f2;
	a = hit_vel - 2 * alpha * (ql.f); // reflected velocity
	double velsign = a.dot((ql.f));

	// This occurs rarely, due to numerical instabilities
	if (velsign < 0) {
	  RESAMPLE = true;
	  break;  // get out of while(1). Resample the velocity and start again.
	}
      }
    } // while(1)

    // Resample a velocity if velsign < 0 / b breaks constraints / getting stuck.
    if (RESAMPLE) {
      RESAMPLE = false;
      continue;
    }
      
    // make last move of time without hitting walls  
    VectorXd final_b = sin(t_left) * a + cos(t_left) * b;

    // verify that we don't violate the constraints due to a numerical instability
    if (_verifyConstraints(final_b)) {
      last_sample = final_b;
      return last_sample;
    }

    Rcpp::Rcout << "BAD FINAL B" << std::endl;
    // at this point we have violated constraints: resample velocity.

  } // while(2)
}



void HmcSampler::addLinearConstraint(const VectorXd & f, const double & g){
  LinearConstraint newConstraint;
  newConstraint.f = f;
  newConstraint.g = g;
  linearConstraints.push_back(newConstraint);
}

void HmcSampler::_getNextLinearHitTime(const VectorXd& a, const VectorXd& b,
				       double& hit_time, int& cn){
  hit_time = -1;
  int prev_cn = cn;		// constraint causing the previous hit

  for (int i = 0; i != linearConstraints.size(); i++){
    LinearConstraint lc = linearConstraints[i];
    double fa = (lc.f).dot(a);
    double fb = (lc.f).dot(b);
    double u = sqrt(fa*fa + fb*fb);
    if (u > lc.g && u > -lc.g) {
      double phi = atan2(-fa, fb); // -PI < phi < PI

      // solve eqn 2.23 for t, first case
      // -PI < t1 + phi < PI
      double t1 = acos(-lc.g/u) - phi; // -PI < t1 < 2PI

      // solve eqn 2.23 for t, second case,
      // since cos(t1 + phi) == cos(-(t2 + phi))
      double t2 = -t1 - 2*phi;	//  -2PI < t2 < PI

      if (t1 < 0){		//  0 < t1 < 2PI
	t1 += 2*M_PI;
      }

      // if t1 is close to 0 or 2pi
      if (abs(t1) < EPS || abs(t1 - 2*M_PI) < EPS) {
	t1 = 0;			//  0 <= t1 < 2PI
      }
      
      if (t2 < 0) {
	t2 += 2*M_PI;		// 0 < t2 < 2PI
      }

      if (abs(t2) < EPS || abs(t2 - 2*M_PI) < EPS) {
	t2 = 0;			// 0 <= t2 < 2PI
      }       
      
      double t = -1;

      if (prev_cn == i) {
	t = (t1 > t2? t1 : t2);
      } else {
	t = (t1 < t2? t1 : t2);
      }
      // if (TRIGGER) {
      // 	Rcpp::Rcout << i << " time " << t1 << " " << t2
      // 		    << " " << t << " " << hit_time << std::endl;
      // }
      if (hit_time < 0 || t < hit_time) {
	hit_time = t;
	cn = i;
	// if (TRIGGER) {
	//   Rcpp::Rcout << "CHANGE CONSTRAINT TO " << i << std::endl;
	// }
      }
    }
    // Rcpp::stop("Invalid hit time."); 
  }
}


bool HmcSampler::_verifyConstraints(const VectorXd& b){

  double constraints_lhs = 1;
  for (int i = 0; i < linearConstraints.size(); i++) { 
    LinearConstraint lc = linearConstraints[i];
    constraints_lhs = (lc.f).dot(b) + lc.g;
    if (!(constraints_lhs > -EPS)) {
      return false;
    }
  }    
  return true;
}
