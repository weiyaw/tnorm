/* 
 * File:   HmcSampler.cpp
 * Author: aripakman
 * 
 * Created on July 4, 2012, 10:44 AM
 */

#define _USE_MATH_DEFINES   // for the constant M_PI
#include <cmath>
#include <random>
#include <RcppEigen.h>
// #include <magnet/math/quartic.hpp>

#include "my_HmcSampler2.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
// using namespace magnet::math;


// const double EPS = 0.000000001;
// const double HmcSampler::min_t = 0.000000001;
const double EPS = 0.00000001;
const double HmcSampler::min_t = 0.0000000000001;


HmcSampler::HmcSampler(const int& d, const MatrixXd& F, const VectorXd& g)
  : dim(d), F_mat(F), g_vec(g) {}

void HmcSampler::setInitialValue(const VectorXd & initial_value){
  // double check =_verifyConstraints( initial_value );
  /*  if (check <0) {
      cout << "Initial condition out of constraint!" << endl;
      exit(1);
      } else */
  lastSample = initial_value;
}

MatrixXd HmcSampler::sampleNext(bool returnTrace ) {  
  COUNTER++;
  bool bad_new_b = false;
  MatrixXd tracePoints = MatrixXd(dim,0);   //this matrix will only be filled if(returnTrace)

  double T = M_PI/2;		// sample how much time T to move
  VectorXd b = lastSample;
  if (!_verifyConstraints(b)) {
    Rcpp::Rcout << "LAST SAMPLE VIOLATES CONSTRAINTS "
		<< (F_mat * lastSample + g_vec).transpose()
		<< std::endl;
  }

  VectorXd a = VectorXd(dim);   // initial velocity 
  int check_interrupt2 = 0;
  TRIGGER = false;
  if (COUNTER == 44967) {TRIGGER = true;}
  
  while(2){
    if (check_interrupt2++ % 50 == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (check_interrupt2 > 1) {
      Rcpp::Rcout << check_interrupt2 << " long outer loop." << std::endl;
      TRIGGER = true;
    }	

    if (check_interrupt2 == 5) {
      Rcpp::Rcout << COUNTER << std::endl;
      Rcpp::stop("too many loops.");
    }

    double velsign = 0;

    // Sample new initial velocity from normal distributions
    for (int i = 0; i < dim; i++) { 
      a(i) = R::norm_rand();
    }

    double tt = T;		// tt is the time left to move 
    double t1;
    double t2; 			// only relevant to quadratic constraint.
    int cn1, cn2;		// constraint number

    if (TRIGGER) {
      Rcpp::Rcout << "init vel:" << a.transpose() << std::endl;
      Rcpp::Rcout << "prev sam:" << b.transpose() << std::endl;
    }

    bool first_bounce = true;    // for the first move, we do not fear that a small t1 or t2 is due to being in the boundary from the previous bounce.
    int check_interrupt1 = 0;

    while (1){
      // advance the particle for tt = PI/2 

      if (check_interrupt1++ % 50 == 0) {
	Rcpp::checkUserInterrupt();
      }
      // if (check_interrupt1 % 100 == 0) {
      // 	TRIGGER = false;
      // 	Rcpp::Rcout << check_interrupt1 << " long inner loop." << std::endl;
      // }

      t1 = 0; 
      if (!linearConstraints.empty()) {
	_getNextLinearHitTime(a, b, t1, cn1);
      }
      if (TRIGGER) {Rcpp::Rcout << "Hit time: " << t1 << std::endl;}

      double t = t1;		// how much time to move. if t==0, move tt 
      bool linear_hit = true;

      // if no wall to be hit (t==0) or not enough time left to hit the wall (tt<t2) 
      if (t == 0 || tt < t) {
      // if (tt < t) {
	break;
      } else {            
	tt = tt - t;
	VectorXd new_b   = sin(t) * a + cos(t) * b;   // hit location 
	VectorXd hit_vel = cos(t) * a - sin(t) * b;   // hit velocity

	if (TRIGGER) {
	Rcpp::Rcout << "New b : "
		    << (F_mat * new_b + g_vec).transpose() << std::endl;}

       	if (!_verifyConstraints(new_b)) {
	  TRIGGER = true;
	  Rcpp::Rcout << "NEW B VIOLATES " << COUNTER << std::endl;
	  Rcpp::Rcout << "t1: " << t1 << std::endl;
	  Rcpp::Rcout << "old b: "
		      << (F_mat * lastSample + g_vec).transpose() << std::endl;
	  Rcpp::Rcout << "test b : "
		      << (F_mat * (sin(t*0.1) * a + cos(t*0.1) * lastSample) + g_vec).transpose()
		      << std::endl;
	  bad_new_b = true;
	  break;
	}
	b = new_b;
	bad_new_b = false;

	// reflect the velocity and verify that it points in the right direction
	LinearConstraint ql = linearConstraints[cn1];
	double f2 = ((ql.f).dot((ql.f)));
	double alpha = ((ql.f).dot(hit_vel))/f2;
	a = hit_vel - 2 * alpha * (ql.f); // reflected velocity
	velsign = a.dot((ql.f));

	// This occurs rarely, due to numerical instabilities
	if (velsign < 0) {
	  Rcpp::Rcout << velsign << " inner velsign less than 0." << std::endl;
	  break; // get out of while(1). Resample the velocity and start again.
	}
      }
    } //while(1)
    if (bad_new_b) {
      b = lastSample;
      Rcpp::Rcout << "Bad new b. Resample" << std::endl;
      continue;
    }

    if (velsign < 0) {		
      Rcpp::Rcout << velsign << " outer velsign less than 0." << std::endl;
      continue;		 // if velocity is negetive, go to beginning of while(2)
    }

    // make last move of time tt without hitting walls  
    VectorXd bb =  sin(tt) * a + cos(tt) * b;

    // verify that we don't violate the constraints due to a numerical instability

    if (TRIGGER) {Rcpp::Rcout << " bbbb : " << (F_mat * bb + g_vec).transpose() << std::endl;}
    if (_verifyConstraints(bb)) {
      // if (!(_verifyConstraints(bb) < 0)) {
      lastSample = bb;      
      return lastSample;
    }
    Rcpp::Rcout << "BAD BB CONSTRAINTS" << std::endl;
    // at this point we have violated constraints: resample. 

  } // while(2)
}



void HmcSampler::addLinearConstraint(const VectorXd & f, const double & g){
  LinearConstraint newConstraint;
  newConstraint.f = f;
  newConstraint.g = g;
  linearConstraints.push_back(newConstraint);
}
void HmcSampler::addQuadraticConstraint(const MatrixXd & A, const VectorXd & B, const double & C){
  QuadraticConstraint newConstraint;
  newConstraint.A = A;
  newConstraint.B = B;
  newConstraint.C = C;
  quadraticConstraints.push_back(newConstraint);
    
}

void HmcSampler::_getNextLinearHitTime(const VectorXd& a, const VectorXd& b,
				       double& hit_time, int& cn ){
  hit_time = 0;
    
  for (int i=0; i != linearConstraints.size(); i++ ){
    LinearConstraint lc = linearConstraints[i];
    double fa = (lc.f).dot(a);
    double fb = (lc.f).dot(b);
    double u = sqrt(fa*fa + fb*fb);
    if (u > lc.g && u > -lc.g){
      double phi = atan2(-fa,fb); // -PI < phi < PI


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
      
      double t = t1;

      if (t1 == 0) {
      	t = t2;
      } else if (t2 == 0) {
      	t = t1;
      } else {
	t = (t1 < t2? t1 : t2);
      }
      if (TRIGGER) {
      	Rcpp::Rcout << i << " time " << t1 << " " << t2 << " " << t << " " << hit_time << std::endl;
      }
      if (t > min_t && (hit_time == 0 || t < hit_time)) {
	hit_time = t;
	cn = i;
	if (TRIGGER) { Rcpp::Rcout << "CHANGE CONSTRAINT TO " << i << std::endl; }
      }
    }       
  }
}

void HmcSampler::_getNextQuadraticHitTime(const VectorXd & a, const VectorXd & b, double & hit_time, int & cn , const bool first_bounce ){
  //     hit_time=0;
  
  //     double mint;
  //     if (first_bounce) {mint=0;}
  //     else {mint=min_t;}
    
  // for (int i=0; i != quadraticConstraints.size(); i++ ){
        
  //     QuadraticConstraint qc = quadraticConstraints[i];
  //     double q1= - ((a.transpose())*(qc.A))*a;
  //     q1 = q1 + ((b.transpose())*(qc.A))*b;
  //     double q2= (qc.B).dot(b);
  //     double q3= qc.C + a.transpose()*(qc.A)*a;
  //     double q4= 2*b.transpose()*(qc.A)*a;
  //     double q5= (qc.B).dot(a);

  //     double r4 = q1*q1 + q4*q4;
  //     double r3 = 2*q1*q2 + 2*q4*q5;
  //     double r2 = q2*q2 + 2*q1*q3 + q5*q5 -q4*q4;
  //     double r1 = 2*q2*q3 - 2*q4*q5;
  //     double r0=  q3*q3 - q5*q5;

  //     double roots[]={0,0,0,0};
  //     double aa = r3/r4;
  //     double bb = r2/r4;
  //     double cc = r1/r4;
  //     double dd = r0/r4;

  //     //Solve quartics of the form x^4 + aa x^3 + bb x^2 + cc x + dd ==0
  //     int sols = quarticSolve(aa, bb, cc, dd, roots[0], roots[1],  roots[2],  roots[3]);
  //     for (int j=0; j<sols; j++){
  //         double r = roots[j];
  //         if (abs(r) <=1 ){               
  //             double l1 = q1*r*r + q2*r + q3;
  //             double l2 = -sqrt(1-r*r)*(q4*r + q5); 
  //             if (l1/l2 > 0){
  //                 double t = acos(r);
  //                 if (   t> mint      && (hit_time == 0 || t < hit_time)){
  //                    hit_time=t;
  //                    cn=i;                                          
  //                 }                    
  //             }
  //         }            
  //     }                
  // }    
    
    
}


// double HmcSampler::_verifyConstraints(const VectorXd & b){
//   double r =0;
    
//   // for (int i=0; i != quadraticConstraints.size(); i++ ){       
//   //   QuadraticConstraint qc = quadraticConstraints[i];
//   //   double check = ((b.transpose())*(qc.A))*b + (qc.B).dot(b) + qc.C;
//   //   if (i==0 || check < r) {
//   //     r = check;
//   //   }
//   // }

//   for (int i=0; i != linearConstraints.size(); i++ ){       
//     LinearConstraint lc = linearConstraints[i];
//     double check = (lc.f).dot(b) + lc.g;
//     if (i==0 || check < r) {
//       r = check;
//     }
//   }    
//   return r;
// }

bool HmcSampler::_verifyConstraints(const VectorXd& b){

  // for (int i=0; i != quadraticConstraints.size(); i++ ){       
  //   QuadraticConstraint qc = quadraticConstraints[i];
  //   double check = ((b.transpose())*(qc.A))*b + (qc.B).dot(b) + qc.C;
  //   if (i==0 || check < r) {
  //     r = check;
  //   }
  // }

  double constraints_lhs = 1;
  for (int i = 0; i < linearConstraints.size(); i++) { 
    LinearConstraint lc = linearConstraints[i];
    constraints_lhs = (lc.f).dot(b) + lc.g;
    // if (TRIGGER) { Rcpp::Rcout << constraints_lhs << " "; }
    if (!(constraints_lhs > -EPS)) {
      // if (TRIGGER) { Rcpp::Rcout << std::endl; }
      return false;
    }
  }    
  // if (TRIGGER) { Rcpp::Rcout << std::endl; }
  return true;
}

void HmcSampler::_updateTrace( VectorXd const & a,  VectorXd const & b, double const & t, MatrixXd & tracePoints){
  double const stepsize = .01;
  int steps = t/stepsize;
    
  int c = tracePoints.cols();
  tracePoints.conservativeResize(Eigen::NoChange, c+steps+1);
  for (int i=0; i<steps; i++){
    VectorXd bb= sin(i*stepsize)*a + cos(i*stepsize)*b;
    //      cout << bb.transpose() << endl;
    tracePoints.col(c+i) = bb;                    
  }
  VectorXd bb= sin(t)*a + cos(t)*b;
  tracePoints.col(c+steps) = bb;
}
