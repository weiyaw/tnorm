#define _USE_MATH_DEFINES   // for the constant M_PI
#include <cmath>
#include <cstddef>
#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <RcppEigen.h>
#ifdef __clang__
#  pragma clang diagnostic pop
#endif
#include "hmc_sampler.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double HmcSampler::EPS = 0.00000001;

HmcSampler::HmcSampler(const VectorXd& init_val,
		       const int& d,
		       const MatrixXd& F,
		       const VectorXd& g)
  : F_mat(F), g_vec(g), dim(d), last_sample(init_val) {}

MatrixXd HmcSampler::sampleNext() {  
  COUNTER++;
  TRIGGER = false;

  constexpr double kTimeHorizon = M_PI / 2.0;
  constexpr int kInterruptInterval = 50;
  constexpr int kMaxOuterIterations = 500;

  VectorXd velocity(dim);
  int interrupt_counter = 0;

  for (int iteration = 0;; ++iteration) {
    if (interrupt_counter++ % kInterruptInterval == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (iteration > kMaxOuterIterations) {
      Rcpp::stop("Too many loops.");
    }

    VectorXd position = last_sample;
    _sampleVelocity(velocity);

    double time_left = kTimeHorizon;
    int hit_constraint = -1;

    if (!_advanceUntilHit(position, velocity, time_left, hit_constraint)) {
      continue;
    }

    VectorXd candidate = std::sin(time_left) * velocity + std::cos(time_left) * position;
    if (_verifyConstraints(candidate)) {
      last_sample = candidate;
      return last_sample;
    }
  }
}



void HmcSampler::addLinearConstraint(const VectorXd& f, const double& g){
  LinearConstraint newConstraint;
  newConstraint.f = f;
  newConstraint.g = g;
  linearConstraints.push_back(newConstraint);
}

void HmcSampler::_getNextLinearHitTime(const VectorXd& a, const VectorXd& b,
				       double& hit_time, int& cn){
  hit_time = -1;
  int prev_cn = cn;		// constraint causing the previous hit

  const double two_pi = 2.0 * M_PI;

  for (std::size_t i = 0; i < linearConstraints.size(); ++i) {
    const LinearConstraint& lc = linearConstraints[i];
    const double fa = lc.f.dot(a);
    const double fb = lc.f.dot(b);
    const double u = std::sqrt(fa * fa + fb * fb);

    if (!(u > lc.g && u > -lc.g)) {
      continue;
    }

    const double phi = std::atan2(-fa, fb); // -PI < phi < PI
    double t1 = std::acos(-lc.g / u) - phi; // -PI < t1 < 2PI
    double t2 = -t1 - 2 * phi;	//  -2PI < t2 < PI

    if (t1 < 0) {
      t1 += two_pi;
    }
    if (std::abs(t1) < EPS || std::abs(t1 - two_pi) < EPS) {
      t1 = 0;
    }

    if (t2 < 0) {
      t2 += two_pi;
    }
    if (std::abs(t2) < EPS || std::abs(t2 - two_pi) < EPS) {
      t2 = 0;
    }

    double t = -1;
    if (prev_cn == static_cast<int>(i)) {
      t = (t1 > t2 ? t1 : t2);
    } else if (t1 == 0 || t2 == 0) {
      const double midpoint = std::fabs(t2 - t1) / 2.0;
      if (u * std::cos(midpoint + phi) > -lc.g) {
        t = (t1 > t2 ? t1 : t2);
      } else {
        t = 0;
      }
    } else {
      t = (t1 < t2 ? t1 : t2);
    }

    if (hit_time < 0 || t < hit_time) {
      hit_time = t;
      cn = static_cast<int>(i);
    }
  }
}


bool HmcSampler::_verifyConstraints(const VectorXd& b){

  for (const auto& lc : linearConstraints) {
    const double constraints_lhs = lc.f.dot(b) + lc.g;
    if (!(constraints_lhs > -EPS)) {
      return false;
    }
  }    
  return true;
}

void HmcSampler::_sampleVelocity(VectorXd& velocity) const {
  for (int i = 0; i < dim; ++i) {
    velocity(i) = R::norm_rand();
  }
}

bool HmcSampler::_advanceUntilHit(VectorXd& position,
                                  VectorXd& velocity,
                                  double& time_left,
                                  int& hit_constraint) {
  constexpr int kInterruptInterval = 50;
  int interrupt_counter = 0;

  while (true) {
    if (interrupt_counter++ % kInterruptInterval == 0) {
      Rcpp::checkUserInterrupt();
    }

    double hit_time = -1;
    if (!linearConstraints.empty()) {
      _getNextLinearHitTime(velocity, position, hit_time, hit_constraint);
    }

    if (hit_time < 0 || time_left < hit_time) {
      hit_constraint = -1;
      break;
    }

    if (hit_time < EPS) {
      return false;
    }

    time_left -= hit_time;
    const VectorXd hit_velocity = std::cos(hit_time) * velocity - std::sin(hit_time) * position;
    position = std::sin(hit_time) * velocity + std::cos(hit_time) * position;

    const LinearConstraint& constraint = linearConstraints.at(hit_constraint);
    if (!_reflectVelocity(constraint, hit_velocity, velocity)) {
      return false;
    }
  }

  return true;
}

bool HmcSampler::_reflectVelocity(const LinearConstraint& constraint,
                                  const VectorXd& hit_velocity,
                                  VectorXd& velocity) {
  const double norm_sq = constraint.f.dot(constraint.f);
  const double alpha = constraint.f.dot(hit_velocity) / norm_sq;
  velocity = hit_velocity - 2 * alpha * constraint.f;

  const double velsign = velocity.dot(constraint.f);
  return velsign >= 0;
}
