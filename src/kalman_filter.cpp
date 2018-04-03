#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // Update for the Laser
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose();
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Update for Radar

  // Finding px and py
  double px = x_(0);
  double py = x_(1);

  // Finding vx and vy
  double vx = x_(2);
  double vy = x_(3);

  // Calculating rho
  float rho = pow(pow(px,2)+pow(py,2),0.5);

  // Checking if rho is zero 
  // if (!rho){
  //   return;
  // }

  // Calculate the Jacobian Matrix
  H_ = tools.CalculateJacobian( x_ );
  
  // Computing h(x)
  VectorXd h = VectorXd(3);
  h << sqrt(rho), atan(py/px), (px * vx + py * vy)/sqrt(rho);

  // Measurement update
  VectorXd y = z - h;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
}
