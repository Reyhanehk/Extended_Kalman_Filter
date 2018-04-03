#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  F_ = MatrixXd(4, 4);
  P_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  //measurement matrix laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //the initial transition matrix F_
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  ekf_.F_ = F_;

  //state covariance matrix P
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  ekf_.P_ = P_;

  ekf_.Q_ = MatrixXd(4, 4);

 //Set the process and measurement noises
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //Convert radar from polar to cartesian coordinates and initialize state.

      //polar measurements
      float rho = measurement_pack.raw_measurements_[0]; //range
      float phi = measurement_pack.raw_measurements_[1]; //bearing
      float rho_dot = measurement_pack.raw_measurements_[2]; //veloicty of rho

      //changing coordinates polar -> cartesian
      float x = rho * cos(phi);
      float y = rho * sin(phi);

      //changing coordinates polar -> cartesian for velocity measurements
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);

      //intializing Radar measurements
      ekf_.x_ << x, y, vx, vy;

    }else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     //Initialize state from Lidar
     ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    //we can't let x and y be zero.
      if ( fabs(ekf_.x_(0) + ekf_.x_(1)) < 1e-4){
          ekf_.x_(0) = 1e-4;
          ekf_.x_(1) = 1e-4;
      }


    // Save the initial timestamp for dt calculation
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    // Update the state transition matrix F according to the new elapsed time.
    // - Time is measured in seconds.
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

    //updating the process covariance matrix Q

    //Modify the Q matrix so that the time is integrated
    ekf_.Q_(0,0) = (dt_4/4)*noise_ax_;
    ekf_.Q_(0,2) = (dt_3/2)*noise_ax_;

    ekf_.Q_(1,1) = (dt_4/4)*noise_ay_;
    ekf_.Q_(1,3) = (dt_3/2)*noise_ay_;

    ekf_.Q_(2,0) = (dt_3/2)*noise_ax_;
    ekf_.Q_(2,2) = dt_2*noise_ax_;

    ekf_.Q_(3,1) = (dt_3/2)*noise_ay_;
    ekf_.Q_(3,3) = dt_2*noise_ay_;


    //Make Prediction only if we the dt is big enough.
    if ( dt > 1e-3 )
    {
        ekf_.Predict();
    }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // Using a Jacobian Matrix instead of H
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
