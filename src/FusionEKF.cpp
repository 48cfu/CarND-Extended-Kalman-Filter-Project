#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */


}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

      float phi = measurement_pack.raw_measurements_[1];
      float px = measurement_pack.raw_measurements_[0] * cos(phi);
      float py = measurement_pack.raw_measurements_[0] * sin(phi);
      float phi_dot = measurement_pack.raw_measurements_[2];
      float angle_between_radial_direction_and_velocity = abs(M_PI/2.0 - 2* phi);
      float v = phi_dot * cos(angle_between_radial_direction_and_velocity);
      float vx = v * sin(phi);
      float vy = v * cos(phi);
      ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    
    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    int n = ekf_.x_.size();
    float sigma = 1;
    MatrixXd P_in = sigma * MatrixXd::Identity(n, n);
    P_in(2, 2) = 10000;
    P_in(3, 3) = 10000;
    MatrixXd Q_in = P_in;
    MatrixXd F_in = MatrixXd::Identity(n, n);
    
    ekf_.Init(ekf_.x_, P_in, F_in, H_laser_, R_laser_, Q_in);
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float noise_ax = 9.0;
  float noise_ay = 9.0;
  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
    // TODO: YOUR CODE HERE
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  (dt_4/4.0)*noise_ax, 0, (dt_3/2.0)*noise_ax, 0,
        0, (dt_4/4.0)*noise_ay, 0, (dt_3/2.0)*noise_ay,
        (dt_3/2.0)*noise_ax, 0, (dt_2)*noise_ax, 0,
        0, (dt_3/2.0)*noise_ay, 0, dt_2*noise_ay;
  cout << ekf_.Q_ << endl;      
  ekf_.Predict();
  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //compute Jacobian at the best estimate
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    // measurement update
    ekf_.R_ = R_laser_;
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
                0, 1, 0, 0;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
