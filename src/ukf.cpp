#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  // P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);  // safer initialization for covariance matrix

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.9;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // time when the state is true, in us
  time_us_ = 0.0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  
  double weight_0 = lambda_/(lambda_+n_aug_);
  double weight = 0.5/(lambda_+n_aug_);

  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  
  if (!is_initialized_) {
		std::cout << "Kalman Filter Initialization " << std::endl;
	  
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      
      double range = meas_package.raw_measurements_(0);
      double az = meas_package.raw_measurements_(1);
      double p_x = range*cos(az);
      double p_y = range*sin(az);
      
      x_ <<  range*cos(az), range*sin(az), 0, 0, 0; 
      
      P_(0,0) = std_radr_*std_radr_*cos(az)*cos(az) + range*range*sin(az)*sin(az)*std_radphi_*std_radphi_;          
      P_(1,1) = std_radr_*std_radr_*sin(az)*sin(az) + range*range*cos(az)*cos(az)*std_radphi_*std_radphi_;   
      P_(2,2) = 10.0;
      P_(3,3) = (2*M_PI)*(2*M_PI);  // no info about the azimuth
      P_(4,4) = 0.5*0.5;            // yaw rate (rad/s)^2  

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      std::cout << "Initialized by lidar\n";
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

      /*
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1000.0, 0, 0,
            0, 0, 0, (2.0*M_PI)*(2.0*M_PI), 0,
            0, 0, 0, 0, 0.1*0.1;
      */
      P_(0,0) = std_laspx_*std_laspx_;
      P_(1,1) = std_laspy_*std_laspy_;
      P_(2,2) = 50.0;     // speed
      P_(3,3) = (2.0*M_PI)*(2.0*M_PI);  // no info about the azimuth
      P_(4,4) = 0.25;            // yaw rate (rad/s)^2 
      
    }

		time_us_ = meas_package.timestamp_;
		is_initialized_ = true;
		return;
	}

	// compute the time elapsed between the current and previous measurements
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;
  
  /*
  std::cout << x_ << std::endl;
  std::cout << P_ << std::endl;
  std::cout << Xsig_pred_ << std::endl;
  std::cout << weights_ << std::endl;
  */
  // std::cout << delta_t << std::endl;
  
  // Prediction step
  // std::cout << "n_aug_ = " << n_aug_ << std::endl;
  // std::cout << "lambda_ = " << lambda_ << std::endl;
  UKF::Prediction(delta_t);

  // Correction step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    std::cout << "Process radar measurement\n";
    UKF::UpdateRadar(meas_package);

  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    std::cout << "Process lidar measurement\n";
    UKF::UpdateLidar(meas_package);

  }
    //  angle normalization
  // while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  // while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
   
}

void UKF::AugmentedSigmaPoints() {
  
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

    // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + sqrt(n_aug_ + lambda_)*A.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(n_aug_ + lambda_)*A.col(i);
  }

  // write result
  Xsig_pred_ = Xsig_aug;    // temp
}

void UKF::SigmaPointPrediction(double delta_t) {

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);

  for (int i=0; i<2*n_aug_ + 1; i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);
    double nu_a = Xsig_pred_(5,i);
    double nu_yawdd = Xsig_pred_(6,i);

    double px_p, py_p;
     if (fabs(yawd) > 0.001)
     {
         px_p = p_x + v/yawd*(sin(yaw+yawd*delta_t) -sin(yaw));
         py_p = p_y + v/yawd*(-cos(yaw+yawd*delta_t) + cos(yaw));
     } else{
         px_p = p_x + v*cos(yaw)*delta_t;
         py_p = p_y + v*sin(yaw)*delta_t;
     }
     double v_p = v;
     double yaw_p = yaw + yawd*delta_t;
     double yawd_p = yawd; 

     // add noise
     px_p = px_p + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
     py_p = py_p + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
     v_p = v_p + delta_t*nu_a;
     yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
     yawd_p = yawd_p + delta_t*nu_yawdd;

    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;

  }

  // write result
  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance() {
  
    // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // predict state mean
  x = Xsig_pred_*weights_;

  // predict state covariance matrix
  P.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {     
      VectorXd x_diff = Xsig_pred_.col(i) - x;

      // angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      P = P + weights_(i) * x_diff * x_diff.transpose();
  }

  // write result
  x_ = x;
  P_ = P;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  UKF::AugmentedSigmaPoints();
  UKF::SigmaPointPrediction(delta_t);
  UKF::PredictMeanAndCovariance();
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  // lidar measurement vector
  VectorXd z = meas_package.raw_measurements_;
  
  //measurement matrix
	MatrixXd H = MatrixXd(2, 5);
	H << 1, 0, 0, 0, 0,
			 0, 1, 0, 0, 0;

  MatrixXd R = MatrixXd(2,2);
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  VectorXd z_pred = H * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H.transpose();
	MatrixXd S = H * P_ * Ht + R;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;
  //std::cout << y << std::endl;
	//new estimate
	x_ = x_ + (K * y);

	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H) * P_;
}

void UKF::PredictRadarMeasurement(Eigen::MatrixXd* Zsig_out, Eigen::VectorXd* z_out, Eigen::MatrixXd* S_out){
  
  int n_z = 3;  // hardcode the dimension of radar measurement

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      // double ywad = Xsig_pred(4,i);
      
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig(1,i) = atan2(p_y, p_x);
      Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/Zsig(0,i);
      
  }
  
  // calculate mean predicted measurement
  z_pred = Zsig*weights_;
  
  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {     
      VectorXd z_diff = Zsig.col(i) - z_pred;

      // angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  S =  S + R;

    // write result
  *Zsig_out = Zsig; 
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  VectorXd z = meas_package.raw_measurements_;
  
  int n_z = 3;

  // create matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  // create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);

  UKF::PredictRadarMeasurement(&Zsig, &z_pred, &S);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {     
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI){
      std::cout << "Norm\n";
      z_diff(1)-=2.*M_PI;
    } 
    while (z_diff(1)<-M_PI){
      std::cout << "Norm\n";
      z_diff(1)+=2.*M_PI;
    }

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
   while (z_diff(1)>= M_PI) {
     z_diff(1)-=2.*M_PI;
  }

  while (z_diff(1)<-M_PI) {
     z_diff(1)+=2.*M_PI;
   }

  // update state mean and covariance matrix
  x_ = x_ + K*z_diff;

  P_ = P_ - K*S*K.transpose();

}