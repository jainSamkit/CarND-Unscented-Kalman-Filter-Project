#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1,0,0,0,0,
            0,1,0,0,0,
            0,0,0.8,0,0,
            0,0,0,1,0,
            0,0,0,0,1;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/8.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:


  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_= false;
  
  time_us_=0;
  
  n_aug=7;
  n_x_=5;

  lambda_pred = 3 - n_x_;
  lambda_aug= 3 - n_aug;


  weights_aug=VectorXd(2*n_aug+1);

  double weight_ = lambda_aug/(lambda_aug+n_aug);
  weights_aug(0)=weight_;

  for(int i=1;i<2*n_aug+1;i++)
  {
    weight_= 0.5/(n_aug+lambda_aug);
    weights_aug(i) = weight_;
  } 

  weight_=lambda_pred/(lambda_pred+n_x_);
  weights_pred=VectorXd(2*n_x_+1);
  weights_pred(0) = weight_;

  for(int i=1;i<2*n_x_+1;i++)
  {
    weight_= 0.5/(n_x_+lambda_pred);
    weights_pred(i) = weight_;
  }   
  X_sig_pred = MatrixXd(n_x_,2*n_aug+1);
  X_sig_pred.fill(0.0);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_)
  {
    if(meas_package.sensor_type_== MeasurementPackage::RADAR)
    {
      float rho=meas_package.raw_measurements_[0];
      float phi=meas_package.raw_measurements_[1];
      while(phi>M_PI) phi-=2*M_PI;
      while(phi<-M_PI) phi+=2*M_PI;
      float rho_dot=meas_package.raw_measurements_[2];

      x_<<rho*cos(phi),rho*sin(phi),0,0,0;
      cout<<"samkit"<<endl;

    }
    else if(meas_package.sensor_type_== MeasurementPackage::LASER)
    {
      x_<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;
      
    }
    if(fabs(x_[0])<0.0001 || fabs(x_[1])<0.0001)
    {
      if(fabs(x_[0])<0.0001)
      {
        x_[0]=0.1;
      }
      else
        x_[1]=0.1;
    }
    is_initialized_=true;
    time_us_ = meas_package.timestamp_;
    
    return;
  }
  double delta_t = (meas_package.timestamp_-time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;
  
  if(meas_package.sensor_type_== MeasurementPackage::LASER)
  {
    use_radar_=false;
    use_laser_=true;

    
 
  }
  else
  {
    use_radar_=true;
    use_laser_=false;

    
  }
  // cout<<meas_package.raw_measurements_<<endl;
  Prediction(delta_t);
  cout<<"After Prediction"<<endl;
  cout<<x_<<endl;
  cout<<endl;
  cout<<P_<<endl;

  cout<<"Printing meas_package"<<endl;
  cout<<meas_package.raw_measurements_<<endl;
  cout<<endl;
  if(use_laser_)
  { 
    UpdateLidar(meas_package);
    cout<<"After Lidar"<<endl;
    cout<<x_<<endl;
    cout<<P_<<endl;


  }
  else
  {
    UpdateRadar(meas_package);
    cout<<"After Radar"<<endl;
    cout<<x_<<endl;
    cout<<endl;
    cout<<P_<<endl;
  }

  
  // var_v.push_back(delta_vel*delta_vel);
  // var_yaw_rate.push_back(delta_yaw_rate*delta_yaw_rate);
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  /*generating sigma points first*/
  MatrixXd A = P_.llt().matrixL();
  MatrixXd X_sig = MatrixXd(n_x_,2*n_x_+1);

  X_sig.col(0)=x_;

  for(int i=0;i<n_x_;i++)
  {
    X_sig.col(i+1) = x_ + sqrt(lambda_pred+n_x_) * A.col(i);
    X_sig.col(i+n_x_+1) = x_ -sqrt(lambda_pred+n_x_) * A.col(i);
  }

  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7,7);

  MatrixXd X_sig_aug = MatrixXd(n_aug,2*n_aug+1);

  x_aug.head(5)=x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;


  MatrixXd L = P_aug.llt().matrixL();

  X_sig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug; i++)
  {
    X_sig_aug.col(i+1)       = x_aug + sqrt(lambda_aug+n_aug) * L.col(i);
    X_sig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_aug+n_aug) * L.col(i);
  }
  
  for (int i = 0; i< 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = X_sig_aug(0,i);
    double p_y = X_sig_aug(1,i);
    double v = X_sig_aug(2,i);
    double yaw = X_sig_aug(3,i);
    double yawd = X_sig_aug(4,i);
    double nu_a = X_sig_aug(5,i);
    double nu_yawdd = X_sig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    X_sig_pred(0,i) = px_p;
    X_sig_pred(1,i) = py_p;
    X_sig_pred(2,i) = v_p;
    X_sig_pred(3,i) = yaw_p;
    X_sig_pred(4,i) = yawd_p;
  }

  

  /*predicting mean state and covariance*/

  x_.fill(0.0);
  for(int i=0;i<2*n_aug+1;i++)
  {
    x_+=weights_aug(i) * X_sig_pred.col(i);
  }
  // cout<<"Here I am printing x_"<<endl;
  // cout<<x_<<endl;

  P_.fill(0.0);
  // cout<<"printing P_"<<endl;
  
  for(int i=0;i<2*n_aug+1;i++)
  {
    VectorXd x_diff = X_sig_pred.col(i)-x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2*M_PI;

    P_+=weights_aug(i) * x_diff * x_diff.transpose();
    
  }
  // cout<<P_<<endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z=2;
  MatrixXd H_laser = MatrixXd(2,5);
  H_laser << 1,0,0,0,0,
             0,1,0,0,0;


  VectorXd z_pred = H_laser * x_;
  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1];

  VectorXd y = z-z_pred;
  MatrixXd Ht = H_laser.transpose();


  MatrixXd R_laser = MatrixXd(n_z,n_z);
  R_laser <<std_laspx_*std_laspx_,0,
      0,std_laspy_*std_laspy_;

  MatrixXd S = H_laser * P_ * Ht + R_laser;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  x_=x_+(K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K*H_laser) * P_;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  

  int n_z = 3;
  MatrixXd Z_sigma_pred = MatrixXd(n_z,2*n_aug+1);
  for(int i=0;i<2*n_aug+1;i++)
  {
    double px = X_sig_pred(0,i);
    double py = X_sig_pred(1,i);
    double v = X_sig_pred(2,i);
    double yaw = X_sig_pred(3,i);

    double v1 = v*cos(yaw);
    double v2 = v*sin(yaw);

    Z_sigma_pred(0,i) = sqrt(px*px + py*py);
    Z_sigma_pred(1,i) = atan2(py,px);
    Z_sigma_pred(2,i) = (px*v1 + py*v2)/sqrt(px*px + py*py);
  }

  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0;i<2*n_aug+1;i++)
  {
    z_pred+=weights_aug(i) * Z_sigma_pred.col(i);
  }

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],meas_package.raw_measurements_[2];

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Z_sigma_pred.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;

    S = S + weights_aug(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R_radar = MatrixXd(n_z,n_z);
  R_radar.fill(0.0);
  R_radar(0,0) = std_radr_*std_radr_;
  R_radar(1,1) = std_radphi_*std_radphi_;
  R_radar(2,2) = std_radrd_*std_radrd_;

  S = S + R_radar;

  MatrixXd Tc = MatrixXd(n_x_,n_z);
  Tc.fill(0.0);

  for(int i=0;i<2*n_aug+1;i++)
  {
    VectorXd x_diff = X_sig_pred.col(i)-x_;
    while(x_diff(3)>M_PI) x_diff(3)-=2*M_PI;
    while(x_diff(3)<-M_PI) x_diff(3)+=2*M_PI;

    // cout<<"Printing predicted sigma points"<<endl;
    // cout<<X_sig_pred<<endl;
    // cout<<endl;

    // cout<<"Printing x_diff"<<endl;
    // cout<<x_diff<<endl;

    // cout<<"Printing x_"<<endl;
    // cout<<x_<<endl;

    VectorXd z_diff = Z_sigma_pred.col(i)-z_pred;
    while(z_diff(1)>M_PI) z_diff(1)-=2*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;

    Tc+=weights_aug(i)*x_diff * z_diff.transpose();

  }
  // cout<<"Printing Tc"<<endl;
  // cout<<Tc<<endl;

  // cout<<"Printing S inverse"<<endl;
  // cout<<S.inverse()<<endl;

  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2*M_PI;

  //update state mean and covariance matrix

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
