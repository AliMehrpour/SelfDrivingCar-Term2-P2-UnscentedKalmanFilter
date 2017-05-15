#include <iostream>
#include "ukf.h"
#include "tools.h"
#include "templates.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::pair;
using std::make_pair;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // Measurement dimensions
  n_laser_z_ = 2;
  n_radar_z_ = 3;

  // State Dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // initial sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // augmented state dimension
  n_aug_ = 7;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.45;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55; //0.1;

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

  // measurement covariance matrices
  R_laser_ = NoiseCovarianceMatrixLaser(n_laser_z_, std_laspx_, std_laspy_);
  R_radar_ = NoiseCovarianceMatrixRadar(n_radar_z_, std_radr_, std_radphi_, std_radrd_);
}

UKF::~UKF() { }

/**
 * @param {Data} measurement The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(Data measurement) {
  if (!is_initialized_) {
    cout << "Initializing Fusion Unscented Kalman Filter" << endl;

    float px, py, v, yaw, yawd;
    if (measurement.sensor_type_ == Data::LASER) {
      px = measurement.data_[0];
      py = measurement.data_[1];
      v = 0.0;
      yaw = 0.0;
      yawd = 0.0;
    }
    else if (measurement.sensor_type_ == Data::RADAR) {
      // Convert polar to cartesian coordinates
      float rho = measurement.data_[0];       // Range - radial distance from origin
      float phi = measurement.data_[1];       // Bearing - angle between rho and x
      float rho_dot = measurement.data_[2];   // Radial Velocity - change of p (range rate)

      px = rho * cos(phi); // meteres
      py = rho * sin(phi); // meteres
      v = rho_dot;         // meteres/sec
      yaw = phi;           // radians
      yawd = 0.0;          // radians/sec
    }

    // Handling special cases initialization issues
    if (fabs(px) < EPS and fabs(py) < EPS) {
      px = EPS;
      py = EPS;
    }

    // Initialize x_ with first measurement
    x_ = VectorXd(5);
    x_ << px, py, v, yaw, yawd;

    // Initial covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    is_initialized_ = true;
    previous_timestamp_ = measurement.timestamp_;

    return;
  }

  // Calculate delta time and update previous timestamp
  double delta_t = (measurement.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement.timestamp_;

  // Prediction
  Prediction(delta_t);

  // Update
  switch (measurement.sensor_type_) {
    case Data::RADAR:
      UpdateRadar(measurement);
      break;
    case Data::LASER:
      UpdateLidar(measurement);
      break;
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug = AugmentedSigmaPoints(n_x_, n_aug_, std_a_, std_yawdd_, x_, P_);

  Xsig_pred_ = PredictSigmaPoints(n_x_, n_aug_, delta_t, Xsig_aug);

  pair<VectorXd, MatrixXd> pred = PredictMeanAndCovariance(n_x_, n_aug_, Xsig_pred_);
  x_ = pred.first;
  P_ = pred.second;
}

/**
 * Generate Sigma Points
 */
MatrixXd UKF::GenerateSigmaPoints(int n_x, const VectorXd& x, const MatrixXd& P) {
  // calculate spreading parameter
  double lambda = 3 - n_x;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  //calculate square root of P
  MatrixXd A = P.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0) = x;

  //set remaining sigma points
  for (int i = 0; i < n_x; i++) {
    Xsig.col(i + 1) = x + sqrt(lambda + n_x) * A.col(i);
    Xsig.col(i + 1 + n_x) = x - sqrt(lambda + n_x) * A.col(i);
  }

  return Xsig;
}

/**
 * Augment Sigma Points
 */
MatrixXd UKF::AugmentedSigmaPoints(int n_x, int n_aug, double std_a,
                                   double std_yawdd, const VectorXd& x,
                                   const MatrixXd& P) {
  //define spreading parameter
  double lambda = 3 - n_aug;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P;
  P_aug(5, 5) = std_a * std_a;
  P_aug(6, 6) = std_yawdd * std_yawdd;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i);
  }

  return Xsig_aug;
}

/**
 * Predict Sigma Points
 */
MatrixXd UKF::PredictSigmaPoints(int n_x, int n_aug, double delta_t,
                                 const MatrixXd& Xsig_aug) {
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //predict sigma points
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }

  return Xsig_pred;
}

/**
 * Predict Mean and Covariance
 */
pair<VectorXd, MatrixXd> UKF::PredictMeanAndCovariance(
    int n_x, int n_aug, const MatrixXd& Xsig_pred) {
  //define spreading parameter
  double lambda = 3 - n_aug;

  //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

  //create vector for weights
  VectorXd weights = VectorXd(2 * n_aug + 1);

  // set weights
  double weight_0 = lambda / (lambda + n_aug);
  weights(0) = weight_0;
  for (int i = 1; i < 2 * n_aug + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug + lambda);
    weights(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
    x = x + weights(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 1; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - Xsig_pred.col(0);

    //angle normalization
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));

    P = P + weights(i) * x_diff * x_diff.transpose();
  }

  return std::make_pair(x, P);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(Data measurement) {

  /*
   * Predict Measurement
   */
  // Measurement dimension, laser can measure px, py
  int n_z = n_laser_z_;
  VectorXd z = VectorXd(n_z);
  z << measurement.data_[0], measurement.data_[1];

  // predict measurement
  VectorXd z_pred;
  MatrixXd S;
  pair<VectorXd, MatrixXd> pred = PredictLaserMeasurement(n_x_, n_aug_, n_z, std_laspx_, std_laspy_, Xsig_pred_);
  z_pred = pred.first;
  S = pred.second;

  /*
   * Update State
   */
  MatrixXd Zsig = SigmaPointsLaserMeasurementSpace(n_z, n_aug_, Xsig_pred_);

  pair<VectorXd, MatrixXd> state = UpdateState(n_x_, n_aug_, n_z, Xsig_pred_, x_, P_, Zsig, z_pred, S, z);
  x_ = state.first;
  P_ = state.second;

  // update NIS
  NIS_laser_ = (measurement.data_-z_pred).transpose()*S.inverse()*(measurement.data_-z_pred);
}

MatrixXd UKF::SigmaPointsLaserMeasurementSpace(int n_z, int n_aug,
                                               const MatrixXd& Xsig_pred) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    // extract values for better readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }
  return Zsig;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(Data measurement) {

  // Measurement dimension, radar can measure r, phi, and r_dot
  int n_z = n_radar_z_;

  VectorXd z = VectorXd(n_z);
  z << measurement.data_[0], measurement.data_[1], measurement.data_[2];

  // predict measurement
  VectorXd z_pred;
  MatrixXd S;
  pair<VectorXd, MatrixXd> pred = PredictRadarMeasurement(n_x_, n_aug_, n_z, std_radr_,
                                                          std_radphi_, std_radrd_, Xsig_pred_);
  z_pred = pred.first;
  S = pred.second;

  MatrixXd Zsig = SigmaPointsRadarMeasurementSpace(n_z, n_aug_, Xsig_pred_);

  pair<VectorXd, MatrixXd> state = UpdateState(n_x_, n_aug_, n_z, Xsig_pred_, x_, P_, Zsig, z_pred, S, z);
  x_ = state.first;
  P_ = state.second;

  NIS_radar_ = (measurement.data_-z_pred).transpose()*S.inverse()*(measurement.data_-z_pred);
}

MatrixXd UKF::SigmaPointsRadarMeasurementSpace(int n_z, int n_aug,
                                               const MatrixXd& Xsig_pred) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    // extract values for better readability
    double p_x = Xsig_pred(0, i);
    double p_y = Xsig_pred(1, i);
    double v = Xsig_pred(2, i);
    double yaw = Xsig_pred(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    // measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);           //r
    Zsig(1, i) = atan2(p_y, p_x);                       //phi
    if (Zsig(0, i) != 0) {
      Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i);  //r_dot
    }
    else {
      Zsig(2, i) = 0;
    }

  }
  return Zsig;
}

MatrixXd UKF::NoiseCovarianceMatrixRadar(int n_z, double std_radr,
                                         double std_radphi, double std_radrd) {
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr * std_radr, 0, 0,
       0, std_radphi * std_radphi, 0,
       0, 0, std_radrd * std_radrd;
  return R;
}

/*
 * Predict Radar Measurement
 */
pair<VectorXd, MatrixXd> UKF::PredictRadarMeasurement(
    int n_x, int n_aug, int n_z, double std_radr, double std_radphi,
    double std_radrd, const MatrixXd& Xsig_pred) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = SigmaPointsRadarMeasurementSpace(n_z, n_aug, Xsig_pred);

  MatrixXd S;
  VectorXd z_pred;
  pair<VectorXd, MatrixXd> pred = PredictMeasurementCovariance(n_x, n_aug, n_z, Zsig, Xsig_pred);
  z_pred = pred.first;
  S = pred.second;

  //add measurement noise covariance matrix
  MatrixXd R = NoiseCovarianceMatrixRadar(n_z, std_radr, std_radphi, std_radrd);
  S = S + R;

  return std::make_pair(z_pred, S);
}

MatrixXd UKF::NoiseCovarianceMatrixLaser(int n_z, double std_laspx,
                                         double std_laspy) {
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx * std_laspx, 0, 0, std_laspy * std_laspy;
  return R;
}

/*
 * Predict Laser Measurement
 */
pair<VectorXd, MatrixXd> UKF::PredictLaserMeasurement(
    int n_x, int n_aug, int n_z, double std_laspx, double std_laspy,
    const MatrixXd& Xsig_pred) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = SigmaPointsLaserMeasurementSpace(n_z, n_aug, Xsig_pred);

  MatrixXd S;
  VectorXd z_pred;
  pair<VectorXd, MatrixXd> pred = PredictMeasurementCovariance(n_x, n_aug, n_z, Zsig, Xsig_pred);
  z_pred = pred.first;
  S = pred.second;

  //add measurement noise covariance matrix
  MatrixXd R = NoiseCovarianceMatrixLaser(n_z, std_laspx, std_laspy);
  S = S + R;

  return std::make_pair(z_pred, S);
}

VectorXd UKF::MeanPredictedMeasurement(int n_z, int n_aug,
                                       const VectorXd& weights,
                                       const MatrixXd& Zsig) {
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }
  return z_pred;
}

MatrixXd UKF::MeasurementCovarianceMatrixS(int n_z, int n_aug,
                                           const MatrixXd& Zsig,
                                           const VectorXd& z_pred,
                                           const VectorXd& weights) {
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 1; i < 2 * n_aug + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - Zsig.col(0);

    //angle normalization
    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

    S = S + weights(i) * z_diff * z_diff.transpose();
  }
  return S;
}

/*
 * Predict Laser Measurement
 */
pair<VectorXd, MatrixXd> UKF::PredictMeasurementCovariance(
    int n_x, int n_aug, int n_z, const MatrixXd& Zsig,
    const MatrixXd& Xsig_pred) {

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = WeightsVector(n_aug, lambda);

  //mean predicted measurement
  VectorXd z_pred = MeanPredictedMeasurement(n_z, n_aug, weights, Zsig);

  //measurement covariance matrix S
  MatrixXd S = MeasurementCovarianceMatrixS(n_z, n_aug, Zsig, z_pred, weights);

  return make_pair(z_pred, S);

}

VectorXd UKF::WeightsVector(int n_aug, double lambda) {
  //set vector for weights
  VectorXd weights = VectorXd(2 * n_aug + 1);
  double weight_0 = lambda / (lambda + n_aug);
  weights(0) = weight_0;
  for (int i = 1; i < 2 * n_aug + 1; i++) {
    //2n+1 weights
    double weight = 0.5 / (n_aug + lambda);
    weights(i) = weight;
  }

  return weights;
}

MatrixXd UKF::CrossCorrelationMatrix(int n_x, int n_aug, int n_z,
                                     const MatrixXd& Zsig,
                                     const VectorXd& z_pred,
                                     const MatrixXd& Xsig_pred,
                                     const VectorXd& x,
                                     const VectorXd& weights) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 1; i < 2 * n_aug + 1; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - Zsig.col(0);

    //angle normalization
    z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - Xsig_pred.col(0);

    //angle normalization
    x_diff(3) = normalizeRadiansPiToMinusPi(x_diff(3));
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  return Tc;
}

/*
 * Update State
 */
pair<VectorXd, MatrixXd> UKF::UpdateState(int n_x, int n_aug, int n_z,
                                           const MatrixXd& Xsig_pred,
                                           const VectorXd& x, const MatrixXd& P,
                                           const MatrixXd& Zsig,
                                           const VectorXd& z_pred,
                                           const MatrixXd& S,
                                           const VectorXd& z) {

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = WeightsVector(n_aug, lambda);

  //calculate cross correlation matrix
  MatrixXd Tc = CrossCorrelationMatrix(n_x, n_aug, n_z, Zsig, z_pred, Xsig_pred,
                                       x, weights);

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = normalizeRadiansPiToMinusPi(z_diff(1));

  //update state mean and covariance matrix
  VectorXd x_update = x + K * z_diff;
  MatrixXd P_update = P - K * S * K.transpose();

  return make_pair(x_update, P_update);
}
