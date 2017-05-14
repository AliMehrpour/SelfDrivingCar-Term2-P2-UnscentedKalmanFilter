#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "data.h"
#include <vector>
#include <stdexcept>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::pair;

class UKF {
 public:
  bool is_initialized_; // Initially false, set to true in first call of ProcessMeasurement

  long previous_timestamp_; // Previous timestamp

  VectorXd x_; // State vector [px, py, vel_abs, yaw_angle, yaw_rate]
  MatrixXd P_; // State covariance matrix

  int n_laser_z_; // Laser measurement dimension
  int n_radar_z_; // Laser measurement dimension

  MatrixXd R_laser_; // Laser measurement covariance matrices
  MatrixXd R_radar_; // Radar measurement covariance matrices

  MatrixXd Xsig_pred_; // Predicted sigma matrix

  long time_us_; // Time when state is true

  double std_a_; // Process noise standard deviation longitudinal acceleration in m/s^2

  double std_yawdd_; // Process noise standard deviation yaw acceleration in rad/s^2

  double std_laspx_; // Laser measurement noise standard deviation position1 in m

  double std_laspy_; // Laser measurement noise standard deviation position2 in m

  double std_radr_; // Radar measurement noise standard deviation radius in m

  double std_radphi_; // Radar measurement noise standard deviation angle in rad

  double std_radrd_; // Radar measurement noise standard deviation radius change in m/s

  int n_x_; // State dimension

  int n_aug_; // Augmented state dimension

  double lambda_; // Sigma point spreading parameter

  double NIS_radar_; // Current NIS for radar
  double NIS_laser_; // Current NIS for laser

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param measurement The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(Data measurement);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param measurement The measurement at k+1
   */
  void UpdateLidar(Data measurement);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param measurement The measurement at k+1
   */
  void UpdateRadar(Data measurement);

 private:
  /*
   * Generate Sigma Points
   * @param n_x State Dimension
   * @param x State Vector
   * @param P State Coveriance Matrix
   */
  MatrixXd GenerateSigmaPoints(int n_x, const VectorXd& x, const MatrixXd& P);

  /*
   * Augment Sigma Points
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param std_a Process noise standard deviation longitudinal acceleration in m/s^2
   * @param std_yawdd Process noise standard deviation yaw acceleration in rad/s^2
   * @param x State vector
   * @param P State covariance Matrix
   */
  MatrixXd AugmentedSigmaPoints(int n_x, int n_aug, double std_a,
                                double std_yawdd, const VectorXd& x,
                                const MatrixXd& P);

  /*
   * Predict Sigma Points
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param delta_t Time between k and k+1 in s
   * @param Xsig_aug Augmented sigma points matrix
   */
  MatrixXd PredictSigmaPoints(int n_x, int n_aug, double delta_t,
                              const MatrixXd& Xsig_aug);

  /*
   * Predict Mean and Covariance
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param Xsig_pred Predicted sigma points matrix
   */
  pair<VectorXd, MatrixXd> PredictMeanAndCovariance(int n_x, int n_aug,
                                                     const MatrixXd& Xsig_pred);

  /*
   * Predict Radar Measurement
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param n_z Measurement dimension, radar can measure r, phi, and r_dot
   * @param std_radr Radar measurement noise standard deviation radius in m
   * @param std_radphi Radar measurement noise standard deviation angle in rad
   * @param std_radrd Radar measurement noise standard deviation radius change in m/s
   * @param Xsig_pred Predicted sigma points matrix
   */
  pair<VectorXd, MatrixXd> PredictRadarMeasurement(int n_x, int n_aug, int n_z,
                                                    double std_radr,
                                                    double std_radphi,
                                                    double std_radrd,
                                                    const MatrixXd& Xsig_pred);

  /*
   * Predict Laser Measurement
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param n_z Measurement dimension, laser can measure px, py
   * @param std_laspx Laser measurement noise standard deviation position1 in m
   * @param std_laspy Laser measurement noise standard deviation position2 in m
   * @param Xsig_pred Predicted sigma points matrix
   */
  pair<VectorXd, MatrixXd> PredictLaserMeasurement(int n_x, int n_aug, int n_z,
                                                    double std_laspx,
                                                    double std_laspy,
                                                    const MatrixXd& Xsig_pred);

  /*
   * Predict Measurement Covariance
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param n_z Measurement dimension, radar can measure r, phi, and r_dot
   * @param Zsig sigma points in measurement space MatrixXd(n_z, 2 * n_aug + 1)
   * @param Xsig_pred Predicted sigma points matrix
   */
  pair<VectorXd, MatrixXd> PredictMeasurementCovariance(
      int n_x, int n_aug, int n_z, const MatrixXd& Zsig,
      const MatrixXd& Xsig_pred);

  /*
   * Update State
   * @param n_x State dimension
   * @param n_aug Augmentation dimension
   * @param n_z Measurement dimension, radar can measure r, phi, and r_dot
   * @param x State vector
   * @param P State covariance Matrix
   * @param Zsig sigma points in measurement space MatrixXd(n_z, 2 * n_aug + 1)
   * @param z_pred vector for mean predicted measurement
   * @param S predicted measurement covariance
   * @param z vector for incoming measurement
   *
   */
  pair<VectorXd, MatrixXd> UpdateState(int n_x, int n_aug, int n_z,
                                        const MatrixXd& Xsig_pred,
                                        const VectorXd& x, const MatrixXd& P,
                                        const MatrixXd& Zsig,
                                        const VectorXd& z_pred,
                                        const MatrixXd& S, const VectorXd& z);
  /*
   * Sigma Points in Radar Measurement Space
   * @param n_z Measurement dimension, radar can measure r, phi, and r_dot
   * @param n_aug Augmentation dimension
   * @param Xsig_pred Predicted sigma points matrix
   */
  MatrixXd SigmaPointsRadarMeasurementSpace(int n_z, int n_aug,
                                            const MatrixXd& Xsig_pred);

  /*
   * Sigma Points in Laser Measurement Space
   * @param n_z Measurement dimension, lidar can measure px, py
   * @param n_aug Augmentation dimension
   * @param Xsig_pred Predicted sigma points matrix
   */
  MatrixXd SigmaPointsLaserMeasurementSpace(int n_z, int n_aug,
                                            const MatrixXd& Xsig_pred);

  VectorXd WeightsVector(int n_aug, double lambda);

  VectorXd MeanPredictedMeasurement(int n_z, int n_aug, const VectorXd& weights,
                                    const MatrixXd& Zsig);

  MatrixXd MeasurementCovarianceMatrixS(int n_z, int n_aug,
                                        const MatrixXd& Zsig,
                                        const VectorXd& z_pred,
                                        const VectorXd& weights);

  MatrixXd CrossCorrelationMatrix(int n_x, int n_aug, int n_z,
                                  const MatrixXd& Zsig, const VectorXd& z_pred,
                                  const MatrixXd& Xsig_pred, const VectorXd& x,
                                  const VectorXd& weights);

  static MatrixXd NoiseCovarianceMatrixRadar(int n_z, double std_radr,
                                             double std_radphi,
                                             double std_radrd);

  static MatrixXd NoiseCovarianceMatrixLaser(int n_z, double std_laspx,
                                             double std_laspy);
};

#endif /* UKF_H */
