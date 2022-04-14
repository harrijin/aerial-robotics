#include "structurecomputer.h"
#include <Eigen/LU>

void pr(Eigen::MatrixXd m) { std::cout << m << std::endl; }
void pr(Eigen::VectorXd m) { std::cout << m << std::endl; }
void pr(Eigen::Matrix3d m) { std::cout << m << std::endl; }
void pr(Eigen::Vector3d m) { std::cout << m << std::endl; }
void pr(Eigen::Vector2d m) { std::cout << m << std::endl; }

Eigen::Vector2d backProject(const Eigen::Matrix3d& RCI,
                            const Eigen::Vector3d& rc_I,
                            const Eigen::Vector3d& X3d) {
  using namespace Eigen;
  Vector3d t = -RCI * rc_I;
  MatrixXd Pextrinsic(3, 4);
  Pextrinsic << RCI, t;
  SensorParams sp;
  MatrixXd Pc = sp.K() * Pextrinsic;
  VectorXd X(4, 1);
  X.head(3) = X3d;
  X(3) = 1;
  Vector3d x = Pc * X;
  Vector2d xc_pixels = (x.head(2) / x(2)) / sp.pixelSize();
  return xc_pixels;
}

Eigen::Vector3d pixelsToUnitVector_C(const Eigen::Vector2d& rPixels) {
  using namespace Eigen;
  SensorParams sp;
  // Convert input vector to meters
  Vector2d rMeters = rPixels * sp.pixelSize();
  // Write as a homogeneous vector, with a 1 in 3rd element
  Vector3d rHomogeneous;
  rHomogeneous.head(2) = rMeters;
  rHomogeneous(2) = 1;
  // Invert the projection operation through the camera intrinsic matrix K to
  // yield a vector rC in the camera coordinate frame that has a Z value of 1
  Vector3d rC = sp.K().lu().solve(rHomogeneous);
  // Normalize rC so that output is a unit vector
  return rC.normalized();
}

void StructureComputer::clear() {
  // Zero out contents of point_
  point_.rXIHat.fill(0);
  point_.Px.fill(0);
  // Clear bundleVec_
  bundleVec_.clear();
}

void StructureComputer::push(std::shared_ptr<const CameraBundle> bundle) {
  bundleVec_.push_back(bundle);
}
Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}
// This function is where the computation is performed to estimate the
// contents of point_.  The function returns a copy of point_.
//
Point StructureComputer::computeStructure() {
  // Throw an error if there are fewer than 2 CameraBundles in bundleVec_,
  // since in this case structure computation is not possible.
  if (bundleVec_.size() < 2) {
    throw std::runtime_error(
        "At least 2 CameraBundle objects are "
        "needed for structure computation.");
  }

  // *********************************************************************
  // Fill in here the required steps to calculate the 3D position of the
  // feature point and its covariance.  Put these respectively in
  // point_.rXIHat and point_.Px
  // *********************************************************************
  Eigen::MatrixXd H(2*bundleVec_.size(), 4);
  // Populate H matrix;
  for(auto &bundle: bundleVec_){
    Eigen::MatrixXd proj_mat(3,4);
    Eigen::MatrixXd intrinsic_mat(3,4);
    intrinsic_mat << bundle->RCI, -1*bundle->RCI*bundle->rc_I;
    proj_mat = sensorParams_.K()*intrinsic_mat;
    double x_tilde = bundle->rx[0]*sensorParams_.pixelSize();
    double y_tilde = bundle->rx[1]*sensorParams_.pixelSize();
    H << x_tilde*proj_mat.block<1,4>(2,0) - proj_mat.block<1,4>(0,0);
    H << y_tilde*proj_mat.block<1,4>(2,0) - proj_mat.block<1,4>(1,0);
  }
  Eigen::MatrixXd Hr(2*bundleVec_.size(), 3);
  Eigen::MatrixXd z(2*bundleVec_.size(), 1);
  Hr << H.block(0,0, 2*bundleVec_.size(),3);
  z << H.block(0,3, 2*bundleVec_.size(),1);
  Eigen::MatrixXd R(2*bundleVec_.size(), 2*bundleVec_.size());
  R << pow(sensorParams_.pixelSize(),2)*blkdiag(sensorParams_.Rc(),bundleVec_.size());
  Eigen::MatrixXd Rinv = R.inverse();
  Eigen::MatrixXd Px_inv = Hr.transpose()*Rinv*Hr;
  point_.Px = Px_inv.inverse();
  point_.rXIHat = Px_inv.ldlt().solve(Hr.transpose());
  return point_;
}
