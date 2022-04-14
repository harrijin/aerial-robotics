// Toolbox of elementary navigation-related functions.  Matrix operations
// employ Eigen matrix objects.

#ifndef __NAVTOOLBOX_EIGEN_H
#define __NAVTOOLBOX_EIGEN_H

#include <Eigen/Dense>
#include "definitions.h"

typedef Eigen::Matrix<f64,4,1> Vector4d;

namespace navtbx {
  Vector4d dc2quat(const Eigen::Matrix3d& dc); 
  Eigen::Matrix3d quat2dc(const Vector4d& q);
  Eigen::Matrix3d euler2dc(const Eigen::Vector3d& e);
  Eigen::Vector3d dc2euler(const Eigen::Matrix3d& RBW);
  Vector4d euler2quat(const Eigen::Vector3d& e);
  Eigen::Vector3d quat2euler(const Vector4d& q);
  f64 heading(const Vector4d& q);
  f64 rotationAngle(const Eigen::Matrix3d& dc);
  Eigen::Matrix3d wahbaSolver(const Eigen::VectorXd& aVec,
                              const Eigen::MatrixXd& vWMat,
                              const Eigen::MatrixXd& vBMat);
  Eigen::Vector3d elAzToUnitEnu(f64 el, f64 az);
  ReturnValue enuVecToElAz(const Eigen::Vector3d& v, f64& el, f64& az);
  ReturnValue vecToElAz(const Eigen::Vector3d& v, f64& el, f64& az);
  ReturnValue ecef2LatLonAlt(const Eigen::Vector3d& pVec,
                             f64& lat, f64& lon, f64& alt);
  ReturnValue latLonAlt2Ecef(const f64 lat, const f64 lon, const f64 alt,
                             Eigen::Vector3d& pVec);
  ReturnValue getElAz(const Eigen::Vector3d& rRx, const Eigen::Vector3d& rSv,
                      f64& el, f64& az);
  ReturnValue getElAz(const Eigen::Vector3d& rRx,const Eigen::Vector3d& rSv,
                      const f64 lat, const f64 lon, const f64 alt,
                      f64& el, f64& az);
  Eigen::Matrix3d ecef2enu(const Eigen::Vector3d& rEcef);
  Vector4d qmult(const Vector4d& q1, const Vector4d& q2);
  Vector4d diffQuat2FullQuat(const Eigen::Vector3d& de);
  Eigen::Vector3d getGravitationalAccelerationVector(const Eigen::Vector3d& rEcef);
  Eigen::Matrix3d Rx(const f64 theta);
  Eigen::Matrix3d Ry(const f64 theta);
  Eigen::Matrix3d Rz(const f64 theta);
  Eigen::Matrix3d cpe(const Eigen::Vector3d& vec);
} // namespace

#endif
