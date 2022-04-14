#include <cmath>
#include "navtoolbox.h"

namespace navtbx {

  // Converts elevation and azimuth angles (in radians) to a unit vector in ENU
  // coordinates.  Azimuth is defined to be zero along the north direction.
  //
  Eigen::Vector3d elAzToUnitEnu(f64 el, f64 az) {
    Eigen::Vector3d u;
    const f64 cEN = std::cos(el);
    u(0) = cEN*std::sin(az);
    u(1) = cEN*std::cos(az);
    u(2) = std::sin(el);
    return u;
  }

  // Returns the elevation and azimuth angles, in radians, corresponding to
  // the input vector.  The vector is assumed to be expressed in the ENU
  // coordinate system, whose X axis corresponds to an azimuth angle of 90 and
  // whose Z axis corresponds to an elevation angle of pi/2.  Azimuth is
  // measured positive clockwise from the Y (north) axis in the XY plane so
  // that the Y axis corresponds to a azimuth angle of 0.
  //
  ReturnValue enuVecToElAz(const Eigen::Vector3d& v, f64& el, f64& az) {
    Eigen::Vector3d u = v.normalized();
    el = std::asin(u(2));
    az = std::atan2(u(0), u(1));
    return RETVAL_SUCCESS;
  }

  // Returns the elevation and azimuth angles, in radians, corresponding to
  // the input vector.  The vector is assumed to be expressed in a
  // right-handed coordinate system whose X axis corresponds to an azimuth
  // angle of 0 and whose Z axis corresponds to an elevation angle of pi/2.
  // Azimuth is measured positive clockwise from the X axis in the XY plane so
  // that the Y axis corresponds to a azimuth angle of -pi/2 = 3*pi/2.
  //
  ReturnValue vecToElAz(const Eigen::Vector3d& v, f64& el, f64& az) {
    Eigen::Vector3d u = v.normalized();
    el = std::asin(u(2));
    az = std::atan2(-u(1), u(0));
    return RETVAL_SUCCESS;
  }

  // Converts ECEF position to WGS-84 latitude, longitude, and altitude.
  //
  // pVec is the 3-by-1 ECEF position in meters.  lat and lon are the output
  // latitude (geodetic) and longitude, in radians. alt is the output
  // altitude, in meters.  The conversion algorithm is derived from the
  // function Geodetic::Geodetic(const ECEF& right, GeoidModel* geo) found in
  // Geodetic.cpp in the GPSTk.  Returns RETVAL_SUCCESS if the lat and alt
  // values converge within maxIter.  Otherwise returns RETVAL_FAILURE.
  //
  ReturnValue ecef2LatLonAlt(const Eigen::Vector3d& pVec, f64& lat, f64& lon, f64& alt) {

    //----- Local constants
    const f64 aa = AA_WGS84;
    const f64 e2 = esquare_WGS84;
    const s32 maxIter = 5;
    const f64 convergenceFactor = 1.0e-9;
    const f64 convergenceFactor_x_aa = convergenceFactor*aa;

    //----- Convert to lat lon alt
    f64 X = pVec(0), Y = pVec(1), Z = pVec(2);
    f64 p = sqrt(X*X + Y*Y);
    lat = atan2(Z, p * (1.0 - e2));
    f64 slat, N, altOld, latOld;
    alt = 0.0;
    ReturnValue retVal = RETVAL_FAILURE;
    for(s32 ii=0;ii<maxIter;ii++) {
      slat = sin(lat);
      N = aa / sqrt(1.0 - e2 * slat * slat);
      altOld = alt;
      alt = p/cos(lat) - N;
      latOld = lat;
      lat = atan2(Z, p * (1.0 - e2 * (N/(N+alt))));
      if(fabs(lat-latOld) < convergenceFactor &&
         fabs(alt-altOld) < convergenceFactor_x_aa) {
        retVal = RETVAL_SUCCESS;
        break;
      }
    }
    lon = atan2(Y,X);
    if(lon < 0.0)
      lon += TWOPI;

    return retVal;
  }

  // Converts WGS-84 latitude, longitude, and altitude to ECEF position.
  //
  // lat and lon are the input latitude and longitude, in radians. alt is the
  // input altitude, in meters. pVec is the 3-by-1 ECEF position in meters.
  //
  ReturnValue latLonAlt2Ecef(const f64 lat, const f64 lon, const f64 alt,
                             Eigen::Vector3d& pVec) {

    //----- Local constants
    const f64 aa = AA_WGS84;
    const f64 bb = BB_WGS84;
    const f64 ee = e_WGS84;

    f64 N = aa / sqrt(1.0 - ee*ee*pow(sin(lat),2));
    pVec(0) = (N + alt)*cos(lat)*cos(lon);
    pVec(1) = (N + alt)*cos(lat)*sin(lon);
    pVec(2) = (bb*bb/(aa*aa)*N + alt)*sin(lat);

    return RETVAL_SUCCESS;
  }


  // Gets the elevation and azimuth angles el and az (in radians) of the SV
  // located at rSv with respect to the observer located at rRx.  3-by-1 vectors
  // rSv and rRx are expressed in meters in the ECEF reference frame.
  //
  ReturnValue getElAz(const Eigen::Vector3d& rRx, const Eigen::Vector3d& rSv,
                      f64& el, f64& az) {
    f64 lat, lon, alt;
    ecef2LatLonAlt(rRx,lat,lon,alt);
    return getElAz(rRx,rSv,lat,lon,alt,el,az);
  }

  // Gets the elevation and azimuth angles el and az (in radians) of the SV
  // located at rSv with respect to the observer located at rRx.  The 3-by-1
  // vectors rRx and rSv are expressed in meters in the ECEF reference frame.
  // The inputs lat (rad geodetic), lon (rad), and alt (meters) are assumed to
  // be the approximate latitude, longitude, and altitude (LLA) corresponding
  // to rRx.  Allowing these to be input in addition to rRx permits el and az
  // to be calculated rapidly in cases where both rRx and LLA have been
  // previously calculated.
  //
  ReturnValue getElAz(const Eigen::Vector3d& rRx,const Eigen::Vector3d& rSv,
                      const f64 lat, const f64 lon, const f64 alt,
                      f64& el, f64& az) {

    //----- Unit vectors of E, N, U local directions, in ECEF coordinates
    Eigen::Vector3d eECEF, nECEF, uECEF;
    f64 sinlat = sin(lat);
    f64 sinlon = sin(lon);
    f64 coslat = cos(lat);
    f64 coslon = cos(lon);
    eECEF(0) = -sinlon;
    eECEF(1) = coslon;
    eECEF(2) = 0.;
    nECEF(0) = -coslon*sinlat;
    nECEF(1) = -sinlon*sinlat;
    nECEF(2) = coslat;
    uECEF(0) = coslon*coslat;
    uECEF(1) = sinlon*coslat;
    uECEF(2) = sinlat;

    //----- Unit vector to satellite
    Eigen::Vector3d rhohat = rSv - rRx;
    rhohat.normalize();

    //----- Elevation and Azimuth angle calculations
    el = PI/2 - acos(uECEF.dot(rhohat));
    az = atan2(eECEF.dot(rhohat),nECEF.dot(rhohat));
    if(az < 0)
      az = az + TWOPI;

    return RETVAL_SUCCESS;
  }

  // Produces the rotation matrix Renu_ecef used to express a vector written
  // in ECEF coordinates relative to the ECEF vector rEcef as a vector written
  // in the local vertical, east, north, up (ENU) coordinate frame with origin
  // at rEcef; e.g., vEnu = Renu_ecef*vEcef.  The transformation is based on
  // the WGS-84 ellipsoid.
  //
  // Notes: See Misra and Enge appendix to chapter 4.
  //
  // INPUTS
  //
  // rECEF ------- 3-by-1 position vector locating the origin of the local ENU
  //               frame, in ECEF meters.
  //
  // OUTPUTS
  //
  // Renu_ecef --- 3-by-3 rotation matrix defined such that a vector vEcef
  //               defined relative to rEcef and expressed in ECEF coordinates
  //               can be expressed in the local ENU coordinate frame with
  //               origin at rEcef by vEnu = Renu_ecef*vEcef.
  //
  Eigen::Matrix3d ecef2enu(const Eigen::Vector3d& rEcef) {

    //----- Define WGS-84 Earth parameters
    const f64 aa = AA_WGS84;
    const f64 bb = BB_WGS84;
    const f64 ee = e_WGS84;
    const f64 ep = sqrt((aa*aa - bb*bb)/(bb*bb));

    //----- Convert to (phi,lambda,h) geodetic coordinates
    f64 x = rEcef(0), y = rEcef(1), z = rEcef(2);
    f64 lambda = atan2(y, x);
    f64 p = sqrt(x*x + y*y);
    f64 theta = atan2(z*aa, p*bb);
    f64 phi = atan2(z + ep*ep*bb*pow(sin(theta),3),
                    p - ee*ee*aa*pow(cos(theta),3));

    //----- Form the rotation matrix
    Eigen::Matrix3d Renu_ecef;
    Renu_ecef(0,0) = -sin(lambda);
    Renu_ecef(0,1) = cos(lambda);
    Renu_ecef(0,2) = 0;
    Renu_ecef(1,0) = -sin(phi)*cos(lambda);
    Renu_ecef(1,1) = -sin(phi)*sin(lambda);
    Renu_ecef(1,2) = cos(phi);
    Renu_ecef(2,0) = cos(phi)*cos(lambda);
    Renu_ecef(2,1) = cos(phi)*sin(lambda);
    Renu_ecef(2,2) = sin(phi);

    return Renu_ecef;
  }

  // Converts the frame-rotation-type direction-cosine matrix dc to a
  // quaternion.
  //
  // Notes: Given that both quat and -quat produce the same direction cosine
  // matrix, the convention adopted here is that the scalar component of the
  // quaternion; i.e., quat(3), will be nonnegative.  
  //
  // INPUTS
  //
  // dc ---------- 3-by-3 direction cosine matrix.
  //
  //
  // OUTPUTS
  //
  // quat -------- 4-by-1 quaternion corresponding to dc.
  //
  Vector4d dc2quat(const Eigen::Matrix3d& dc) {

    Vector4d quat = Eigen::MatrixXd::Zero(4,1);
    f64 tr = dc.trace();
    f64 s;
    if(tr > 0) {
      // trace is positive
      s = sqrt(tr+1.0);
      quat(3) = s/2.0;
      s = 0.5/s;
      quat(0) = (dc(1,2)-dc(2,1))*s;
      quat(1) = (dc(2,0)-dc(0,2))*s;
      quat(2) = (dc(0,1)-dc(1,0))*s;
    }
    else {
      // trace is negative or zero
      s32 nxt[3] = {1,2,0}, i,j,k;
      i = 0;
      if (dc(1,1)>dc(0,0))
        i=1;
      if (dc(2,2)>dc(i,i))
        i=2;
      j = nxt[i];
      k = nxt[j];
      s = sqrt((dc(i,i)-(dc(j,j)+dc(k,k)))+1.0);
      quat(i) = s*0.5;
      if (s != 0.0)
        s=0.5/s;
      quat(3) = (dc(j,k)-dc(k,j))*s;
      quat(j) = (dc(i,j)+dc(j,i))*s;
      quat(k) = (dc(i,k)+dc(k,i))*s;

      if (quat(3)<0.0)
        quat = quat*(-1);
    }

    return quat;
  }

  // Converts the quaternion quat into frame-rotation-type direction-cosine
  // matrix dc.
  //
  //
  // INPUTS
  //
  // quat -------- 4-by-1 quaternion vector quat = (q1;q2;q3;q4) = (qbar;q4)
  //               in which the vector component qbar = (q1;q2;q3) is the
  //               Euler axis of rotation scaled by sin(theta/2) and the 4th
  //               component is cos(theta/2), where theta is the rotation
  //               angle.
  //
  //
  // OUTPUTS
  //
  // dc ---------- 3-by-3 direction cosine matrix corresponding to quat.
  //
  Eigen::Matrix3d quat2dc(const Vector4d& quat) {

    f64 q0, q1, q2, q3;
    q0 = quat(0);
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);

    Eigen::Matrix3d dc = Eigen::MatrixXd::Zero(3,3);
    dc(0,0) = q0*q0 - q1*q1 - q2*q2 + q3*q3;
    dc(0,1) = 2.0*(q0*q1 + q2*q3);
    dc(0,2) = 2.0*(q0*q2 - q1*q3);
    dc(1,0) = 2.0*(q0*q1 - q2*q3);
    dc(1,1) = -q0*q0 + q1*q1 - q2*q2 + q3*q3;
    dc(1,2) = 2.0*(q1*q2 + q0*q3);
    dc(2,0) = 2.0*(q0*q2 + q1*q3);
    dc(2,1) = 2.0*(q1*q2 - q0*q3);
    dc(2,2) = -q0*q0 - q1*q1 + q2*q2 + q3*q3;

    return dc;
  }

  // q = qmult(q1,q2) performs quaternion multiplication on 4-by-1 quaternions
  // q1 and q2.  The multiplication is in the same order as the equivalent
  // direction-cosine multiplication A(q1)*A(q2).
  //
  // INPUTS
  //
  // q1 -------- 4-by-1 left quaternion operand.
  //
  // q2 -------- 4-by-1 right quaternion operand.
  //
  // OUTPUTS
  //
  // q --------- (return value) 4-by-1 quaternion product of q1 and q2.
  //
  Vector4d qmult(const Vector4d& q1, const Vector4d& q2) {

    Eigen::Matrix4d Q;

    Q(0,0) = q1(3);
    Q(0,1) = q1(2);
    Q(0,2) = -q1(1);
    Q(0,3) = q1(0);

    Q(1,0) = -q1(2);
    Q(1,1) = q1(3);
    Q(1,2) = q1(0);
    Q(1,3) = q1(1);

    Q(2,0) = q1(1);
    Q(2,1) = -q1(0);
    Q(2,2) = q1(3);
    Q(2,3) = q1(2);

    Q(3,0) = -q1(0);
    Q(3,1) = -q1(1);
    Q(3,2) = -q1(2);
    Q(3,3) = q1(3);

    return Q*q2;
  }

  // Converts the input differential quaternion in reduced dimensional form to
  // a full quaternion.
  //
  // INPUTS:
  // de -- 3-by-1 differential quaternion
  //
  // OUTPUTS:
  // dq -- 4-by-1 quaternion representation of the differential quaternion
  //
  Vector4d diffQuat2FullQuat(const Eigen::Vector3d& de) {
    Vector4d dq;
    f64 deNorm = de.norm();
    dq.head(3) = de;
    if(deNorm <= 1)
      dq(3) = sqrt(1 - deNorm*deNorm);
    else {
      dq(3) = 1;
      dq /= sqrt(1 + deNorm*deNorm);
    }

    return dq;
  }

  // Converts Euler angles phi = e(0), theta = e(1), and psi = e(2) (in
  // radians) into a direction cosine matrix for a 3-1-2 rotation.
  //
  // Let the world (W) and body (B) reference frames be initially aligned.  In
  // a 3-1-2 order, rotate B away from W by angles psi (yaw, about the body Z
  // axis), phi (roll, about the body X axis), and theta (pitch, about the
  // body Y axis).  RBW can then be used to cast a vector expressed in W
  // coordinates as a vector in B coordinates: vB = RBW * vW.
  //
  Eigen::Matrix3d euler2dc(const Eigen::Vector3d& e) {
    const f64 cPhi = cos(e(0)), sPhi = sin(e(0)), cThe = cos(e(1)),
      sThe = sin(e(1)), cPsi = cos(e(2)), sPsi = sin(e(2));
    Eigen::Matrix3d RBW;
    RBW <<
      cPsi*cThe - sPhi*sPsi*sThe, cThe*sPsi + cPsi*sPhi*sThe, -cPhi*sThe,
      -cPhi*sPsi,                                  cPhi*cPsi,       sPhi,
      cPsi*sThe + cThe*sPhi*sPsi, sPsi*sThe - cPsi*cThe*sPhi,  cPhi*cThe;
    return RBW;
  }

  // Converts a direction cosine matrix RBW to Euler angles phi = e(0), theta =
  // e(1), and psi = e(2) (in radians) for a 3-1-2 rotation.
  //
  // Let the world (W) and body (B) reference frames be initially aligned.  In a
  // 3-1-2 order, rotate B away from W by angles psi (yaw, about the body Z
  // axis), phi (roll, about the body X axis), and theta (pitch, about the body Y
  // axis).  RBW can then be used to cast a vector expressed in W coordinates as
  // a vector in B coordinates: vB = RBW * vW
  //
  // Guarantees (ignoring round-off):
  // 1. If |e(0)| < pi/2 - 1.4e-7, then dc2euler(euler2dc(e)) == e, modulo 2pi.
  // 2. If det(M) = +1, M is orthonormal, and no singularities are
  //    encountered, then euler2dc(dc2euler(M)) == M.
  Eigen::Vector3d dc2euler(const Eigen::Matrix3d& RBW) {
    constexpr f64 epsilon = 1e-14;
    // Note that if both arguments of the atan2 below are zero, as occurs when
    // cos(phi) = 0, then the result is undefined. Thus, for phi = pi/2 + n*pi, for
    // n an integer, the 3-1-2 Euler representation is singular.  The intuition
    // here is that when the roll angle is pi/2, the first and third rotations
    // (about the z and y axes, respectively), have exactly the same effect, and so
    // can't be distinguished, leading to a non-unique representation for psi and
    // theta; only psi + theta is constrained, but not psi and theta individually.
    // This function detects this singularity and returns a zero vector.
    if(std::abs(RBW(1,2)) > 1 - epsilon)
      return Eigen::Vector3d::Zero();
    const f64 phi = asin(RBW(1,2));
    const f64 theta = atan2(-RBW(0,2),RBW(2,2));
    const f64 psi = atan2(-RBW(1,0),RBW(1,1));
    Eigen::Vector3d e;
    e << phi, theta, psi;
    return e;
  }

  Eigen::Vector4d euler2quat(const Eigen::Vector3d& e) {
    return dc2quat(euler2dc(e));
  }
  
  Eigen::Vector3d quat2euler(const Vector4d& q) {
    return dc2euler(quat2dc(q));
  }

  // Returns the yaw angle corresponding to the input quaternion
  f64 heading(const Vector4d& q) {
    Eigen::Vector3d e = quat2euler(q);
    return e(2);    
  }

  // Returns the full rotation angle (in radians) of the direction cosine
  // matrix dc
  f64 rotationAngle(const Eigen::Matrix3d& dc) {
    // Let a be the full rotation angle.  Then the trace of dc is equal to 1 +
    // 2*cos(a).
    return acos((dc.trace() - 1)/2);
  }

  // Outputs the gravitational acceleration vector at the input ECEF position.
  //
  // Notes: This version is a simple second-order WGS 84 model (includes J2 but
  // nothing higher).
  //
  // INPUTS
  //
  // rECEF ------- 3-by-1 position vector in ECEF meters.
  //
  //
  // OUTPUTS
  //
  // gVec -------- 3-by-1 gravitational acceleration vector expressed in ECEF
  //               coordinates in units of m/s^2.
  //
  Eigen::Vector3d getGravitationalAccelerationVector(const Eigen::Vector3d& rEcef) {

    //----- WGS84 constants
    const f64 aa = AA_WGS84;       // m, ellipsoid semi-major axis
    const f64 GM = 3986004.418e8;  // m^3/s^2, Earth's gravitational constant
    const f64 J2 = 1.08263e-3;     // J2 coefficient

    //----- Intermediate constants
    const f64 R2 = rEcef.dot(rEcef);
    const f64 R  = sqrt(R2);
    const f64 aaoR2 = aa*aa/(R*R);
    const f64 zoR2 = 5.0*rEcef(2)*rEcef(2)/(R*R);

    //----- Calculation of gVec
    Eigen::Vector3d gVec;
    gVec(0) = -(GM/R2)*(1.0 + 1.5*J2*aaoR2*(1.0-zoR2))*(rEcef(0)/R);
    gVec(1) = -(GM/R2)*(1.0 + 1.5*J2*aaoR2*(1.0-zoR2))*(rEcef(1)/R);
    gVec(2) = -(GM/R2)*(1.0 + 1.5*J2*aaoR2*(3.0-zoR2))*(rEcef(2)/R);

    return gVec;
  }

  // Returns the direction cosine matrix corresponding to a rotation of theta
  // radians about the x-axis.
  //
  Eigen::Matrix3d Rx(const f64 theta) {
    Eigen::Matrix3d R = Eigen::MatrixXd::Identity(3,3);
    f64 c = cos(theta);
    f64 s = sin(theta);
    R(1,1) =  c; R(1,2) = s;
    R(2,1) = -s; R(2,2) = c;
    return R;
  }

  // Returns the direction cosine matrix corresponding to a rotation of theta
  // radians about the y-axis.
  //
  Eigen::Matrix3d Ry(const f64 theta) {
    Eigen::Matrix3d R = Eigen::MatrixXd::Identity(3,3);
    f64 c = cos(theta);
    f64 s = sin(theta);
    R(0,0) =  c; R(0,2) = -s;
    R(2,0) =  s; R(2,2) =  c;
    return R;
  }

  // Returns the direction cosine matrix corresponding to a rotation of theta
  // radians about the z-axis.
  //
  Eigen::Matrix3d Rz(const f64 theta) {
    Eigen::Matrix3d R = Eigen::MatrixXd::Identity(3,3);
    f64 c = cos(theta);
    f64 s = sin(theta);
    R(0,0) =  c; R(0,1) = s;
    R(1,0) = -s; R(1,1) = c;
    return R;
  }

  // Returns the cross-product-equivalent matrix for the input vec.
  //
  Eigen::Matrix3d cpe(const Eigen::Vector3d& vec) {
    Eigen::Matrix3d crossMat;
    crossMat(0, 1) = -vec(2);
    crossMat(0, 2) = vec(1);
    crossMat(1, 0) = vec(2);
    crossMat(1, 2) = -vec(0);
    crossMat(2, 0) = -vec(1);
    crossMat(2, 1) = vec(0);
    return crossMat;
  }

  // Solves Wahba's problem via SVD.  In other words, this function finds the
  // rotation matrix RBW that minimizes the cost Jw:
  //
  //                     N
  //    Jw(RBW) = (1/2) sum ai*||viB - RBW*viW||^2
  //                    i=1
  //
  // where viB and viW are the ith vector representations in body and world
  // frame, respectively, and where ai is the ith weight.  Note that a
  // solution only exists for N > 1 (at least two pairs of vectors are
  // required).
  //
  // INPUTS
  //
  // aVec ------- Nx1 vector of least-squares weights.  aVec(i-1) is the weight
  //              corresponding to the ith pair of vectors.  
  //
  // vWMat ------ 3xN matrix of 3x1 vectors expressed in the W frame.
  //              vWMat.col(i-1) is the ith 3x1 vector.
  //
  // vBMat ------ 3xN matrix of 3x1 vectors expressed in the B
  //              frame. vBMat.col(i-1) is the ith 3x1 vector, which
  //              corresponds to vWMat.col(i-1).
  //
  // OUTPUTS
  // 
  // RBW -------- 3x3 direction cosine matrix indicating the attitude of the
  //              B frame relative to the W frame, such that vB = RBW*vW.
  //
  Eigen::Matrix3d wahbaSolver(const Eigen::VectorXd& aVec,
                              const Eigen::MatrixXd& vWMat,
                              const Eigen::MatrixXd& vBMat) {

    using namespace Eigen;
    
    Matrix3d B = MatrixXd::Zero(3,3);
    const u32 N = aVec.size();
    for(u32 ii=0; ii<N; ii++) 
      B += aVec(ii)*vBMat.col(ii)*vWMat.col(ii).transpose();
    JacobiSVD<Matrix3d> svd(B,ComputeFullU | ComputeFullV);
    Matrix3d U = svd.matrixU();
    Matrix3d V = svd.matrixV();
    Vector3d m;
    m << 1, 1, U.determinant()*V.determinant();
    return U*m.asDiagonal()*V.transpose();
  }

} // namespace
