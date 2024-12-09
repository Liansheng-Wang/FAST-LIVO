#ifndef Manifold_MATH_H_H_
#define Manifold_MATH_H_H_


#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>


#include "alias.h"


namespace LI2Sup {

template <typename Derived>
inline Eigen::Matrix<typename Derived::Scalar, 3, 3> SKEW(const Eigen::MatrixBase<Derived> &v) {
  Eigen::Matrix<typename Derived::Scalar, 3, 3>  m;
  m << typename Derived::Scalar(0), -v[2], v[1], 
      v[2], typename Derived::Scalar(0), -v[0], 
      -v[1], v[0], typename Derived::Scalar(0);
  return m;
}


template <typename T>
inline Eigen::Matrix<T, 3, 3> Exp(const Eigen::Matrix<T, 3, 1> &&ang) {
  T ang_norm = ang.norm();
  if (ang_norm > 1e-10) {
    Eigen::Matrix<T, 3, 1> r_axis = ang / ang_norm;
    Eigen::Matrix<T, 3, 3> K;
    K = SKEW(r_axis);
    /// Roderigous Tranformation
    return Eye3 + std::sin(ang_norm) * K + (1.0 - std::cos(ang_norm)) * K * K;
  } else {
    return Eye3;
  }
}

// left Jacobian
template <typename T>
inline Eigen::Matrix<T, 3, 3> Gamma_1(const Eigen::Matrix<T, 3, 1>& ang){
  T ang_norm = ang.norm();
  if (ang_norm > 1e-10) {
    Eigen::Matrix<T, 3, 1> r_axis = ang / ang_norm;
    Eigen::Matrix<T, 3, 3> K;
    K = SKEW(r_axis);
    return Eye3 + ((1.0 - std::cos(ang_norm)) * K + (ang_norm - std::sin(ang_norm))* K * K) / ang_norm;
  } else {
    return Eye3;
  }
}


// Second-order differential variable
template <typename T>
inline Eigen::Matrix<T, 3, 3> Gamma_2(const Eigen::Matrix<T, 3, 1>& ang){
  T ang_norm = ang.norm();
  if (ang_norm > 1e-5) {
    Eigen::Matrix<T, 3, 1> r_axis = ang / ang_norm;
    Eigen::Matrix<T, 3, 3> K;
    K = SKEW(r_axis);
    /// Roderigous Tranformation
    return 0.5 * Eye3 + ( 2.0 * (ang_norm - std::sin(ang_norm)) * K 
                          + ((ang_norm * ang_norm + 2 * std::cos(ang_norm) - 2) * K ) ) 
                        / ( 2 * ang_norm * ang_norm);
  } else {
    return 0.5 * Eye3;
  }
}


template <typename T, typename Ts>
inline Eigen::Matrix<T, 3, 3> Exp(const Eigen::Matrix<T, 3, 1> &ang_vel, const Ts &dt) {
  T ang_vel_norm = ang_vel.norm();
  Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();

  if (ang_vel_norm > 0.0000001) {
    Eigen::Matrix<T, 3, 1> r_axis = ang_vel / ang_vel_norm;
    Eigen::Matrix<T, 3, 3> K;

    K = SKEW(r_axis);

    T r_ang = ang_vel_norm * dt;

    return Eye3 + std::sin(r_ang) * K + (1.0 - std::cos(r_ang)) * K * K;
  } else {
    return Eye3;
  }
}

template <typename T>
inline Eigen::Matrix<T, 3, 3> Exp(const T &v1, const T &v2, const T &v3) {
  T &&norm = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
  Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();
  if (norm > 1e-10) {
    Eigen::Matrix<T, 3, 1> r_ang = {v1 / norm, v2 / norm, v3 / norm};
    Eigen::Matrix<T, 3, 3> K;
    K = SKEW(r_ang);
    return Eye3 + std::sin(norm) * K + (1.0 - std::cos(norm)) * K * K;
  } else {
    return Eye3;
  }
}

/* Logrithm of a Rotation Matrix */
template <typename T>
inline Eigen::Matrix<T, 3, 1> Log(const Eigen::Matrix<T, 3, 3> &R) {
  T theta = (R.trace() > 3.0 - 1e-6) ? 0.0 : std::acos(0.5 * (R.trace() - 1));
  Eigen::Matrix<T, 3, 1> K(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
  return (std::abs(theta) < 0.001) ? (0.5 * K) : (0.5 * theta / std::sin(theta) * K);
}

template <typename T>
inline Eigen::Matrix<T, 3, 1> RotMtoEuler(const Eigen::Matrix<T, 3, 3> &rot) {
  T sy = sqrt(rot(0, 0) * rot(0, 0) + rot(1, 0) * rot(1, 0));
  bool singular = sy < 1e-6;
  T x, y, z;
  if (!singular) {
    x = atan2(rot(2, 1), rot(2, 2));
    y = atan2(-rot(2, 0), sy);
    z = atan2(rot(1, 0), rot(0, 0));
  } else {
    x = atan2(-rot(1, 2), rot(1, 1));
    y = atan2(-rot(2, 0), sy);
    z = 0;
  }
  Eigen::Matrix<T, 3, 1> ang(x, y, z);
  return ang;
}


inline M4 hat6(const V6 &xi){
  M4 xi_hat;
  xi_hat  <<  0.0,   -xi(2),  xi(1), xi(3),
              xi(2),  0.0,   -xi(0), xi(4),
             -xi(1),  xi(0),  0.0,   xi(5),
              0.0,    0.0,    0.0,   0.0;
  return xi_hat;
}

inline V6 vee6(const M4 &xi_hat){
  V6 xi;
  xi << -xi_hat(1,2), xi_hat(0,2), -xi_hat(0,1),
          xi_hat(0,3), xi_hat(1,3), xi_hat(2,3);
  return xi;
}





/**
 * SO3 and SE3 ref to https://github.com/hku-mars/BALM/blob/master/src/compare_test/SE3/SE3.hpp
 */
class SO3
{
public:
  SO3(const M3 &R = Eye3):R_(R){}

  SO3(const V3 &w):R_(Eye3){
    this->exp(hat(w));
  }

  SO3(const SO3 &R):R_(R.R()){}

  template<typename OtherDerived>
  SO3(const Eigen::MatrixBase<OtherDerived>& rhs):
  R_(rhs)
  {}


  static M3 hat(const V3 &v){
    M3 m;
    m << 0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0;
    return m;
  }

  static V3 vee(const M3 &w_hat){
    V3 w;
    w << -w_hat(1,2), w_hat(0,2), -w_hat(0,1);
    return w;
  }

  static SO3 Exp(const V3 &ang) {
    scalar &&norm = ang.norm();
    M3 I = Eye3;
    if (norm > 1e-10) {
      V3 r_ang = ang / norm;
      M3 K;
      K = hat(r_ang);
      return SO3(I + std::sin(norm) * K + (1.0 - std::cos(norm)) * K * K);
    } else {
      return SO3(I);
    }
  }

  static SO3 Exp(const V3 &ang_vel, const scalar &dt) {
    scalar ang_vel_norm = ang_vel.norm();
    M3 I = Eye3;
    if (ang_vel_norm > 1e-12) {
      V3 r_axis = ang_vel / ang_vel_norm;
      M3 K;
      K = hat(r_axis);
      scalar r_ang = ang_vel_norm * dt;
      return SO3(I + std::sin(r_ang) * K + (1.0 - std::cos(r_ang)) * K * K);
    } else {
      return SO3(I);
    }
  }

  SO3& operator=(const SO3& rhs){
    if (this == &rhs)
        return *this;
    R_ = rhs.R();
    return *this;
  }

  SO3 operator*(const SO3& rhs) const{
    M3 res = R_ * rhs.R();
    return SO3(res);
  }

  SO3 operator*(scalar s) const {
    M3 res = R_ * s;
    return SO3(res);
  }

  V3 operator*(const V3& rhs) const{
    return R_ * rhs;
  }

  template<typename OtherDerived>
  V3 operator*(const Eigen::MatrixBase<OtherDerived>& rhs) const{
    return R_ * rhs;
  }

  void update(const V3 &dw){
    SO3 dR(dw);
    R_ = dR.R() * R_;
  }

  void updateRhs(const V3 &dw){
    SO3 dR(dw);
    R_ = R_ * dR.R();
  }

  void exp(const M3 &w_hat){
    V3 w = vee(w_hat);
    double o = w.norm();
    if ( o < 1e-12){
        R_ << Eye3 + w_hat;
        return;
    }
    double c1 = std::sin(o)/o;
    double c2 = (1 - std::cos(o))/o/o;
    R_ << Eye3 + c1 * w_hat + c2 * w_hat *w_hat;
  }

  V4 coeffs() const{
    return Quat(R_).coeffs();
  }

  M3 log(double *ro = nullptr) const{
    M3 res;
    double tr = (R_.trace()-1)*0.5;
    double o;
    if (tr  < 1.0 - 1e-9 && tr > -1.0 + 1e-9 )
    {
      o = std::fabs(std::acos(tr));
      res << 0.5 * o / std::sin(o) * ( R_ - R_.transpose());
    }
    else if (tr >= 1.0 - 1e-9 )
    {
      o = 0.0;
      res << M3::Zero();
    }
    else
    {
      o = M_PI;
      V3 w;
      if( R_(0,0) > R_(1,1) && R_(0,0) > R_(2,2) )
      {
        w << R_(0,0) + 1.0,
            0.5 * ( R_(0,1) + R_(1,0)),
            0.5 * ( R_(0,2) + R_(2,0));
      }
      else if( R_(1,1) > R_(0,0) && R_(1,1) > R_(2,2) )
      {
        w << 0.5 * ( R_(1,0) + R_(0,1)),
              R_(1,1) + 1.0,
              0.5 * ( R_(1,2) + R_(2,1));
      }
      else
      {
        w << 0.5 * ( R_(2,0) + R_(0,2)),
              0.5 * ( R_(2,1) + R_(1,2)),
              R_(2,2) + 1.0;
      }
      double length = w.norm();
      if (length > 0.0)
      {
        w *= M_PI / length;
      }
      else
      {
        w << 0.0, 0.0, 0.0;
      }
      res = hat(w);
    }
    if (ro != nullptr) *ro = o;
    return res;
  }

  V3 log_vee() const{
    M3 w_hat = this->log();
    return vee(w_hat);
  }

  SO3 inverse(void) const{
    return SO3(R_.transpose());
  }

  M3 adjoint() const{
    return R_.transpose();
  }

  M3 R() const{
    return R_;
  }

  Quat quaternion() const{
    Quat res = Quat(R_);
    return res.normalized();
  }

  M3 matrix() const{
    return R_;
  }

  M3& ref2R(){
    return R_;
  }

  double distance(const SO3 &rhs) const{
    return (*this * rhs.inverse()).log_vee().norm();
  }

  void print(void) const{
    std::cout << R_ << std::endl;
  }

  void print_lie(void) const{
    V3 w =  this->log_vee();
    std::cout << w << std::endl;
  }

public:
  M3 R_;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};



class SE3
{
public:
  SE3(const M4 &T = M4::Identity()): T_(T) {
    R_ = T_.topLeftCorner<3,3>();
    t_ = T_.topRightCorner<3,1>();
  }

  SE3(const V6 &xi){
    this->exp(hat6(xi));
    R_ = T_.topLeftCorner<3,3>();
    t_ = T_.topRightCorner<3,1>();
  }
  
  SE3(const SE3 &T): T_(T.T()),R_(T.R()),t_(T.t()){}

  SE3(const SO3 &R, const V3 &t):
  T_(M4::Identity())
  {
    T_.block<3,3>(0,0) = R.R();
    T_.block<3,1>(0,3) = t;
    R_ = R.R();
    t_ = t;
  }

  SE3(const M3 &R, const V3 &t):
  T_(M4::Identity())
  {
    T_.block<3,3>(0,0) = R;
    T_.block<3,1>(0,3) = t;
    R_ = R;
    t_ = t;
  }

  SE3(const Quat &q, const V3 &t):
  T_(M4::Identity())
  {
    T_.block<3,3>(0,0) = q.toRotationMatrix();
    T_.block<3,1>(0,3) = t;
    R_ = T_.topLeftCorner<3,3>();
    t_ = t;
  }

  template<typename OtherDerived>
  SE3(const Eigen::MatrixBase<OtherDerived>& rhs)  :
  T_(rhs)
  { }

  SE3& operator=(const SE3& rhs){
    if (this == &rhs)
        return *this;
    T_ = rhs.T();
    R_ = rhs.R();
    t_ = rhs.t();
    return *this;
  }

  SE3 operator*(const SE3& rhs) const{
    M4 res = T_ * rhs.T();
    return SE3(res);
  }

  V3  operator*(const V3& point) const{
    return R_*point + t_;
  }

  V4  operator*(const V4& point) const{
    return T_ * point;
  }

  void update(const V6 &dxi){
    SE3 dT(dxi);
    T_ = dT.T() * T_;
  }

  void updateRhs(const V6 &dxi){
    SE3 dT(dxi);
    T_ = T_ * dT.T();
    R_ = T_.topLeftCorner<3,3>();
    t_ = T_.topRightCorner<3,1>();
  }

  void exp(const M4 &xi_hat){
    V6 xi = vee6(xi_hat);
    V3 w = xi.head<3>();
    V3 v = xi.tail<3>();
    SO3 rotation(w);
    M3 w_hat = xi_hat.topLeftCorner<3,3>();

    M3 V = M3::Identity();
    double o = w.norm();
    if ( o > 1e-12){
      double c2 = (1 - std::cos(o))/o/o;
      double c3 = (o - std::sin(o))/o/o/o;
      V += c2*w_hat + c3*w_hat*w_hat;
    }
    V3 t = V * v;

    T_  << rotation.R(), t,
           0,0,0,1;
    R_ = rotation.R();
    t_ = t;
  }

  M4 log(void) const{
    SO3 rotation(this->R());
    double o;
    M3 w_hat = rotation.log(&o);
    M3 Vinv = M3::Identity();
    if (o > 1e-12)
    {
      double c1 = std::sin(o);
      double c2 = (1 - std::cos(o))/o;
      double k1 = 1/o/o*(1 - 0.5*c1/c2);
      Vinv += -0.5*w_hat + k1* w_hat*w_hat;
    }
    V3 v = Vinv * T_.topRightCorner<3,1>();

    M4 xi_hat = M4::Zero();
    xi_hat << w_hat, v,
              0,0,0,0;
    return xi_hat;
  }

  V6 log_vee() const{
    M4 xi_hat = this->log();
    return vee6(xi_hat);
  }

  V3 transform(const V3 & p) const{
    return R_*p + t_;
  }

  SE3 inverse(void) const{
    M4 inv;
    M3 R = R_;
    R.transposeInPlace();
    inv << R, -R * this->t(),
           0,0,0,1;
    return SE3(inv);
  }

  M6 adjoint() const{
    M6 res(M6::Zero());
    M3 tx = SO3::hat( this->t() );
    res.topLeftCorner<3,3>() << R();
    res.bottomRightCorner<3,3>() << R();
    res.bottomLeftCorner<3,3>() << tx*R();
    return res;
  }

  M4 T() const{
    return T_;
  }

  M4 matrix() const{
    return T_;
  }

  M4& ref2T(){
    return T_;
  }

  SO3 so3() const{
    return SO3(R_);
  }

  M3 R() const{
    return R_;
  }

  V3 translation() const{
    return t_;
  }

  V3 t() const{
    return t_;
  }

  Quat quaternion() const{
    Quat res = Quat(R());
    return res.normalized();
  }

  double distance(const SE3 &rhs) const{
    return (*this * rhs.inverse()).log_vee().norm();
  }

  void print(void) const{
    std::cout << T_ << std::endl;
  }

  void print_lie(void) const{
    std::cout << this->log_vee() << std::endl;
  }

public:
  M4 T_;
  M3 R_;
  V3 t_;
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


/// only can be used for gravity vector.
class S2
{
public:
  V3 dir_;
  scalar length_;

public:
  S2(){
   dir_ = {0,0,1};
   length_ = 1; 
  }
  
  S2(V3 input): dir_(input) {
    length_ = input.norm();
  }
  
  S2(scalar d1, scalar d2, scalar d3){
    V3 input = V3(d1,d2,d3);
    dir_ = input;
    length_ = input.norm();
  }

  // static M3_2 S2_Bx(const V3& __g, const scalar length) {
  //   M3_2 res;
  //   if(__g[0] + length > 1e-10)
  //   {
  //     res << -__g[1], -__g[2],
  //         length - __g[1]*__g[1]/(length+__g[0]), -__g[2]*__g[1]/(length+__g[0]),
  //         -__g[2]*__g[1]/(length+__g[0]), length-__g[2]*__g[2]/(length+__g[0]);
  //     res /= length;
  //   }
  //   else
  //   {
  //     res = Eigen::Matrix<scalar, 3, 2>::Zero();
  //     res(1, 1) = -1;
  //     res(2, 0) = 1;
  //   }
  //   return res;
  // }

  static M3_2 S2_Bx(const V3& __g) {
    
    M3_2 B_x;
    V3 g = __g;
    g.normalize();

    B_x(0, 0) = -g[1];
    B_x(0, 1) = -g[2];
    B_x(1, 0) = 1 - g[1]*g[1]/(1+g[0]);
    B_x(1, 1) = - g[2]*g[1]/(1+g[0]);
    B_x(2, 0) = B_x(1, 1);
    B_x(2, 1) = 1-g[2]*g[2]/(1+g[0]);

    return B_x;
  }
  

  // static M3_2 S2_Bx(const V3& __g) {
  //   M3_2 res;
  //   V3 basis = V3::UnitX();
  //   res.col(0) = __g.cross(basis).normalized();
  //   res.col(1) = __g.cross(res.col(0)).normalized();
  //   return res;
  // }

  void update (const V2& delta){
    scalar one_z = length_ + dir_[2];
    scalar xx = dir_[0] * dir_[0];
    scalar xy = dir_[0] * dir_[1];
    scalar yy = dir_[1] * dir_[1];

    Eigen::Matrix<scalar, 3, 2> b;

    b << length_ - xx / one_z, - xy / one_z,
         - xy / one_z, length_ - yy / one_z,
         - dir_[0], - dir_[1];

    dir_ = SO3::Exp(b * delta) * dir_;
  }

};

}

#endif

