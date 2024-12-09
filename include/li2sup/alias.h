/**
 * @file eigen_types.h
 * @author WangLiansheng (lswang@mail.ecust.edu.cn)
 * @date 2023-10-06
 * @copyright Copyright (c) 2023
 */

#ifndef ALIAS_TYPES_H_H_
#define ALIAS_TYPES_H_H_

#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>


namespace LI2Sup {

using scalar = double;
// using scalar = float;

/* alias for eigen */
using V2 = Eigen::Matrix<scalar, 2, 1>;
using V3 = Eigen::Matrix<scalar, 3, 1>;
using V4 = Eigen::Matrix<scalar, 4, 1>;
using V5 = Eigen::Matrix<scalar, 5, 1>;
using V6 = Eigen::Matrix<scalar, 6, 1>;
using V7 = Eigen::Matrix<scalar, 7, 1>;
using V8 = Eigen::Matrix<scalar, 8, 1>;
using V9 = Eigen::Matrix<scalar, 9, 1>;
using V12 = Eigen::Matrix<scalar, 12, 1>;
using V15 = Eigen::Matrix<scalar, 15, 1>;
using V18 = Eigen::Matrix<scalar, 18, 1>;
using V24 = Eigen::Matrix<scalar, 24, 1>;
using V30 = Eigen::Matrix<scalar, 30, 1>;
using V36 = Eigen::Matrix<scalar, 36, 1>;
using V40 = Eigen::Matrix<scalar, 40, 1>;
using VX  = Eigen::Matrix<scalar, -1, 1>;

using VV3 = std::vector<V3, Eigen::aligned_allocator<V3>>;
using VV4 = std::vector<V4, Eigen::aligned_allocator<V4>>;
using VV5 = std::vector<V5, Eigen::aligned_allocator<V5>>;

using M1 = Eigen::Matrix<scalar, 1, 1>;
using M2 = Eigen::Matrix<scalar, 2, 2>;
using M3 = Eigen::Matrix<scalar, 3, 3>;
using M4 = Eigen::Matrix<scalar, 4, 4>;
using M5 = Eigen::Matrix<scalar, 5, 5>;
using M6 = Eigen::Matrix<scalar, 6, 6>;
using M7 = Eigen::Matrix<scalar, 7, 7>;
using M8 = Eigen::Matrix<scalar, 8, 8>;
using M9 = Eigen::Matrix<scalar, 9, 9>;
using M12 = Eigen::Matrix<scalar, 12, 12>;
using M15 = Eigen::Matrix<scalar, 15, 15>;
using M18 = Eigen::Matrix<scalar, 18, 18>;
using M24 = Eigen::Matrix<scalar, 24, 24>;

using M3_2 = Eigen::Matrix<scalar, 3, 2>;
using M2_3 = Eigen::Matrix<scalar, 2, 3>;

using MX3 = Eigen::Matrix<scalar, Eigen::Dynamic, 3>;

using Quat = Eigen::Quaternion<scalar>;

const M3 Eye3 = M3::Identity();
const V3 ZeroV3(0, 0, 0);
const V3 EyeV3(1,1,1);

using V2i = Eigen::Vector2i;
using V3i = Eigen::Vector3i;
using VXi = Eigen::VectorXi;

}
#endif
