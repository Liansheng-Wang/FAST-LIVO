cmake_minimum_required(VERSION 2.8.3)
project(fast_livo)

# 编译设置
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -O3")

# 定义根目录
add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

# 检测多核处理器
message(STATUS "Current CPU architecture: ${CMAKE_SYSTEM_PROCESSOR}")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
  include(ProcessorCount)
  ProcessorCount(N)
  if(N GREATER 5)
    add_definitions(-DMP_EN -DMP_PROC_NUM=4)
  elseif(N GREATER 3)
    math(EXPR PROC_NUM "${N} - 2")
    add_definitions(-DMP_EN -DMP_PROC_NUM="${PROC_NUM}")
  else()
    add_definitions(-DMP_PROC_NUM=1)
  endif()
else()
  add_definitions(-DMP_PROC_NUM=1)
endif()

# 查找依赖包
find_package(OpenMP QUIET)
find_package(PythonLibs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED COMPONENTS core imgproc highgui)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  sensor_msgs
  roscpp
  rospy
  std_msgs
  pcl_ros
  tf
  livox_ros_driver
  message_generation
  eigen_conversions
  vikit_common
  vikit_ros
  cv_bridge
  image_transport
)

# 输出调试信息
message(STATUS "Eigen include dirs: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "PCL libraries: ${PCL_LIBRARIES}")
message(STATUS "OpenCV libraries: ${OpenCV_LIBRARIES}")

# 包含路径
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  include
)

# 消息定义
add_message_files(FILES Pose6D.msg States.msg)
generate_messages(DEPENDENCIES geometry_msgs)

# 定义 catkin 包
catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs message_runtime cv_bridge image_transport vikit_common vikit_ros
  DEPENDS EIGEN3 PCL OpenCV Sophus
  INCLUDE_DIRS include
)

# 库和可执行文件
add_library(ikdtree include/ikd-Tree/ikd_Tree.cpp include/FOV_Checker/FOV_Checker.cpp)
add_library(vio src/lidar_selection.cpp src/frame.cpp src/point.cpp src/map.cpp)
add_executable(fastlivo_mapping src/laserMapping.cpp src/IMU_Processing.cpp src/preprocess.cpp)

# 链接库
target_link_libraries(ikdtree ${PCL_LIBRARIES})
target_link_libraries(fastlivo_mapping ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${PYTHON_LIBRARIES} vio ikdtree)
target_link_libraries(fastlivo_mapping PUBLIC OpenMP::OpenMP_CXX)