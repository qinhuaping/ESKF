#ifndef NDEBUG
#define NDEBUG
#endif

#include "ESKF.hpp"
#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace Eigen;

namespace eskf {
  
  std::ofstream debugFile("/home/elia/Documents/eskf.csv");
  static double curr_time_sec = 0.0; 
  
  template <typename T> inline T sq(T var) {
    return var * var;
  }
	
  template<typename Scalar>
  static inline constexpr const Scalar &constrain(const Scalar &val, const Scalar &min_val, const Scalar &max_val) {
    return (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
  }
  
  template<typename Type>
  Type wrap_pi(Type x) {
      while (x >= Type(M_PI)) {
          x -= Type(2.0 * M_PI);
      }

      while (x < Type(-M_PI)) {
          x += Type(2.0 * M_PI);
      }
      return x;
  }
  
  /**
   * Rotation quaternion from vector
   *
   * The axis of rotation is given by vector direction and
   * the angle is given by the norm.
   *
   * @param vec rotation vector
   * @return quaternion representing the rotation
   */
  ESKF::quat ESKF::from_axis_angle(ESKF::vec3 vec) {
    quat q;
    scalar_t theta = vec.norm();

    if (theta < scalar_t(1e-10)) {
      q.w() = scalar_t(1.0);
      q.x() = q.y() = q.z() = 0;
      return q;
    }

    //vec /= theta;
    vec3 tmp = vec / theta;
    return from_axis_angle(tmp, theta);
  }

  /**
   * Rotation quaternion from axis and angle
   * XXX DEPRECATED, use AxisAngle class
   *
   * @param axis axis of rotation
   * @param theta scalar describing angle of rotation
   * @return quaternion representing the rotation
   */
  ESKF::quat ESKF::from_axis_angle(const ESKF::vec3 &axis, scalar_t theta) {
    quat q;

    if (theta < scalar_t(1e-10)) {
      q.w() = scalar_t(1.0);
      q.x() = q.y() = q.z() = 0;
    }

    scalar_t magnitude = sin(theta / 2.0f);

    q.w() = cos(theta / 2.0f);
    q.x() = axis(0) * magnitude;
    q.y() = axis(1) * magnitude;
    q.z() = axis(2) * magnitude;
    
    return q;
  }
  
  ESKF::ESKF() {
    debugFile << "time_" << "," << "yaw_" << "," << "pitch_" << "," << "roll_" << std::endl;
    //debugFile << "time_" << "," << "dq1_" << "," << "dq2_" << "," << "dq3_" << "," << "dq4_" << std::endl;
    //debugFile << "time_" << "," << "q1_" << "," << "q2_" << "," << "q3_" << "," << "q4_" << std::endl;
    //debugFile << "time_" << "," << "corrected_delta_ang(0)_" << "," << "corrected_delta_ang(1)_" << "," << "corrected_delta_ang(2)_" << std::endl;
    // zeros state_
    state_.quat_nominal = quat(1, 0, 0, 0);
    state_.gyro_bias = vec3(0, 0, 0);
        
    //  zeros P_
    for (unsigned i = 0; i < k_num_states_; i++) {
      for (unsigned j = 0; j < k_num_states_; j++) {
        P_[i][j] = 0.0f;
      }
	  }
	
    dx_.setZero();
    
    rosb2px4b(0,0) = 1.0;
    rosb2px4b(0,1) = 0.0;
    rosb2px4b(0,2) = 0.0;
    
    rosb2px4b(1,0) =  0.0;
    rosb2px4b(1,1) = -1.0;
    rosb2px4b(1,2) =  0.0;
    
    rosb2px4b(2,0) =  0.0;
    rosb2px4b(2,1) =  0.0;
    rosb2px4b(2,2) = -1.0;
    
    q_ne = quat(0, 0.7071, 0.7071, 0); // rotation from ned to enu
    q_rb = quat(0, 1, 0, 0); // rotation ros body to px4 body
  }
  
  bool ESKF::initialize(scalar_t dt) {
    initialiseCovariance(dt);
    return true;
  }
	
  void ESKF::initialiseCovariance(scalar_t dt) {
     // define the initial angle uncertainty as variances for a rotation vector
	  vec3 rot_vec_var;
	  rot_vec_var(2) = rot_vec_var(1) = rot_vec_var(0) = sq(initial_tilt_err);

	  // update the quaternion state covariances
	  initialiseQuatCovariances(rot_vec_var);

    // gyro bias
	  P_[4][4] = sq(switch_on_gyro_bias * dt);
	  P_[5][5] = P_[4][4];
	  P_[6][6] = P_[4][4];  
  }
  
  void ESKF::predict(const ESKF::vec3 &w, ESKF::scalar_t dt) {
    static bool isFirstTime = false;
    if(!isFirstTime)
      isFirstTime = true;
    if(isFirstTime)
      firstPredict = true;
    
    // convert ROS body to PX4 body frame IMU data
    vec3 px4body_w = rosb2px4b * w;
    
    vec3 delta_angle = px4body_w * dt; // current delta angle  (rad)
    //debugFile << curr_time_sec << "," << delta_angle(0) << "," << delta_angle(1) << "," << delta_angle(2) << std::endl;
    // apply imu bias corrections
    vec3 corrected_delta_ang = delta_angle;//- vec3(3.8785e-05,3.8785e-05,3.8785e-05);//state_.gyro_bias;
    //debugFile << curr_time_sec << "," << corrected_delta_ang(0) << "," << corrected_delta_ang(1) << "," << corrected_delta_ang(2) << std::endl;
    // convert the delta angle to a delta quaternion
    quat dq;
    //corrected_delta_ang = vec3(0.1, 0.2, 0.3);
    dq = from_axis_angle(corrected_delta_ang);
    //debugFile << curr_time_sec << "," << dq.w() << "," << dq.x() << "," << dq.y() << "," << dq.z() << std::endl;
    // rotate the previous quaternion by the delta quaternion using a quaternion multiplication
    state_.quat_nominal = state_.quat_nominal * dq;
    // quaternions must be normalised whenever they are modified
    state_.quat_nominal.normalize();
    //debugFile << curr_time_sec << "," << state_.quat_nominal.w() << "," << state_.quat_nominal.x() << "," << state_.quat_nominal.y() << "," << state_.quat_nominal.z() << std::endl;
    constrainStates(dt);
    mat3 R_to_earth = quat_to_invrotmat(state_.quat_nominal);
    scalar_t yaw = atan2f(-R_to_earth(0, 1), R_to_earth(1, 1)); // first rotation (yaw)
    scalar_t pitch = atan2f(-R_to_earth(2, 0), R_to_earth(2, 2)); // third rotation (pitch)
    scalar_t roll = asinf(R_to_earth(2, 1)); // second rotation (roll)
    debugFile << curr_time_sec << "," << yaw << "," << pitch << "," << roll << std::endl;
    
    // error-state jacobian
    // assign intermediate state variables
    scalar_t q0 = state_.quat_nominal.w();
    scalar_t q1 = state_.quat_nominal.x();
    scalar_t q2 = state_.quat_nominal.y();
    scalar_t q3 = state_.quat_nominal.z();

    scalar_t dax = delta_angle(0);
    scalar_t day = delta_angle(1);
    scalar_t daz = delta_angle(2);

    scalar_t dax_b = state_.gyro_bias(0);
    scalar_t day_b = state_.gyro_bias(1);
    scalar_t daz_b = state_.gyro_bias(2);

    // compute noise variance for stationary processes
    scalar_t process_noise[k_num_states_] = {};

    // convert rate of change of rate gyro bias (rad/s**2) as specified by the parameter to an expected change in delta angle (rad) since the last update
    scalar_t d_ang_bias_sig = dt * dt * constrain(gyro_bias_p_noise, 0.0, 1.0);

    // Construct the process noise variance diagonal for those states with a stationary process model
    // These are kinematic states and their error growth is controlled separately by the IMU noise variances
    for (unsigned i = 0; i <= 3; i++) {
      process_noise[i] = 0.0;
    }
    
    // delta angle bias states
    process_noise[4] = process_noise[5] = process_noise[6] = sq(d_ang_bias_sig);
        
    // assign IMU noise variances
    // inputs to the system are 3 delta angles and 3 delta velocities
    scalar_t daxVar, dayVar, dazVar;
    gyro_noise = constrain(gyro_noise, 0.0, 1.0);
    daxVar = dayVar = dazVar = sq(dt * gyro_noise); // gyro prediction variance TODO get variance from sensor
    
    // intermediate calculations
    scalar_t SF[9];
    SF[0] = day/2 - day_b/2;
    SF[1] = daz/2 - daz_b/2;
    SF[2] = dax/2 - dax_b/2;
    SF[3] = dax_b/2 - dax/2;
    SF[4] = daz_b/2 - daz/2;
    SF[5] = day_b/2 - day/2;
    SF[6] = q1/2;
    SF[7] = q2/2;
    SF[8] = q3/2;

    scalar_t SG[2];
    SG[0] = q0/2;
    SG[1] = q2/2;

    scalar_t SQ[8];
    SQ[0] = (dayVar*q1*SG[0])/2 - (daxVar*q3*SG[1])/2 - (dazVar*q1*SG[0])/2;
    SQ[1] = (daxVar*q3*SG[0])/2 - (dayVar*q3*SG[0])/2 - (dazVar*q1*SG[1])/2;
    SQ[2] = (daxVar*q1*SG[1])/2 - (dazVar*q3*SG[0])/2 - (dayVar*q1*q2)/4;
    SQ[3] = (dayVar*q2*q3)/4 - (dazVar*q3*SG[1])/2 - (daxVar*q1*SG[0])/2;
    SQ[4] = (dazVar*q1*q3)/4 - (daxVar*q1*q3)/4 - (dayVar*q2*SG[0])/2;
    SQ[5] = dazVar*SG[0]*SG[1] - daxVar*SG[0]*SG[1] - (dayVar*q1*q3)/4;
    SQ[6] = sq(SG[0]);
    SQ[7] = sq(q1);

    scalar_t SPP[14];
    SPP[0] = SF[8];
    SPP[1] = SF[7];
    SPP[2] = P_[3][4] + P_[0][4]*SF[1] + P_[1][4]*SF[0] + P_[2][4]*SF[3] - P_[5][4]*SF[6] + P_[4][4]*SPP[1] - (P_[6][4]*q0)/2;
    SPP[3] = P_[3][6] + P_[0][6]*SF[1] + P_[1][6]*SF[0] + P_[2][6]*SF[3] - P_[5][6]*SF[6] + P_[4][6]*SPP[1] - (P_[6][6]*q0)/2;
    SPP[4] = P_[3][5] + P_[0][5]*SF[1] + P_[1][5]*SF[0] + P_[2][5]*SF[3] - P_[5][5]*SF[6] + P_[4][5]*SPP[1] - (P_[6][5]*q0)/2;
    SPP[5] = P_[2][4] + P_[0][4]*SF[0] + P_[1][4]*SF[4] + P_[3][4]*SF[2] + P_[6][4]*SF[6] - P_[4][4]*SPP[0] - (P_[5][4]*q0)/2;
    SPP[6] = P_[2][6] + P_[0][6]*SF[0] + P_[1][6]*SF[4] + P_[3][6]*SF[2] + P_[6][6]*SF[6] - P_[4][6]*SPP[0] - (P_[5][6]*q0)/2;
    SPP[7] = P_[2][5] + P_[0][5]*SF[0] + P_[1][5]*SF[4] + P_[3][5]*SF[2] + P_[6][5]*SF[6] - P_[4][5]*SPP[0] - (P_[5][5]*q0)/2;
    SPP[8] = P_[1][4] + P_[0][4]*SF[2] + P_[2][4]*SF[1] + P_[3][4]*SF[5] + P_[5][4]*SPP[0] - P_[6][4]*SPP[1] - (P_[4][4]*q0)/2;
    SPP[9] = P_[1][6] + P_[0][6]*SF[2] + P_[2][6]*SF[1] + P_[3][6]*SF[5] + P_[5][6]*SPP[0] - P_[6][6]*SPP[1] - (P_[4][6]*q0)/2;
    SPP[10] = P_[1][5] + P_[0][5]*SF[2] + P_[2][5]*SF[1] + P_[3][5]*SF[5] + P_[5][5]*SPP[0] - P_[6][5]*SPP[1] - (P_[4][5]*q0)/2;
    SPP[11] = P_[0][4] + P_[1][4]*SF[3] + P_[2][4]*SF[5] + P_[3][4]*SF[4] + P_[4][4]*SF[6] + P_[5][4]*SPP[1] + P_[6][4]*SPP[0];
    SPP[12] = P_[0][6] + P_[1][6]*SF[3] + P_[2][6]*SF[5] + P_[3][6]*SF[4] + P_[4][6]*SF[6] + P_[5][6]*SPP[1] + P_[6][6]*SPP[0];
    SPP[13] = P_[0][5] + P_[1][5]*SF[3] + P_[2][5]*SF[5] + P_[3][5]*SF[4] + P_[4][5]*SF[6] + P_[5][5]*SPP[1] + P_[6][5]*SPP[0];
    
    // covariance update
    // calculate variances and upper diagonal covariances for quaternion and gyro bias states
    scalar_t nextP[k_num_states_][k_num_states_];
    nextP[0][0] = P_[0][0] + P_[1][0]*SF[3] + P_[2][0]*SF[5] + P_[3][0]*SF[4] + P_[4][0]*SF[6] + P_[5][0]*SPP[1] + P_[6][0]*SPP[0] + SF[6]*SPP[11] + SPP[0]*SPP[12] + SPP[1]*SPP[13] + (daxVar*SQ[7])/4 + SF[3]*(P_[0][1] + P_[1][1]*SF[3] + P_[2][1]*SF[5] + P_[3][1]*SF[4] + P_[4][1]*SF[6] + P_[5][1]*SPP[1] + P_[6][1]*SPP[0]) + SF[5]*(P_[0][2] + P_[1][2]*SF[3] + P_[2][2]*SF[5] + P_[3][2]*SF[4] + P_[4][2]*SF[6] + P_[5][2]*SPP[1] + P_[6][2]*SPP[0]) + SF[4]*(P_[0][3] + P_[1][3]*SF[3] + P_[2][3]*SF[5] + P_[3][3]*SF[4] + P_[4][3]*SF[6] + P_[5][3]*SPP[1] + P_[6][3]*SPP[0]) + (dayVar*sq(q2))/4 + (dazVar*sq(q3))/4;
    nextP[0][1] = P_[0][1] + SQ[3] + P_[1][1]*SF[3] + P_[2][1]*SF[5] + P_[3][1]*SF[4] + P_[4][1]*SF[6] + P_[5][1]*SPP[1] + P_[6][1]*SPP[0] + SPP[0]*SPP[13] - SPP[1]*SPP[12] - (q0*SPP[11])/2 + SF[2]*(P_[0][0] + P_[1][0]*SF[3] + P_[2][0]*SF[5] + P_[3][0]*SF[4] + P_[4][0]*SF[6] + P_[5][0]*SPP[1] + P_[6][0]*SPP[0]) + SF[1]*(P_[0][2] + P_[1][2]*SF[3] + P_[2][2]*SF[5] + P_[3][2]*SF[4] + P_[4][2]*SF[6] + P_[5][2]*SPP[1] + P_[6][2]*SPP[0]) + SF[5]*(P_[0][3] + P_[1][3]*SF[3] + P_[2][3]*SF[5] + P_[3][3]*SF[4] + P_[4][3]*SF[6] + P_[5][3]*SPP[1] + P_[6][3]*SPP[0]);
    nextP[1][1] = P_[1][1] + P_[0][1]*SF[2] + P_[2][1]*SF[1] + P_[3][1]*SF[5] + P_[5][1]*SPP[0] - P_[6][1]*SPP[1] + SPP[0]*SPP[10] - SPP[1]*SPP[9] + daxVar*SQ[6] - (P_[4][1]*q0)/2 - (q0*SPP[8])/2 + dazVar*sq(SG[1]) + SF[2]*(P_[1][0] + P_[0][0]*SF[2] + P_[2][0]*SF[1] + P_[3][0]*SF[5] + P_[5][0]*SPP[0] - P_[6][0]*SPP[1] - (P_[4][0]*q0)/2) + SF[1]*(P_[1][2] + P_[0][2]*SF[2] + P_[2][2]*SF[1] + P_[3][2]*SF[5] + P_[5][2]*SPP[0] - P_[6][2]*SPP[1] - (P_[4][2]*q0)/2) + SF[5]*(P_[1][3] + P_[0][3]*SF[2] + P_[2][3]*SF[1] + P_[3][3]*SF[5] + P_[5][3]*SPP[0] - P_[6][3]*SPP[1] - (P_[4][3]*q0)/2) + (dayVar*sq(q3))/4;
    nextP[0][2] = P_[0][2] + SQ[4] + P_[1][2]*SF[3] + P_[2][2]*SF[5] + P_[3][2]*SF[4] + P_[4][2]*SF[6] + P_[5][2]*SPP[1] + P_[6][2]*SPP[0] + SF[6]*SPP[12] - SPP[0]*SPP[11] - (q0*SPP[13])/2 + SF[0]*(P_[0][0] + P_[1][0]*SF[3] + P_[2][0]*SF[5] + P_[3][0]*SF[4] + P_[4][0]*SF[6] + P_[5][0]*SPP[1] + P_[6][0]*SPP[0]) + SF[4]*(P_[0][1] + P_[1][1]*SF[3] + P_[2][1]*SF[5] + P_[3][1]*SF[4] + P_[4][1]*SF[6] + P_[5][1]*SPP[1] + P_[6][1]*SPP[0]) + SF[2]*(P_[0][3] + P_[1][3]*SF[3] + P_[2][3]*SF[5] + P_[3][3]*SF[4] + P_[4][3]*SF[6] + P_[5][3]*SPP[1] + P_[6][3]*SPP[0]);
    nextP[1][2] = P_[1][2] + SQ[1] + P_[0][2]*SF[2] + P_[2][2]*SF[1] + P_[3][2]*SF[5] + P_[5][2]*SPP[0] - P_[6][2]*SPP[1] + SF[6]*SPP[9] - SPP[0]*SPP[8] - (P_[4][2]*q0)/2 - (q0*SPP[10])/2 + SF[0]*(P_[1][0] + P_[0][0]*SF[2] + P_[2][0]*SF[1] + P_[3][0]*SF[5] + P_[5][0]*SPP[0] - P_[6][0]*SPP[1] - (P_[4][0]*q0)/2) + SF[4]*(P_[1][1] + P_[0][1]*SF[2] + P_[2][1]*SF[1] + P_[3][1]*SF[5] + P_[5][1]*SPP[0] - P_[6][1]*SPP[1] - (P_[4][1]*q0)/2) + SF[2]*(P_[1][3] + P_[0][3]*SF[2] + P_[2][3]*SF[1] + P_[3][3]*SF[5] + P_[5][3]*SPP[0] - P_[6][3]*SPP[1] - (P_[4][3]*q0)/2);
    nextP[2][2] = P_[2][2] + P_[0][2]*SF[0] + P_[1][2]*SF[4] + P_[3][2]*SF[2] + P_[6][2]*SF[6] - P_[4][2]*SPP[0] + SF[6]*SPP[6] - SPP[0]*SPP[5] + dayVar*SQ[6] + (dazVar*SQ[7])/4 - (P_[5][2]*q0)/2 - (q0*SPP[7])/2 + SF[0]*(P_[2][0] + P_[0][0]*SF[0] + P_[1][0]*SF[4] + P_[3][0]*SF[2] + P_[6][0]*SF[6] - P_[4][0]*SPP[0] - (P_[5][0]*q0)/2) + SF[4]*(P_[2][1] + P_[0][1]*SF[0] + P_[1][1]*SF[4] + P_[3][1]*SF[2] + P_[6][1]*SF[6] - P_[4][1]*SPP[0] - (P_[5][1]*q0)/2) + SF[2]*(P_[2][3] + P_[0][3]*SF[0] + P_[1][3]*SF[4] + P_[3][3]*SF[2] + P_[6][3]*SF[6] - P_[4][3]*SPP[0] - (P_[5][3]*q0)/2) + (daxVar*sq(q3))/4;
    nextP[0][3] = P_[0][3] + SQ[2] + P_[1][3]*SF[3] + P_[2][3]*SF[5] + P_[3][3]*SF[4] + P_[4][3]*SF[6] + P_[5][3]*SPP[1] + P_[6][3]*SPP[0] - SF[6]*SPP[13] + SPP[1]*SPP[11] - (q0*SPP[12])/2 + SF[1]*(P_[0][0] + P_[1][0]*SF[3] + P_[2][0]*SF[5] + P_[3][0]*SF[4] + P_[4][0]*SF[6] + P_[5][0]*SPP[1] + P_[6][0]*SPP[0]) + SF[0]*(P_[0][1] + P_[1][1]*SF[3] + P_[2][1]*SF[5] + P_[3][1]*SF[4] + P_[4][1]*SF[6] + P_[5][1]*SPP[1] + P_[6][1]*SPP[0]) + SF[3]*(P_[0][2] + P_[1][2]*SF[3] + P_[2][2]*SF[5] + P_[3][2]*SF[4] + P_[4][2]*SF[6] + P_[5][2]*SPP[1] + P_[6][2]*SPP[0]);
    nextP[1][3] = P_[1][3] + SQ[5] + P_[0][3]*SF[2] + P_[2][3]*SF[1] + P_[3][3]*SF[5] + P_[5][3]*SPP[0] - P_[6][3]*SPP[1] - SF[6]*SPP[10] + SPP[1]*SPP[8] - (P_[4][3]*q0)/2 - (q0*SPP[9])/2 + SF[1]*(P_[1][0] + P_[0][0]*SF[2] + P_[2][0]*SF[1] + P_[3][0]*SF[5] + P_[5][0]*SPP[0] - P_[6][0]*SPP[1] - (P_[4][0]*q0)/2) + SF[0]*(P_[1][1] + P_[0][1]*SF[2] + P_[2][1]*SF[1] + P_[3][1]*SF[5] + P_[5][1]*SPP[0] - P_[6][1]*SPP[1] - (P_[4][1]*q0)/2) + SF[3]*(P_[1][2] + P_[0][2]*SF[2] + P_[2][2]*SF[1] + P_[3][2]*SF[5] + P_[5][2]*SPP[0] - P_[6][2]*SPP[1] - (P_[4][2]*q0)/2);
    nextP[2][3] = P_[2][3] + SQ[0] + P_[0][3]*SF[0] + P_[1][3]*SF[4] + P_[3][3]*SF[2] + P_[6][3]*SF[6] - P_[4][3]*SPP[0] - SF[6]*SPP[7] + SPP[1]*SPP[5] - (P_[5][3]*q0)/2 - (q0*SPP[6])/2 + SF[1]*(P_[2][0] + P_[0][0]*SF[0] + P_[1][0]*SF[4] + P_[3][0]*SF[2] + P_[6][0]*SF[6] - P_[4][0]*SPP[0] - (P_[5][0]*q0)/2) + SF[0]*(P_[2][1] + P_[0][1]*SF[0] + P_[1][1]*SF[4] + P_[3][1]*SF[2] + P_[6][1]*SF[6] - P_[4][1]*SPP[0] - (P_[5][1]*q0)/2) + SF[3]*(P_[2][2] + P_[0][2]*SF[0] + P_[1][2]*SF[4] + P_[3][2]*SF[2] + P_[6][2]*SF[6] - P_[4][2]*SPP[0] - (P_[5][2]*q0)/2);
    nextP[3][3] = P_[3][3] + P_[0][3]*SF[1] + P_[1][3]*SF[0] + P_[2][3]*SF[3] - P_[5][3]*SF[6] + P_[4][3]*SPP[1] - SF[6]*SPP[4] + SPP[1]*SPP[2] + (dayVar*SQ[7])/4 + dazVar*SQ[6] - (P_[6][3]*q0)/2 - (q0*SPP[3])/2 + daxVar*sq(SG[1]) + SF[1]*(P_[3][0] + P_[0][0]*SF[1] + P_[1][0]*SF[0] + P_[2][0]*SF[3] - P_[5][0]*SF[6] + P_[4][0]*SPP[1] - (P_[6][0]*q0)/2) + SF[0]*(P_[3][1] + P_[0][1]*SF[1] + P_[1][1]*SF[0] + P_[2][1]*SF[3] - P_[5][1]*SF[6] + P_[4][1]*SPP[1] - (P_[6][1]*q0)/2) + SF[3]*(P_[3][2] + P_[0][2]*SF[1] + P_[1][2]*SF[0] + P_[2][2]*SF[3] - P_[5][2]*SF[6] + P_[4][2]*SPP[1] - (P_[6][2]*q0)/2);
    nextP[0][4] = SPP[11];
    nextP[1][4] = SPP[8];
    nextP[2][4] = SPP[5];
    nextP[3][4] = SPP[2];
    nextP[4][4] = P_[4][4];
    nextP[0][5] = SPP[13];
    nextP[1][5] = SPP[10];
    nextP[2][5] = SPP[7];
    nextP[3][5] = SPP[4];
    nextP[4][5] = P_[4][5];
    nextP[5][5] = P_[5][5];
    nextP[0][6] = SPP[12];
    nextP[1][6] = SPP[9];
    nextP[2][6] = SPP[6];
    nextP[3][6] = SPP[3];
    nextP[4][6] = P_[4][6];
    nextP[5][6] = P_[5][6];
    nextP[6][6] = P_[6][6];
    
    // add process noise that is not from the IMU
    for (unsigned i = 0; i < k_num_states_; i++) {
      nextP[i][i] += process_noise[i];
    }
    
    // covariance matrix is symmetrical, so copy upper half to lower half
    for (unsigned row = 1; row < k_num_states_; row++) {
      for (unsigned column = 0 ; column < row; column++) {
	      P_[row][column] = P_[column][row] = nextP[column][row];
      }
    }

    // copy variances (diagonals)
    for (unsigned i = 0; i < k_num_states_; i++) {
      P_[i][i] = nextP[i][i];
    }
    
    curr_time_sec += dt;
    
    if(no_measurement) {
      
      scalar_t R[3] = {}; // observation variances for [PN,PE,PD]
      scalar_t gate_size[3] = {}; // innovation consistency check gate sizes for [PN,PE,PD] observations
      scalar_t Kfusion[k_num_states_] = {}; // Kalman gain vector for any single observation - sequential fusion is used
      scalar_t att_innov[3] = {}; //
      scalar_t att_innov_var[3] = {}; //
      scalar_t att_test_ratio[3] = {}; //
      bool innov_check_pass_map[3] = {}; //
      
      R[0] = 0.5f;
      att_innov[0] = yaw - last_known_yaw;
      att_innov[1] = pitch - last_known_pitch;

      // glitch protection is not required so set gate to a large value
      gate_size[0] = 100.0f;
      R[0] = R[0] * R[0];

      // copy North axis values to East axis
      R[1] = R[0];
      gate_size[1] = gate_size[0];
      
      att_innov[2] = roll - last_known_roll;
      // observation variance - user parameter defined
      R[2] = 2.0f;
      R[2] = R[2] * R[2];
      // innovation gate size
      gate_size[2] = 5.0f;
      
      // calculate innovation test ratios
      for (unsigned obs_index = 0; obs_index < 3; obs_index++) {
      // compute the innovation variance SK = HPH + R
      unsigned state_index = obs_index + 4;	// we start with gyro_bias and this is the 4. state
      att_innov_var[obs_index] = P_[state_index][state_index] + R[obs_index];
        // Compute the ratio of innovation to gate size
        att_test_ratio[obs_index] = sq(att_innov[obs_index]) / (sq(gate_size[obs_index]) * att_innov_var[obs_index]);
      }
      
      bool att_check_pass = ((att_test_ratio[0] <= 1.0f) && (att_test_ratio[1] <= 1.0f) && (att_test_ratio[2] <= 1.0f));
      innov_check_pass_map[0] = innov_check_pass_map[1] = innov_check_pass_map[2] = att_check_pass;
      
      for (unsigned obs_index = 0; obs_index < 3; obs_index++) {
	      // skip fusion if not requested or checks have failed
	      if (!innov_check_pass_map[obs_index]) {
	        continue;
	      }

        unsigned state_index = obs_index + 4;	// we start with gyro_bias and this is the 4. state

        // calculate kalman gain K = PHS, where S = 1/innovation variance
        for (int row = 0; row < k_num_states_; row++) {
          Kfusion[row] = P_[row][state_index] / att_innov_var[obs_index];
        }

        // update covarinace matrix via Pnew = (I - KH)P
        scalar_t KHP[k_num_states_][k_num_states_];
        for (unsigned row = 0; row < k_num_states_; row++) {
          for (unsigned column = 0; column < k_num_states_; column++) {
            KHP[row][column] = Kfusion[row] * P_[state_index][column];
          }
        }

        // if the covariance correction will result in a negative variance, then
        // the covariance marix is unhealthy and must be corrected
        bool healthy = true;
        for (int i = 0; i < k_num_states_; i++) {
          if (P_[i][i] < KHP[i][i]) {
            // zero rows and columns
            zeroRows(P_,i,i);
            zeroCols(P_,i,i);

            //flag as unhealthy
            healthy = false;
          }
        }

        // only apply covariance and state corrrections if healthy
        if (healthy) {
          // apply the covariance corrections
          for (unsigned row = 0; row < k_num_states_; row++) {
            for (unsigned column = 0; column < k_num_states_; column++) {
              P_[row][column] = P_[row][column] - KHP[row][column];
            }
          }

          // correct the covariance marix for gross errors
          fixCovarianceErrors(dt);

          // apply the state corrections
          fuse(Kfusion, att_innov[obs_index]);
        }
      }
      R_to_earth = quat_to_invrotmat(state_.quat_nominal);
      constrainStates(dt);
      last_known_yaw = atan2f(-R_to_earth(0, 1), R_to_earth(1, 1)); // first rotation (yaw)
      last_known_pitch = atan2f(-R_to_earth(2, 0), R_to_earth(2, 2)); // third rotation (pitch)
      last_known_roll = asinf(R_to_earth(2, 1)); // second rotation (roll)
    }
  }
  
  void ESKF::update(const ESKF::quat& q, scalar_t dt) {
    if(!firstPredict) {
      no_measurement = false; // we have measurement now
      return; // if there hasn't been any prediction yet then there should be no correction
    }
    // q here is rotation from enu to ros body
    updateYaw(q, dt);
  }
  
  void ESKF::updateYaw(const quat& q, scalar_t dt) {
    // assign intermediate state variables
    scalar_t q0 = state_.quat_nominal.w();
    scalar_t q1 = state_.quat_nominal.x();
    scalar_t q2 = state_.quat_nominal.y();
    scalar_t q3 = state_.quat_nominal.z();
    
    scalar_t predicted_hdg, measured_hdg;
    scalar_t H_YAW[4];
    
    q_nb = (q_rb.conjugate() * q.conjugate() * q_ne.conjugate()).conjugate();
    q_nb.normalize();
    //0.7071 0 0 0.7071
    //std::cout << "q_nb = " << q_nb.w() << " " << q_nb.x() << " " << q_nb.y() << " " << q_nb.z() << std::endl;
    //std::cout << "q_er = " << std::endl << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << std::endl;
    
    //std::cout << "n2b = " << std::endl << q_nb.toRotationMatrix() << std::endl;
    //std::cout << "e2r = " << std::endl << q.conjugate().toRotationMatrix() << std::endl;
    // update transformation matrix from body to world frame
    //std::cout << "quat_nominal = " << state_.quat_nominal.w() << " " << state_.quat_nominal.x() << " " << state_.quat_nominal.y() << " " << state_.quat_nominal.z() << std::endl;
    mat3 R_to_earth = quat_to_invrotmat(state_.quat_nominal);
    //std::cout << "R_to_earth = " << std::endl << R_to_earth << std::endl;
	  {
      // calculate observaton jacobian when we are observing a rotation in a 312 sequence
      scalar_t t9 = q0*q3;
      scalar_t t10 = q1*q2;
      scalar_t t2 = t9-t10;
      scalar_t t3 = q0*q0;
      scalar_t t4 = q1*q1;
      scalar_t t5 = q2*q2;
      scalar_t t6 = q3*q3;
      scalar_t t7 = t3-t4+t5-t6;
      scalar_t t8 = t7*t7;
      if (t8 > 1e-6f) {
        t8 = 1.0f/t8;
      } else {
        return;
      }
      scalar_t t11 = t2*t2;
      scalar_t t12 = t8*t11*4.0f;
      scalar_t t13 = t12+1.0f;
      scalar_t t14;
      if (fabsf(t13) > 1e-6f) {
        t14 = 1.0f/t13;
      } else {
        return;
      }

      H_YAW[0] = t8*t14*(q3*t3+q3*t4-q3*t5+q3*t6-q0*q1*q2*2.0f)*-2.0f;
      H_YAW[1] = t8*t14*(q2*t3+q2*t4+q2*t5-q2*t6-q0*q1*q3*2.0f)*-2.0f;
      H_YAW[2] = t8*t14*(-q1*t3+q1*t4+q1*t5+q1*t6-q0*q2*q3*2.0f)*2.0f;
      H_YAW[3] = t8*t14*(q0*t3-q0*t4+q0*t5+q0*t6-q1*q2*q3*2.0f)*2.0f;

      /* Calculate the 312 sequence euler angles that rotate from earth to body frame
       * Derived from https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
       * Body to nav frame transformation using a yaw-roll-pitch rotation sequence is given by:
       *
      [ cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw), -cos(roll)*sin(yaw), cos(yaw)*sin(pitch) + cos(pitch)*sin(roll)*sin(yaw)]
      [ cos(pitch)*sin(yaw) + cos(yaw)*sin(pitch)*sin(roll),  cos(roll)*cos(yaw), sin(pitch)*sin(yaw) - cos(pitch)*cos(yaw)*sin(roll)]
      [                               -cos(roll)*sin(pitch),           sin(roll),                                cos(pitch)*cos(roll)]
      */
      scalar_t yaw = atan2f(-R_to_earth(0, 1), R_to_earth(1, 1)); // first rotation (yaw)
      //float roll = asinf(R_to_earth(2, 1)); // second rotation (roll)
      //float pitch = atan2f(-R_to_earth(2, 0), R_to_earth(2, 2)); // third rotation (pitch)

      predicted_hdg = yaw; // we will need the predicted heading to calculate the innovation
      //std::cout << "predicted_hdg = " << predicted_hdg*57.3 << std::endl;
			// calculate the yaw angle for a 312 sequence
			// Values from yaw_input_312.c file produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
			float Tbn_0_1_neg = 2.0f * (q_nb.w() * q_nb.z() - q_nb.x() * q_nb.y());
			float Tbn_1_1 = sq(q_nb.w()) - sq(q_nb.x()) + sq(q_nb.y()) - sq(q_nb.z());
			measured_hdg = atan2f(Tbn_0_1_neg, Tbn_1_1);
      //std::cout << "measured_hdg = " << measured_hdg*57.3 << std::endl;
	  }

    scalar_t R_YAW = sq(fmaxf(yaw_err, 1.0e-2f));
    // Calculate innovation variance and Kalman gains, taking advantage of the fact that only the first 3 elements in H are non zero
    // calculate the innovaton variance
    scalar_t PH[4];
    scalar_t heading_innov_var = R_YAW;
    for (unsigned row = 0; row <= 3; row++) {
		  PH[row] = 0.0f;
      
      for (uint8_t col = 0; col <= 3; col++) {
			  PH[row] += P_[row][col] * H_YAW[col];
		  }
  		heading_innov_var += H_YAW[row] * PH[row];
	  }

	  scalar_t heading_innov_var_inv;

    // check if the innovation variance calculation is badly conditioned
    if (heading_innov_var >= R_YAW) {
      // the innovation variance contribution from the state covariances is not negative, no fault
      heading_innov_var_inv = 1.0f / heading_innov_var;
    } else {
      // the innovation variance contribution from the state covariances is negative which means the covariance matrix is badly conditioned
      // we reinitialise the covariance matrix and abort this fusion step
      initialiseCovariance(dt);
      printf("EKF mag yaw fusion numerical error - covariance reset\n");
      return;
    }

    // calculate the Kalman gains
    // only calculate gains for states we are using
    scalar_t Kfusion[k_num_states_] = {};

    for (uint8_t row = 0; row < k_num_states_; row++) {
      Kfusion[row] = 0.0f;

      for (uint8_t col = 0; col <= 3; col++) {
        Kfusion[row] += P_[row][col] * H_YAW[col];
      }

      Kfusion[row] *= heading_innov_var_inv;
    }

	  // wrap the heading to the interval between +-pi
	  measured_hdg = wrap_pi(measured_hdg);

	  // calculate the innovation
	  scalar_t heading_innov = predicted_hdg - measured_hdg;

	  // wrap the innovation to the interval between +-pi
	  heading_innov = wrap_pi(heading_innov);

	  // innovation test ratio
	  scalar_t yaw_test_ratio = sq(heading_innov) / (sq(heading_innov_gate) * heading_innov_var);
	  
    // set the vision yaw unhealthy if the test fails
	  if (yaw_test_ratio > 1.0f) {
      printf("!!!!!\n");
      scalar_t gate_limit = sqrtf(sq(heading_innov_gate) * heading_innov_var);
			heading_innov = constrain(heading_innov, -gate_limit, gate_limit);
		  //return;
		} /*else {
			// constrain the innovation to the maximum set by the gate
			scalar_t gate_limit = sqrtf(sq(heading_innov_gate) * heading_innov_var);
			heading_innov = constrain(heading_innov, -gate_limit, gate_limit);
		}*/
	
    // apply covariance correction via P_new = (I -K*H)*P
    // first calculate expression for KHP
    // then calculate P - KHP
    float KHP[k_num_states_][k_num_states_];
    float KH[4];
    for (unsigned row = 0; row < k_num_states_; row++) {

      KH[0] = Kfusion[row] * H_YAW[0];
      KH[1] = Kfusion[row] * H_YAW[1];
      KH[2] = Kfusion[row] * H_YAW[2];
      KH[3] = Kfusion[row] * H_YAW[3];

      for (unsigned column = 0; column < k_num_states_; column++) {
        float tmp = KH[0] * P_[0][column];
        tmp += KH[1] * P_[1][column];
        tmp += KH[2] * P_[2][column];
        tmp += KH[3] * P_[3][column];
        KHP[row][column] = tmp;
      }
    }

    // if the covariance correction will result in a negative variance, then
    // the covariance marix is unhealthy and must be corrected
    bool healthy = true;
    for (int i = 0; i < k_num_states_; i++) {
      if (P_[i][i] < KHP[i][i]) {
        // zero rows and columns
        zeroRows(P_,i,i);
        zeroCols(P_,i,i);

        //flag as unhealthy
        healthy = false;
      }
    }

    // only apply covariance and state corrrections if healthy
    if (healthy) {
      // apply the covariance corrections
      for (unsigned row = 0; row < k_num_states_; row++) {
        for (unsigned column = 0; column < k_num_states_; column++) {
          P_[row][column] = P_[row][column] - KHP[row][column];
        }
      }

      // correct the covariance marix for gross errors
      fixCovarianceErrors(dt);

      // apply the state corrections
      fuse(Kfusion, heading_innov);
    }
  }
  
  // calculate the inverse rotation matrix from a quaternion rotation
  ESKF::mat3 ESKF::quat_to_invrotmat(const quat &q) {
    scalar_t q00 = q.w() * q.w();
    scalar_t q11 = q.x() * q.x();
    scalar_t q22 = q.y() * q.y();
    scalar_t q33 = q.z() * q.z();
    scalar_t q01 = q.w() * q.x();
    scalar_t q02 = q.w() * q.y();
    scalar_t q03 = q.w() * q.z();
    scalar_t q12 = q.x() * q.y();
    scalar_t q13 = q.x() * q.z();
    scalar_t q23 = q.y() * q.z();

    mat3 dcm;
    dcm(0, 0) = q00 + q11 - q22 - q33;
    dcm(1, 1) = q00 - q11 + q22 - q33;
    dcm(2, 2) = q00 - q11 - q22 + q33;
    dcm(0, 1) = 2.0f * (q12 - q03);
    dcm(0, 2) = 2.0f * (q13 + q02);
    dcm(1, 0) = 2.0f * (q12 + q03);
    dcm(1, 2) = 2.0f * (q23 - q01);
    dcm(2, 0) = 2.0f * (q13 - q02);
    dcm(2, 1) = 2.0f * (q23 + q01);
    
	  return dcm;
  }
  
  void ESKF::fixCovarianceErrors(scalar_t dt) {
    // NOTE: This limiting is a last resort and should not be relied on
    // TODO: Split covariance prediction into separate F*P*transpose(F) and Q contributions
    // and set corresponding entries in Q to zero when states exceed 50% of the limit
    // Covariance diagonal limits. Use same values for states which
    // belong to the same group (e.g. vel_x, vel_y, vel_z)
    scalar_t P_lim[2] = {};
    P_lim[0] = 1.0f;		// quaternion max var
    P_lim[1] = 1.0f;		// gyro bias max var
    
    for (int i = 0; i <= 3; i++) {
		  // quaternion states
		  P_[i][i] = constrain(P_[i][i], 0.0, P_lim[0]);
	  }
    
	  for (int i = 4; i <= 6; i++) {
		  // gyro bias states
		  P_[i][i] = constrain(P_[i][i], 0.0, P_lim[1]);
	  }
    
	  // force symmetry on the quaternion and gyro bias state covariances
	  makeSymmetrical(P_, 0, k_num_states_ - 1);
  }
  
  // This function forces the covariance matrix to be symmetric
  void ESKF::makeSymmetrical(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last) {
    for (unsigned row = first; row <= last; row++) {
      for (unsigned column = 0; column < row; column++) {
        float tmp = (cov_mat[row][column] + cov_mat[column][row]) / 2;
        cov_mat[row][column] = tmp;
        cov_mat[column][row] = tmp;
      }
    }
  }
  
  // fuse measurement
  void ESKF::fuse(scalar_t *K, scalar_t innovation) {
        
    state_.quat_nominal.w() = state_.quat_nominal.w() - K[0] * innovation;
    state_.quat_nominal.x() = state_.quat_nominal.x() - K[1] * innovation;
    state_.quat_nominal.y() = state_.quat_nominal.y() - K[2] * innovation;
    state_.quat_nominal.z() = state_.quat_nominal.z() - K[3] * innovation;
    
    state_.quat_nominal.normalize();

    for (unsigned i = 0; i < 3; i++) {
      state_.gyro_bias(i) = state_.gyro_bias(i) - K[i + 4] * innovation;
    }
  }
  
  ESKF::vec3 ESKF::getRPY(const mat3& R) {
    vec3 rpy;
	  scalar_t sth = -R(2, 0);
    if (sth > 1) {
      sth = 1;
    } else if (sth < -1) {
      sth = -1;
    }
    
    const scalar_t theta = std::asin(sth);
    const scalar_t cth = std::sqrt(1 - sth*sth);
    
    scalar_t phi, psi;
    if (cth < static_cast<scalar_t>(1.0e-6)) {
      phi = std::atan2(R(0, 1), R(1, 1));
      psi = 0;
    } else {
      phi = std::atan2(R(2, 1), R(2, 2));
      psi = std::atan2(R(1, 0), R(0, 0));
    }
    
    rpy[0] = phi;    //  x, [-pi,pi]
    rpy[1] = theta;  //  y, [-pi/2,pi/2]
    rpy[2] = psi;    //  z, [-pi,pi]
    return rpy;
  }
    
  const ESKF::quat& ESKF::getQuat() const { 
    return state_.quat_nominal; 
  }
  
  // initialise the quaternion covariances using rotation vector variances
  void ESKF::initialiseQuatCovariances(const vec3& rot_vec_var) {
    // calculate an equivalent rotation vector from the quaternion
    scalar_t q0,q1,q2,q3;
    if (state_.quat_nominal.w() >= 0.0f) {
      q0 = state_.quat_nominal.w();
      q1 = state_.quat_nominal.x();
      q2 = state_.quat_nominal.y();
      q3 = state_.quat_nominal.z();
    } else {
      q0 = -state_.quat_nominal.w();
      q1 = -state_.quat_nominal.x();
      q2 = -state_.quat_nominal.y();
      q3 = -state_.quat_nominal.z();
    }
    scalar_t delta = 2.0f*acosf(q0);
    scalar_t scaler = (delta/sinf(delta*0.5f));
    scalar_t rotX = scaler*q1;
    scalar_t rotY = scaler*q2;
    scalar_t rotZ = scaler*q3;

    // autocode generated using matlab symbolic toolbox
    scalar_t t2 = rotX*rotX;
    scalar_t t4 = rotY*rotY;
    scalar_t t5 = rotZ*rotZ;
    scalar_t t6 = t2+t4+t5;
    if (t6 > 1e-9f) {
      scalar_t t7 = sqrtf(t6);
      scalar_t t8 = t7*0.5f;
      scalar_t t3 = sinf(t8);
      scalar_t t9 = t3*t3;
      scalar_t t10 = 1.0f/t6;
      scalar_t t11 = 1.0f/sqrtf(t6);
      scalar_t t12 = cosf(t8);
      scalar_t t13 = 1.0f/powf(t6,1.5f);
      scalar_t t14 = t3*t11;
      scalar_t t15 = rotX*rotY*t3*t13;
      scalar_t t16 = rotX*rotZ*t3*t13;
      scalar_t t17 = rotY*rotZ*t3*t13;
      scalar_t t18 = t2*t10*t12*0.5f;
      scalar_t t27 = t2*t3*t13;
      scalar_t t19 = t14+t18-t27;
      scalar_t t23 = rotX*rotY*t10*t12*0.5f;
      scalar_t t28 = t15-t23;
      scalar_t t20 = rotY*rot_vec_var(1)*t3*t11*t28*0.5f;
      scalar_t t25 = rotX*rotZ*t10*t12*0.5f;
      scalar_t t31 = t16-t25;
      scalar_t t21 = rotZ*rot_vec_var(2)*t3*t11*t31*0.5f;
      scalar_t t22 = t20+t21-rotX*rot_vec_var(0)*t3*t11*t19*0.5f;
      scalar_t t24 = t15-t23;
      scalar_t t26 = t16-t25;
      scalar_t t29 = t4*t10*t12*0.5f;
      scalar_t t34 = t3*t4*t13;
      scalar_t t30 = t14+t29-t34;
      scalar_t t32 = t5*t10*t12*0.5f;
      scalar_t t40 = t3*t5*t13;
      scalar_t t33 = t14+t32-t40;
      scalar_t t36 = rotY*rotZ*t10*t12*0.5f;
      scalar_t t39 = t17-t36;
      scalar_t t35 = rotZ*rot_vec_var(2)*t3*t11*t39*0.5f;
      scalar_t t37 = t15-t23;
      scalar_t t38 = t17-t36;
      scalar_t t41 = rot_vec_var(0)*(t15-t23)*(t16-t25);
      scalar_t t42 = t41-rot_vec_var(1)*t30*t39-rot_vec_var(2)*t33*t39;
      scalar_t t43 = t16-t25;
      scalar_t t44 = t17-t36;

      // zero all the quaternion covariances
      zeroRows(P_,0,3);
      zeroCols(P_,0,3);

      // Update the quaternion internal covariances using auto-code generated using matlab symbolic toolbox
      P_[0][0] = rot_vec_var(0)*t2*t9*t10*0.25f+rot_vec_var(1)*t4*t9*t10*0.25f+rot_vec_var(2)*t5*t9*t10*0.25f;
      P_[0][1] = t22;
      P_[0][2] = t35+rotX*rot_vec_var(0)*t3*t11*(t15-rotX*rotY*t10*t12*0.5f)*0.5f-rotY*rot_vec_var(1)*t3*t11*t30*0.5f;
      P_[0][3] = rotX*rot_vec_var(0)*t3*t11*(t16-rotX*rotZ*t10*t12*0.5f)*0.5f+rotY*rot_vec_var(1)*t3*t11*(t17-rotY*rotZ*t10*t12*0.5f)*0.5f-rotZ*rot_vec_var(2)*t3*t11*t33*0.5f;
      P_[1][0] = t22;
      P_[1][1] = rot_vec_var(0)*(t19*t19)+rot_vec_var(1)*(t24*t24)+rot_vec_var(2)*(t26*t26);
      P_[1][2] = rot_vec_var(2)*(t16-t25)*(t17-rotY*rotZ*t10*t12*0.5f)-rot_vec_var(0)*t19*t28-rot_vec_var(1)*t28*t30;
      P_[1][3] = rot_vec_var(1)*(t15-t23)*(t17-rotY*rotZ*t10*t12*0.5f)-rot_vec_var(0)*t19*t31-rot_vec_var(2)*t31*t33;
      P_[2][0] = t35-rotY*rot_vec_var(1)*t3*t11*t30*0.5f+rotX*rot_vec_var(0)*t3*t11*(t15-t23)*0.5f;
      P_[2][1] = rot_vec_var(2)*(t16-t25)*(t17-t36)-rot_vec_var(0)*t19*t28-rot_vec_var(1)*t28*t30;
      P_[2][2] = rot_vec_var(1)*(t30*t30)+rot_vec_var(0)*(t37*t37)+rot_vec_var(2)*(t38*t38);
      P_[2][3] = t42;
      P_[3][0] = rotZ*rot_vec_var(2)*t3*t11*t33*(-0.5f)+rotX*rot_vec_var(0)*t3*t11*(t16-t25)*0.5f+rotY*rot_vec_var(1)*t3*t11*(t17-t36)*0.5f;
      P_[3][1] = rot_vec_var(1)*(t15-t23)*(t17-t36)-rot_vec_var(0)*t19*t31-rot_vec_var(2)*t31*t33;
      P_[3][2] = t42;
      P_[3][3] = rot_vec_var(2)*(t33*t33)+rot_vec_var(0)*(t43*t43)+rot_vec_var(1)*(t44*t44);
    } else {
      // the equations are badly conditioned so use a small angle approximation
      P_[0][0] = 0.0f;
      P_[0][1] = 0.0f;
      P_[0][2] = 0.0f;
      P_[0][3] = 0.0f;
      P_[1][0] = 0.0f;
      P_[1][1] = 0.25f * rot_vec_var(0);
      P_[1][2] = 0.0f;
      P_[1][3] = 0.0f;
      P_[2][0] = 0.0f;
      P_[2][1] = 0.0f;
      P_[2][2] = 0.25f * rot_vec_var(1);
      P_[2][3] = 0.0f;
      P_[3][0] = 0.0f;
      P_[3][1] = 0.0f;
      P_[3][2] = 0.0f;
      P_[3][3] = 0.25f * rot_vec_var(2);
    }
  }
  
  // zero specified range of rows in the state covariance matrix
  void ESKF::zeroRows(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last) {
    uint8_t row;
    for (row = first; row <= last; row++) {
      memset(&cov_mat[row][0], 0, sizeof(cov_mat[0][0]) * k_num_states_);
    }
  }

  // zero specified range of columns in the state covariance matrix
  void ESKF::zeroCols(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last) {
    uint8_t row;
    for (row = 0; row <= k_num_states_-1; row++) {
      memset(&cov_mat[row][first], 0, sizeof(cov_mat[0][0]) * (1 + last - first));
    }
  }
	
  void ESKF::constrainStates(scalar_t dt) {
	  state_.quat_nominal.w() = constrain(state_.quat_nominal.w(), -1.0, 1.0);
    state_.quat_nominal.x() = constrain(state_.quat_nominal.x(), -1.0, 1.0);
	  state_.quat_nominal.y() = constrain(state_.quat_nominal.y(), -1.0, 1.0);
	  state_.quat_nominal.z() = constrain(state_.quat_nominal.z(), -1.0, 1.0);
	  
    for (int i = 0; i < 3; i++) {
      state_.gyro_bias(i) = constrain(state_.gyro_bias(i), -0.349066 * dt, 0.349066 * dt);
    }
  }
} //  namespace eskf
