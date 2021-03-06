#ifndef NDEBUG
#define NDEBUG
#endif

#include "ESKF.hpp"
#include <cmath>

using namespace Eigen;

namespace eskf {

  template <typename T> inline T sq(T var) {
    return var * var;
  }
  
  template <typename T> inline T max(T val1, T val2) {
    return (val1 > val2) ? val1 : val2;
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

  quat ESKF::from_axis_angle(vec3 vec) {
    quat q;
    scalar_t theta = vec.norm();

    if (theta < scalar_t(1e-10)) {
      q.w() = scalar_t(1.0);
      q.x() = q.y() = q.z() = 0;
      return q;
    }

    vec3 tmp = vec / theta;
    return from_axis_angle(tmp, theta);
  }

  quat ESKF::from_axis_angle(const vec3 &axis, scalar_t theta) {
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

  vec3 ESKF::to_axis_angle(const quat& q) {
    scalar_t axis_magnitude = scalar_t(sqrt(q.x() * q.x() + q.y() * q.y() + q.z() * q.z()));
    vec3 vec;
    vec(0) = q.x();
    vec(1) = q.y();
    vec(2) = q.z();

    if (axis_magnitude >= scalar_t(1e-10)) {
      vec = vec / axis_magnitude;
      vec = vec * wrap_pi(scalar_t(2.0) * atan2(axis_magnitude, q.w()));
    }

    return vec;
  }

  mat3 ESKF::quat2dcm(const quat& q) {
    mat3 dcm;
    scalar_t a = q.w();
    scalar_t b = q.x();
    scalar_t c = q.y();
    scalar_t d = q.z();
    scalar_t aSq = a * a;
    scalar_t bSq = b * b;
    scalar_t cSq = c * c;
    scalar_t dSq = d * d;
    dcm(0, 0) = aSq + bSq - cSq - dSq;
    dcm(0, 1) = 2 * (b * c - a * d);
    dcm(0, 2) = 2 * (a * c + b * d);
    dcm(1, 0) = 2 * (b * c + a * d);
    dcm(1, 1) = aSq - bSq + cSq - dSq;
    dcm(1, 2) = 2 * (c * d - a * b);
    dcm(2, 0) = 2 * (b * d - a * c);
    dcm(2, 1) = 2 * (a * b + c * d);
    dcm(2, 2) = aSq - bSq - cSq + dSq;
    return dcm;
  }

  vec3 ESKF::dcm2vec(const mat3& dcm) {
    scalar_t phi_val = atan2(dcm(2, 1), dcm(2, 2));
    scalar_t theta_val = asin(-dcm(2, 0));
    scalar_t psi_val = atan2(dcm(1, 0), dcm(0, 0));
    scalar_t pi = M_PI;

    if (fabs(theta_val - pi / 2) < 1.0e-3) {
      phi_val = 0.0;
      psi_val = atan2(dcm(1, 2), dcm(0, 2));
    } else if (fabs(theta_val + pi / 2) < 1.0e-3) {
      phi_val = 0.0;
      psi_val = atan2(-dcm(1, 2), -dcm(0, 2));
    }
    return vec3(phi_val, theta_val, psi_val);
  }

  ESKF::ESKF() {
    // zeros state_
    state_.quat_nominal = quat(1, 0, 0, 0);
    state_.vel = vec3(0, 0, 0);
    state_.pos = vec3(0, 0, 0);
    state_.gyro_bias = vec3(0, 0, 0);
    state_.accel_bias = vec3(0, 0, 0);
    
    //  zeros P_
    for (unsigned i = 0; i < k_num_states_; i++) {
      for (unsigned j = 0; j < k_num_states_; j++) {
        P_[i][j] = 0.0f;
      }
    }

    q_ne = quat(0, 0.7071, 0.7071, 0); // rotation from ned to enu
    q_rb = quat(0, 1, 0, 0); // rotation ros body to px4 body

    _imu_down_sampled.delta_ang.setZero();
    _imu_down_sampled.delta_vel.setZero();
    _imu_down_sampled.delta_ang_dt = 0.0f;
    _imu_down_sampled.delta_vel_dt = 0.0f;

    _q_down_sampled.w() = 1.0f;
    _q_down_sampled.x() = 0.0f;
    _q_down_sampled.y() = 0.0f;
    _q_down_sampled.z() = 0.0f;

    _imu_buffer.allocate(_imu_buffer_length);
    for (int index = 0; index < _imu_buffer_length; index++) {
      imuSample imu_sample_init = {};
      _imu_buffer.push(imu_sample_init);
    }

    _ext_vision_buffer.allocate(_obs_buffer_length);
    for (int index = 0; index < _obs_buffer_length; index++) {
      extVisionSample ext_vision_sample_init = {};
      _ext_vision_buffer.push(ext_vision_sample_init);
    }
    
    _gps_buffer.allocate(_obs_buffer_length);
    for (int index = 0; index < _obs_buffer_length; index++) {
      gpsSample gps_sample_init = {};
      _gps_buffer.push(gps_sample_init);
    }
    
    last_known_posNED_ = vec3(0, 0, 0);
    _delVel_sum = vec3(0, 0, 0);

    _dt_ekf_avg = 0.001f * (scalar_t)(FILTER_UPDATE_PERIOD_MS);

    filter_initialised_ = false;
    imu_updated_ = false;
    memset(vel_pos_innov_,0,6*sizeof(float));
    time_last_hgt_fuse_ = time_last_imu_;
  }

  void ESKF::initialiseCovariance() {
    // define the initial angle uncertainty as variances for a rotation vector

    for (unsigned i = 0; i < k_num_states_; i++) {
      for (unsigned j = 0; j < k_num_states_; j++) {
	P_[i][j] = 0.0f;
      }
    }

    // calculate average prediction time step in sec
    float dt = 0.001f * (float)FILTER_UPDATE_PERIOD_MS;

    vec3 rot_vec_var;
    rot_vec_var(2) = rot_vec_var(1) = rot_vec_var(0) = sq(initial_tilt_err);

    // update the quaternion state covariances
    initialiseQuatCovariances(rot_vec_var);

    // velocity
    P_[4][4] = sq(fmaxf(vel_noise, 0.01f));
    P_[5][5] = P_[4][4];
    P_[6][6] = sq(1.5f) * P_[4][4];

    // position
    P_[7][7] = sq(fmaxf(pos_noise, 0.01f));
    P_[8][8] = P_[7][7];
    P_[9][9] = sq(fmaxf(baro_noise, 0.01f));

    // gyro bias
    P_[10][10] = sq(switch_on_gyro_bias * dt);
    P_[11][11] = P_[10][10];
    P_[12][12] = P_[10][10];
  }
  
  bool ESKF::initializeFilter() {
    scalar_t pitch = 0.0;
    scalar_t roll = 0.0;
    scalar_t yaw = 0.0;
    imuSample imu_init = _imu_buffer.get_newest();
    _delVel_sum += imu_init.delta_vel;
    if (_delVel_sum.norm() > 0.001) {
      _delVel_sum.normalize();
      pitch = asin(_delVel_sum(0));
      roll = atan2(-_delVel_sum(1), -_delVel_sum(2));
    } else {
      return false;
    }
    // calculate initial tilt alignment
    state_.quat_nominal = AngleAxis<scalar_t>(yaw, vec3::UnitZ()) * AngleAxis<scalar_t>(pitch, vec3::UnitY()) * AngleAxis<scalar_t>(roll, vec3::UnitX());
    initialiseCovariance();
    return true;    
  }
  
  bool ESKF::collect_imu(imuSample &imu) {
    // accumulate and downsample IMU data across a period FILTER_UPDATE_PERIOD_MS long

    // copy imu data to local variables
    _imu_sample_new.delta_ang	= imu.delta_ang;
    _imu_sample_new.delta_vel	= imu.delta_vel;
    _imu_sample_new.delta_ang_dt = imu.delta_ang_dt;
    _imu_sample_new.delta_vel_dt = imu.delta_vel_dt;
    _imu_sample_new.time_us	= imu.time_us;

    // accumulate the time deltas
    _imu_down_sampled.delta_ang_dt += imu.delta_ang_dt;
    _imu_down_sampled.delta_vel_dt += imu.delta_vel_dt;

    // use a quaternion to accumulate delta angle data
    // this quaternion represents the rotation from the start to end of the accumulation period
    quat delta_q(1, 0, 0, 0);
    quat res = from_axis_angle(imu.delta_ang);
    delta_q = delta_q * res;
    _q_down_sampled = _q_down_sampled * delta_q;
    _q_down_sampled.normalize();

    // rotate the accumulated delta velocity data forward each time so it is always in the updated rotation frame
    mat3 delta_R = quat2dcm(delta_q.inverse());
    _imu_down_sampled.delta_vel = delta_R * _imu_down_sampled.delta_vel;

    // accumulate the most recent delta velocity data at the updated rotation frame
    // assume effective sample time is halfway between the previous and current rotation frame
    _imu_down_sampled.delta_vel += (_imu_sample_new.delta_vel + delta_R * _imu_sample_new.delta_vel) * 0.5f;

    // if the target time delta between filter prediction steps has been exceeded
    // write the accumulated IMU data to the ring buffer
    float target_dt = (float)(FILTER_UPDATE_PERIOD_MS) / 1000;

    if (_imu_down_sampled.delta_ang_dt >= target_dt - _imu_collection_time_adj) {

      // accumulate the amount of time to advance the IMU collection time so that we meet the
      // average EKF update rate requirement
      _imu_collection_time_adj += 0.01f * (_imu_down_sampled.delta_ang_dt - target_dt);
      _imu_collection_time_adj = constrain(_imu_collection_time_adj, -0.5f * target_dt, 0.5f * target_dt);

      imu.delta_ang     = to_axis_angle(_q_down_sampled);
      imu.delta_vel     = _imu_down_sampled.delta_vel;
      imu.delta_ang_dt  = _imu_down_sampled.delta_ang_dt;
      imu.delta_vel_dt  = _imu_down_sampled.delta_vel_dt;

      _imu_down_sampled.delta_ang.setZero();
      _imu_down_sampled.delta_vel.setZero();
      _imu_down_sampled.delta_ang_dt = 0.0f;
      _imu_down_sampled.delta_vel_dt = 0.0f;
      _q_down_sampled.w() = 1.0f;
      _q_down_sampled.x() = _q_down_sampled.y() = _q_down_sampled.z() = 0.0f;

      return true;
    }

    min_obs_interval_us = (_imu_sample_new.time_us - _imu_sample_delayed.time_us) / (_obs_buffer_length - 1);

    return false;
  }
  
  void ESKF::predict(const vec3 &w, const vec3 &a, uint64_t time_us, scalar_t dt) {
    
    // convert ROS body to PX4 body frame IMU data
    vec3 px4body_w = q_rb.toRotationMatrix() * w;
    vec3 px4body_a = q_rb.toRotationMatrix() * a;
    
    vec3 delta_ang = px4body_w * dt; // current delta angle  (rad)
    vec3 delta_vel = px4body_a * dt; //current delta velocity (m/s)
    
    // copy data
    imuSample imu_sample_new = {};
    imu_sample_new.delta_ang = delta_ang;
    imu_sample_new.delta_vel = delta_vel;
    imu_sample_new.delta_ang_dt = dt;
    imu_sample_new.delta_vel_dt = dt;
    imu_sample_new.time_us = time_us;
    
    time_last_imu_ = time_us;
        
    if(collect_imu(imu_sample_new)) {
      _imu_buffer.push(imu_sample_new);
      imu_updated_ = true;
      // get the oldest data from the buffer
      _imu_sample_delayed = _imu_buffer.get_oldest();
    } else {
      imu_updated_ = false;
      return;
    }
    
    if (!filter_initialised_) {
      filter_initialised_ = initializeFilter();

      if (!filter_initialised_) {
        return;
      }
    }
    
    if(!imu_updated_) return;
    
    // apply imu bias corrections
    vec3 corrected_delta_ang = _imu_sample_delayed.delta_ang - state_.gyro_bias;
    vec3 corrected_delta_vel = _imu_sample_delayed.delta_vel - state_.accel_bias; 
    
    // convert the delta angle to a delta quaternion
    quat dq;
    dq = from_axis_angle(corrected_delta_ang);
    // rotate the previous quaternion by the delta quaternion using a quaternion multiplication
    state_.quat_nominal = state_.quat_nominal * dq;
    // quaternions must be normalised whenever they are modified
    state_.quat_nominal.normalize();
        
    // save the previous value of velocity so we can use trapezoidal integration
    vec3 vel_last = state_.vel;
    
    // update transformation matrix from body to world frame
    mat3 R_to_earth = quat_to_invrotmat(state_.quat_nominal);
    
    // Calculate an earth frame delta velocity
    vec3 corrected_delta_vel_ef = R_to_earth * corrected_delta_vel;
        
    // calculate the increment in velocity using the current orientation
    state_.vel += corrected_delta_vel_ef;

    // compensate for acceleration due to gravity
    state_.vel(2) += kOneG * _imu_sample_delayed.delta_vel_dt;
        
    // predict position states via trapezoidal integration of velocity
    state_.pos += (vel_last + state_.vel) * _imu_sample_delayed.delta_vel_dt * 0.5f;
        
    constrainStates();
        
    // calculate an average filter update time
    scalar_t input = 0.5f * (_imu_sample_delayed.delta_vel_dt + _imu_sample_delayed.delta_ang_dt);

    // filter and limit input between -50% and +100% of nominal value
    input = constrain(input, 0.0005f * (scalar_t)(FILTER_UPDATE_PERIOD_MS), 0.002f * (scalar_t)(FILTER_UPDATE_PERIOD_MS));
    _dt_ekf_avg = 0.99f * _dt_ekf_avg + 0.01f * input;
    
    predictCovariance();
    fuseVelPosHeight();
    controlHeightSensorTimeouts();
    fuseHeading();
  }
  
  void ESKF::predictCovariance() {
    // error-state jacobian
    // assign intermediate state variables
    scalar_t q0 = state_.quat_nominal.w();
    scalar_t q1 = state_.quat_nominal.x();
    scalar_t q2 = state_.quat_nominal.y();
    scalar_t q3 = state_.quat_nominal.z();

    scalar_t dax = _imu_sample_delayed.delta_ang(0);
    scalar_t day = _imu_sample_delayed.delta_ang(1);
    scalar_t daz = _imu_sample_delayed.delta_ang(2);

    scalar_t dvx = _imu_sample_delayed.delta_vel(0);
    scalar_t dvy = _imu_sample_delayed.delta_vel(1);
    scalar_t dvz = _imu_sample_delayed.delta_vel(2);

    scalar_t dax_b = state_.gyro_bias(0);
    scalar_t day_b = state_.gyro_bias(1);
    scalar_t daz_b = state_.gyro_bias(2);

    scalar_t dvx_b = state_.accel_bias(0);
    scalar_t dvy_b = state_.accel_bias(1);
    scalar_t dvz_b = state_.accel_bias(2);
	  
    // compute noise variance for stationary processes
    scalar_t process_noise[k_num_states_] = {};
    
    float dt = constrain(_imu_sample_delayed.delta_ang_dt, 0.0005f * (scalar_t)(FILTER_UPDATE_PERIOD_MS), 0.002f * (scalar_t)(FILTER_UPDATE_PERIOD_MS));
    
    // convert rate of change of rate gyro bias (rad/s**2) as specified by the parameter to an expected change in delta angle (rad) since the last update
    scalar_t d_ang_bias_sig = dt * dt * constrain(gyro_bias_p_noise, 0.0f, 1.0f);

    // convert rate of change of accelerometer bias (m/s**3) as specified by the parameter to an expected change in delta velocity (m/s) since the last update
    scalar_t d_vel_bias_sig = dt * dt * constrain(accel_bias_p_noise, 0.0f, 1.0f);
        
    // Construct the process noise variance diagonal for those states with a stationary process model
    // These are kinematic states and their error growth is controlled separately by the IMU noise variances
    for (unsigned i = 0; i <= 9; i++) {
      process_noise[i] = 0.0;
    }
    
    // delta angle bias states
    process_noise[12] = process_noise[11] = process_noise[10] = sq(d_ang_bias_sig);
    // delta_velocity bias states
    process_noise[15] = process_noise[14] = process_noise[13] = sq(d_vel_bias_sig);
        
    // assign IMU noise variances
    // inputs to the system are 3 delta angles and 3 delta velocities
    scalar_t daxVar, dayVar, dazVar;
    scalar_t dvxVar, dvyVar, dvzVar;
    gyro_noise = constrain(gyro_noise, 0.0f, 1.0f);
    daxVar = dayVar = dazVar = sq(dt * gyro_noise); // gyro prediction variance TODO get variance from sensor
    accel_noise = constrain(accel_noise, 0.0f, 1.0f);
    dvxVar = dvyVar = dvzVar = sq(dt * accel_noise); //accel prediction variance TODO get variance from sensor
    
    // intermediate calculations
    scalar_t SF[21];
    SF[0] = dvz - dvz_b;
    SF[1] = dvy - dvy_b;
    SF[2] = dvx - dvx_b;
    SF[3] = 2*q1*SF[2] + 2*q2*SF[1] + 2*q3*SF[0];
    SF[4] = 2*q0*SF[1] - 2*q1*SF[0] + 2*q3*SF[2];
    SF[5] = 2*q0*SF[2] + 2*q2*SF[0] - 2*q3*SF[1];
    SF[6] = day/2 - day_b/2;
    SF[7] = daz/2 - daz_b/2;
    SF[8] = dax/2 - dax_b/2;
    SF[9] = dax_b/2 - dax/2;
    SF[10] = daz_b/2 - daz/2;
    SF[11] = day_b/2 - day/2;
    SF[12] = 2*q1*SF[1];
    SF[13] = 2*q0*SF[0];
    SF[14] = q1/2;
    SF[15] = q2/2;
    SF[16] = q3/2;
    SF[17] = sq(q3);
    SF[18] = sq(q2);
    SF[19] = sq(q1);
    SF[20] = sq(q0);
    
    scalar_t SG[8];
    SG[0] = q0/2;
    SG[1] = sq(q3);
    SG[2] = sq(q2);
    SG[3] = sq(q1);
    SG[4] = sq(q0);
    SG[5] = 2*q2*q3;
    SG[6] = 2*q1*q3;
    SG[7] = 2*q1*q2;
    
    scalar_t SQ[11];
    SQ[0] = dvzVar*(SG[5] - 2*q0*q1)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvyVar*(SG[5] + 2*q0*q1)*(SG[1] - SG[2] + SG[3] - SG[4]) + dvxVar*(SG[6] - 2*q0*q2)*(SG[7] + 2*q0*q3);
    SQ[1] = dvzVar*(SG[6] + 2*q0*q2)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvxVar*(SG[6] - 2*q0*q2)*(SG[1] + SG[2] - SG[3] - SG[4]) + dvyVar*(SG[5] + 2*q0*q1)*(SG[7] - 2*q0*q3);
    SQ[2] = dvzVar*(SG[5] - 2*q0*q1)*(SG[6] + 2*q0*q2) - dvyVar*(SG[7] - 2*q0*q3)*(SG[1] - SG[2] + SG[3] - SG[4]) - dvxVar*(SG[7] + 2*q0*q3)*(SG[1] + SG[2] - SG[3] - SG[4]);
    SQ[3] = (dayVar*q1*SG[0])/2 - (dazVar*q1*SG[0])/2 - (daxVar*q2*q3)/4;
    SQ[4] = (dazVar*q2*SG[0])/2 - (daxVar*q2*SG[0])/2 - (dayVar*q1*q3)/4;
    SQ[5] = (daxVar*q3*SG[0])/2 - (dayVar*q3*SG[0])/2 - (dazVar*q1*q2)/4;
    SQ[6] = (daxVar*q1*q2)/4 - (dazVar*q3*SG[0])/2 - (dayVar*q1*q2)/4;
    SQ[7] = (dazVar*q1*q3)/4 - (daxVar*q1*q3)/4 - (dayVar*q2*SG[0])/2;
    SQ[8] = (dayVar*q2*q3)/4 - (daxVar*q1*SG[0])/2 - (dazVar*q2*q3)/4;
    SQ[9] = sq(SG[0]);
    SQ[10] = sq(q1);
    
    scalar_t SPP[11];
    SPP[0] = SF[12] + SF[13] - 2*q2*SF[2];
    SPP[1] = SF[17] - SF[18] - SF[19] + SF[20];
    SPP[2] = SF[17] - SF[18] + SF[19] - SF[20];
    SPP[3] = SF[17] + SF[18] - SF[19] - SF[20];
    SPP[4] = 2*q0*q2 - 2*q1*q3;
    SPP[5] = 2*q0*q1 - 2*q2*q3;
    SPP[6] = 2*q0*q3 - 2*q1*q2;
    SPP[7] = 2*q0*q1 + 2*q2*q3;
    SPP[8] = 2*q0*q3 + 2*q1*q2;
    SPP[9] = 2*q0*q2 + 2*q1*q3;
    SPP[10] = SF[16];
    
    // covariance update
    // calculate variances and upper diagonal covariances for quaternion, velocity, position and gyro bias states
    scalar_t nextP[k_num_states_][k_num_states_];
    nextP[0][0] = P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10] + (daxVar*SQ[10])/4 + SF[9]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) + SF[11]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SF[10]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) + SF[14]*(P_[0][10] + P_[1][10]*SF[9] + P_[2][10]*SF[11] + P_[3][10]*SF[10] + P_[10][10]*SF[14] + P_[11][10]*SF[15] + P_[12][10]*SPP[10]) + SF[15]*(P_[0][11] + P_[1][11]*SF[9] + P_[2][11]*SF[11] + P_[3][11]*SF[10] + P_[10][11]*SF[14] + P_[11][11]*SF[15] + P_[12][11]*SPP[10]) + SPP[10]*(P_[0][12] + P_[1][12]*SF[9] + P_[2][12]*SF[11] + P_[3][12]*SF[10] + P_[10][12]*SF[14] + P_[11][12]*SF[15] + P_[12][12]*SPP[10]) + (dayVar*sq(q2))/4 + (dazVar*sq(q3))/4;
    nextP[0][1] = P_[0][1] + SQ[8] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10] + SF[8]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SF[7]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SF[11]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) - SF[15]*(P_[0][12] + P_[1][12]*SF[9] + P_[2][12]*SF[11] + P_[3][12]*SF[10] + P_[10][12]*SF[14] + P_[11][12]*SF[15] + P_[12][12]*SPP[10]) + SPP[10]*(P_[0][11] + P_[1][11]*SF[9] + P_[2][11]*SF[11] + P_[3][11]*SF[10] + P_[10][11]*SF[14] + P_[11][11]*SF[15] + P_[12][11]*SPP[10]) - (q0*(P_[0][10] + P_[1][10]*SF[9] + P_[2][10]*SF[11] + P_[3][10]*SF[10] + P_[10][10]*SF[14] + P_[11][10]*SF[15] + P_[12][10]*SPP[10]))/2;
    nextP[1][1] = P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] + daxVar*SQ[9] - (P_[10][1]*q0)/2 + SF[8]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SF[7]*(P_[1][2] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2) + SF[11]*(P_[1][3] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2) - SF[15]*(P_[1][12] + P_[0][12]*SF[8] + P_[2][12]*SF[7] + P_[3][12]*SF[11] - P_[12][12]*SF[15] + P_[11][12]*SPP[10] - (P_[10][12]*q0)/2) + SPP[10]*(P_[1][11] + P_[0][11]*SF[8] + P_[2][11]*SF[7] + P_[3][11]*SF[11] - P_[12][11]*SF[15] + P_[11][11]*SPP[10] - (P_[10][11]*q0)/2) + (dayVar*sq(q3))/4 + (dazVar*sq(q2))/4 - (q0*(P_[1][10] + P_[0][10]*SF[8] + P_[2][10]*SF[7] + P_[3][10]*SF[11] - P_[12][10]*SF[15] + P_[11][10]*SPP[10] - (P_[10][10]*q0)/2))/2;
    nextP[0][2] = P_[0][2] + SQ[7] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10] + SF[6]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SF[10]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) + SF[8]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) + SF[14]*(P_[0][12] + P_[1][12]*SF[9] + P_[2][12]*SF[11] + P_[3][12]*SF[10] + P_[10][12]*SF[14] + P_[11][12]*SF[15] + P_[12][12]*SPP[10]) - SPP[10]*(P_[0][10] + P_[1][10]*SF[9] + P_[2][10]*SF[11] + P_[3][10]*SF[10] + P_[10][10]*SF[14] + P_[11][10]*SF[15] + P_[12][10]*SPP[10]) - (q0*(P_[0][11] + P_[1][11]*SF[9] + P_[2][11]*SF[11] + P_[3][11]*SF[10] + P_[10][11]*SF[14] + P_[11][11]*SF[15] + P_[12][11]*SPP[10]))/2;
    nextP[1][2] = P_[1][2] + SQ[5] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2 + SF[6]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SF[10]*(P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] - (P_[10][1]*q0)/2) + SF[8]*(P_[1][3] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2) + SF[14]*(P_[1][12] + P_[0][12]*SF[8] + P_[2][12]*SF[7] + P_[3][12]*SF[11] - P_[12][12]*SF[15] + P_[11][12]*SPP[10] - (P_[10][12]*q0)/2) - SPP[10]*(P_[1][10] + P_[0][10]*SF[8] + P_[2][10]*SF[7] + P_[3][10]*SF[11] - P_[12][10]*SF[15] + P_[11][10]*SPP[10] - (P_[10][10]*q0)/2) - (q0*(P_[1][11] + P_[0][11]*SF[8] + P_[2][11]*SF[7] + P_[3][11]*SF[11] - P_[12][11]*SF[15] + P_[11][11]*SPP[10] - (P_[10][11]*q0)/2))/2;
    nextP[2][2] = P_[2][2] + P_[0][2]*SF[6] + P_[1][2]*SF[10] + P_[3][2]*SF[8] + P_[12][2]*SF[14] - P_[10][2]*SPP[10] + dayVar*SQ[9] + (dazVar*SQ[10])/4 - (P_[11][2]*q0)/2 + SF[6]*(P_[2][0] + P_[0][0]*SF[6] + P_[1][0]*SF[10] + P_[3][0]*SF[8] + P_[12][0]*SF[14] - P_[10][0]*SPP[10] - (P_[11][0]*q0)/2) + SF[10]*(P_[2][1] + P_[0][1]*SF[6] + P_[1][1]*SF[10] + P_[3][1]*SF[8] + P_[12][1]*SF[14] - P_[10][1]*SPP[10] - (P_[11][1]*q0)/2) + SF[8]*(P_[2][3] + P_[0][3]*SF[6] + P_[1][3]*SF[10] + P_[3][3]*SF[8] + P_[12][3]*SF[14] - P_[10][3]*SPP[10] - (P_[11][3]*q0)/2) + SF[14]*(P_[2][12] + P_[0][12]*SF[6] + P_[1][12]*SF[10] + P_[3][12]*SF[8] + P_[12][12]*SF[14] - P_[10][12]*SPP[10] - (P_[11][12]*q0)/2) - SPP[10]*(P_[2][10] + P_[0][10]*SF[6] + P_[1][10]*SF[10] + P_[3][10]*SF[8] + P_[12][10]*SF[14] - P_[10][10]*SPP[10] - (P_[11][10]*q0)/2) + (daxVar*sq(q3))/4 - (q0*(P_[2][11] + P_[0][11]*SF[6] + P_[1][11]*SF[10] + P_[3][11]*SF[8] + P_[12][11]*SF[14] - P_[10][11]*SPP[10] - (P_[11][11]*q0)/2))/2;
    nextP[0][3] = P_[0][3] + SQ[6] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10] + SF[7]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SF[6]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) + SF[9]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SF[15]*(P_[0][10] + P_[1][10]*SF[9] + P_[2][10]*SF[11] + P_[3][10]*SF[10] + P_[10][10]*SF[14] + P_[11][10]*SF[15] + P_[12][10]*SPP[10]) - SF[14]*(P_[0][11] + P_[1][11]*SF[9] + P_[2][11]*SF[11] + P_[3][11]*SF[10] + P_[10][11]*SF[14] + P_[11][11]*SF[15] + P_[12][11]*SPP[10]) - (q0*(P_[0][12] + P_[1][12]*SF[9] + P_[2][12]*SF[11] + P_[3][12]*SF[10] + P_[10][12]*SF[14] + P_[11][12]*SF[15] + P_[12][12]*SPP[10]))/2;
    nextP[1][3] = P_[1][3] + SQ[4] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2 + SF[7]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SF[6]*(P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] - (P_[10][1]*q0)/2) + SF[9]*(P_[1][2] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2) + SF[15]*(P_[1][10] + P_[0][10]*SF[8] + P_[2][10]*SF[7] + P_[3][10]*SF[11] - P_[12][10]*SF[15] + P_[11][10]*SPP[10] - (P_[10][10]*q0)/2) - SF[14]*(P_[1][11] + P_[0][11]*SF[8] + P_[2][11]*SF[7] + P_[3][11]*SF[11] - P_[12][11]*SF[15] + P_[11][11]*SPP[10] - (P_[10][11]*q0)/2) - (q0*(P_[1][12] + P_[0][12]*SF[8] + P_[2][12]*SF[7] + P_[3][12]*SF[11] - P_[12][12]*SF[15] + P_[11][12]*SPP[10] - (P_[10][12]*q0)/2))/2;
    nextP[2][3] = P_[2][3] + SQ[3] + P_[0][3]*SF[6] + P_[1][3]*SF[10] + P_[3][3]*SF[8] + P_[12][3]*SF[14] - P_[10][3]*SPP[10] - (P_[11][3]*q0)/2 + SF[7]*(P_[2][0] + P_[0][0]*SF[6] + P_[1][0]*SF[10] + P_[3][0]*SF[8] + P_[12][0]*SF[14] - P_[10][0]*SPP[10] - (P_[11][0]*q0)/2) + SF[6]*(P_[2][1] + P_[0][1]*SF[6] + P_[1][1]*SF[10] + P_[3][1]*SF[8] + P_[12][1]*SF[14] - P_[10][1]*SPP[10] - (P_[11][1]*q0)/2) + SF[9]*(P_[2][2] + P_[0][2]*SF[6] + P_[1][2]*SF[10] + P_[3][2]*SF[8] + P_[12][2]*SF[14] - P_[10][2]*SPP[10] - (P_[11][2]*q0)/2) + SF[15]*(P_[2][10] + P_[0][10]*SF[6] + P_[1][10]*SF[10] + P_[3][10]*SF[8] + P_[12][10]*SF[14] - P_[10][10]*SPP[10] - (P_[11][10]*q0)/2) - SF[14]*(P_[2][11] + P_[0][11]*SF[6] + P_[1][11]*SF[10] + P_[3][11]*SF[8] + P_[12][11]*SF[14] - P_[10][11]*SPP[10] - (P_[11][11]*q0)/2) - (q0*(P_[2][12] + P_[0][12]*SF[6] + P_[1][12]*SF[10] + P_[3][12]*SF[8] + P_[12][12]*SF[14] - P_[10][12]*SPP[10] - (P_[11][12]*q0)/2))/2;
    nextP[3][3] = P_[3][3] + P_[0][3]*SF[7] + P_[1][3]*SF[6] + P_[2][3]*SF[9] + P_[10][3]*SF[15] - P_[11][3]*SF[14] + (dayVar*SQ[10])/4 + dazVar*SQ[9] - (P_[12][3]*q0)/2 + SF[7]*(P_[3][0] + P_[0][0]*SF[7] + P_[1][0]*SF[6] + P_[2][0]*SF[9] + P_[10][0]*SF[15] - P_[11][0]*SF[14] - (P_[12][0]*q0)/2) + SF[6]*(P_[3][1] + P_[0][1]*SF[7] + P_[1][1]*SF[6] + P_[2][1]*SF[9] + P_[10][1]*SF[15] - P_[11][1]*SF[14] - (P_[12][1]*q0)/2) + SF[9]*(P_[3][2] + P_[0][2]*SF[7] + P_[1][2]*SF[6] + P_[2][2]*SF[9] + P_[10][2]*SF[15] - P_[11][2]*SF[14] - (P_[12][2]*q0)/2) + SF[15]*(P_[3][10] + P_[0][10]*SF[7] + P_[1][10]*SF[6] + P_[2][10]*SF[9] + P_[10][10]*SF[15] - P_[11][10]*SF[14] - (P_[12][10]*q0)/2) - SF[14]*(P_[3][11] + P_[0][11]*SF[7] + P_[1][11]*SF[6] + P_[2][11]*SF[9] + P_[10][11]*SF[15] - P_[11][11]*SF[14] - (P_[12][11]*q0)/2) + (daxVar*sq(q2))/4 - (q0*(P_[3][12] + P_[0][12]*SF[7] + P_[1][12]*SF[6] + P_[2][12]*SF[9] + P_[10][12]*SF[15] - P_[11][12]*SF[14] - (P_[12][12]*q0)/2))/2;
    nextP[0][4] = P_[0][4] + P_[1][4]*SF[9] + P_[2][4]*SF[11] + P_[3][4]*SF[10] + P_[10][4]*SF[14] + P_[11][4]*SF[15] + P_[12][4]*SPP[10] + SF[5]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SF[3]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) - SF[4]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) + SPP[0]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SPP[3]*(P_[0][13] + P_[1][13]*SF[9] + P_[2][13]*SF[11] + P_[3][13]*SF[10] + P_[10][13]*SF[14] + P_[11][13]*SF[15] + P_[12][13]*SPP[10]) + SPP[6]*(P_[0][14] + P_[1][14]*SF[9] + P_[2][14]*SF[11] + P_[3][14]*SF[10] + P_[10][14]*SF[14] + P_[11][14]*SF[15] + P_[12][14]*SPP[10]) - SPP[9]*(P_[0][15] + P_[1][15]*SF[9] + P_[2][15]*SF[11] + P_[3][15]*SF[10] + P_[10][15]*SF[14] + P_[11][15]*SF[15] + P_[12][15]*SPP[10]);
    nextP[1][4] = P_[1][4] + P_[0][4]*SF[8] + P_[2][4]*SF[7] + P_[3][4]*SF[11] - P_[12][4]*SF[15] + P_[11][4]*SPP[10] - (P_[10][4]*q0)/2 + SF[5]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SF[3]*(P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] - (P_[10][1]*q0)/2) - SF[4]*(P_[1][3] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2) + SPP[0]*(P_[1][2] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2) + SPP[3]*(P_[1][13] + P_[0][13]*SF[8] + P_[2][13]*SF[7] + P_[3][13]*SF[11] - P_[12][13]*SF[15] + P_[11][13]*SPP[10] - (P_[10][13]*q0)/2) + SPP[6]*(P_[1][14] + P_[0][14]*SF[8] + P_[2][14]*SF[7] + P_[3][14]*SF[11] - P_[12][14]*SF[15] + P_[11][14]*SPP[10] - (P_[10][14]*q0)/2) - SPP[9]*(P_[1][15] + P_[0][15]*SF[8] + P_[2][15]*SF[7] + P_[3][15]*SF[11] - P_[12][15]*SF[15] + P_[11][15]*SPP[10] - (P_[10][15]*q0)/2);
    nextP[2][4] = P_[2][4] + P_[0][4]*SF[6] + P_[1][4]*SF[10] + P_[3][4]*SF[8] + P_[12][4]*SF[14] - P_[10][4]*SPP[10] - (P_[11][4]*q0)/2 + SF[5]*(P_[2][0] + P_[0][0]*SF[6] + P_[1][0]*SF[10] + P_[3][0]*SF[8] + P_[12][0]*SF[14] - P_[10][0]*SPP[10] - (P_[11][0]*q0)/2) + SF[3]*(P_[2][1] + P_[0][1]*SF[6] + P_[1][1]*SF[10] + P_[3][1]*SF[8] + P_[12][1]*SF[14] - P_[10][1]*SPP[10] - (P_[11][1]*q0)/2) - SF[4]*(P_[2][3] + P_[0][3]*SF[6] + P_[1][3]*SF[10] + P_[3][3]*SF[8] + P_[12][3]*SF[14] - P_[10][3]*SPP[10] - (P_[11][3]*q0)/2) + SPP[0]*(P_[2][2] + P_[0][2]*SF[6] + P_[1][2]*SF[10] + P_[3][2]*SF[8] + P_[12][2]*SF[14] - P_[10][2]*SPP[10] - (P_[11][2]*q0)/2) + SPP[3]*(P_[2][13] + P_[0][13]*SF[6] + P_[1][13]*SF[10] + P_[3][13]*SF[8] + P_[12][13]*SF[14] - P_[10][13]*SPP[10] - (P_[11][13]*q0)/2) + SPP[6]*(P_[2][14] + P_[0][14]*SF[6] + P_[1][14]*SF[10] + P_[3][14]*SF[8] + P_[12][14]*SF[14] - P_[10][14]*SPP[10] - (P_[11][14]*q0)/2) - SPP[9]*(P_[2][15] + P_[0][15]*SF[6] + P_[1][15]*SF[10] + P_[3][15]*SF[8] + P_[12][15]*SF[14] - P_[10][15]*SPP[10] - (P_[11][15]*q0)/2);
    nextP[3][4] = P_[3][4] + P_[0][4]*SF[7] + P_[1][4]*SF[6] + P_[2][4]*SF[9] + P_[10][4]*SF[15] - P_[11][4]*SF[14] - (P_[12][4]*q0)/2 + SF[5]*(P_[3][0] + P_[0][0]*SF[7] + P_[1][0]*SF[6] + P_[2][0]*SF[9] + P_[10][0]*SF[15] - P_[11][0]*SF[14] - (P_[12][0]*q0)/2) + SF[3]*(P_[3][1] + P_[0][1]*SF[7] + P_[1][1]*SF[6] + P_[2][1]*SF[9] + P_[10][1]*SF[15] - P_[11][1]*SF[14] - (P_[12][1]*q0)/2) - SF[4]*(P_[3][3] + P_[0][3]*SF[7] + P_[1][3]*SF[6] + P_[2][3]*SF[9] + P_[10][3]*SF[15] - P_[11][3]*SF[14] - (P_[12][3]*q0)/2) + SPP[0]*(P_[3][2] + P_[0][2]*SF[7] + P_[1][2]*SF[6] + P_[2][2]*SF[9] + P_[10][2]*SF[15] - P_[11][2]*SF[14] - (P_[12][2]*q0)/2) + SPP[3]*(P_[3][13] + P_[0][13]*SF[7] + P_[1][13]*SF[6] + P_[2][13]*SF[9] + P_[10][13]*SF[15] - P_[11][13]*SF[14] - (P_[12][13]*q0)/2) + SPP[6]*(P_[3][14] + P_[0][14]*SF[7] + P_[1][14]*SF[6] + P_[2][14]*SF[9] + P_[10][14]*SF[15] - P_[11][14]*SF[14] - (P_[12][14]*q0)/2) - SPP[9]*(P_[3][15] + P_[0][15]*SF[7] + P_[1][15]*SF[6] + P_[2][15]*SF[9] + P_[10][15]*SF[15] - P_[11][15]*SF[14] - (P_[12][15]*q0)/2);
    nextP[4][4] = P_[4][4] + P_[0][4]*SF[5] + P_[1][4]*SF[3] - P_[3][4]*SF[4] + P_[2][4]*SPP[0] + P_[13][4]*SPP[3] + P_[14][4]*SPP[6] - P_[15][4]*SPP[9] + dvyVar*sq(SG[7] - 2*q0*q3) + dvzVar*sq(SG[6] + 2*q0*q2) + SF[5]*(P_[4][0] + P_[0][0]*SF[5] + P_[1][0]*SF[3] - P_[3][0]*SF[4] + P_[2][0]*SPP[0] + P_[13][0]*SPP[3] + P_[14][0]*SPP[6] - P_[15][0]*SPP[9]) + SF[3]*(P_[4][1] + P_[0][1]*SF[5] + P_[1][1]*SF[3] - P_[3][1]*SF[4] + P_[2][1]*SPP[0] + P_[13][1]*SPP[3] + P_[14][1]*SPP[6] - P_[15][1]*SPP[9]) - SF[4]*(P_[4][3] + P_[0][3]*SF[5] + P_[1][3]*SF[3] - P_[3][3]*SF[4] + P_[2][3]*SPP[0] + P_[13][3]*SPP[3] + P_[14][3]*SPP[6] - P_[15][3]*SPP[9]) + SPP[0]*(P_[4][2] + P_[0][2]*SF[5] + P_[1][2]*SF[3] - P_[3][2]*SF[4] + P_[2][2]*SPP[0] + P_[13][2]*SPP[3] + P_[14][2]*SPP[6] - P_[15][2]*SPP[9]) + SPP[3]*(P_[4][13] + P_[0][13]*SF[5] + P_[1][13]*SF[3] - P_[3][13]*SF[4] + P_[2][13]*SPP[0] + P_[13][13]*SPP[3] + P_[14][13]*SPP[6] - P_[15][13]*SPP[9]) + SPP[6]*(P_[4][14] + P_[0][14]*SF[5] + P_[1][14]*SF[3] - P_[3][14]*SF[4] + P_[2][14]*SPP[0] + P_[13][14]*SPP[3] + P_[14][14]*SPP[6] - P_[15][14]*SPP[9]) - SPP[9]*(P_[4][15] + P_[0][15]*SF[5] + P_[1][15]*SF[3] - P_[3][15]*SF[4] + P_[2][15]*SPP[0] + P_[13][15]*SPP[3] + P_[14][15]*SPP[6] - P_[15][15]*SPP[9]) + dvxVar*sq(SG[1] + SG[2] - SG[3] - SG[4]);
    nextP[0][5] = P_[0][5] + P_[1][5]*SF[9] + P_[2][5]*SF[11] + P_[3][5]*SF[10] + P_[10][5]*SF[14] + P_[11][5]*SF[15] + P_[12][5]*SPP[10] + SF[4]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SF[3]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SF[5]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) - SPP[0]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) - SPP[8]*(P_[0][13] + P_[1][13]*SF[9] + P_[2][13]*SF[11] + P_[3][13]*SF[10] + P_[10][13]*SF[14] + P_[11][13]*SF[15] + P_[12][13]*SPP[10]) + SPP[2]*(P_[0][14] + P_[1][14]*SF[9] + P_[2][14]*SF[11] + P_[3][14]*SF[10] + P_[10][14]*SF[14] + P_[11][14]*SF[15] + P_[12][14]*SPP[10]) + SPP[5]*(P_[0][15] + P_[1][15]*SF[9] + P_[2][15]*SF[11] + P_[3][15]*SF[10] + P_[10][15]*SF[14] + P_[11][15]*SF[15] + P_[12][15]*SPP[10]);
    nextP[1][5] = P_[1][5] + P_[0][5]*SF[8] + P_[2][5]*SF[7] + P_[3][5]*SF[11] - P_[12][5]*SF[15] + P_[11][5]*SPP[10] - (P_[10][5]*q0)/2 + SF[4]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SF[3]*(P_[1][2] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2) + SF[5]*(P_[1][3] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2) - SPP[0]*(P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] - (P_[10][1]*q0)/2) - SPP[8]*(P_[1][13] + P_[0][13]*SF[8] + P_[2][13]*SF[7] + P_[3][13]*SF[11] - P_[12][13]*SF[15] + P_[11][13]*SPP[10] - (P_[10][13]*q0)/2) + SPP[2]*(P_[1][14] + P_[0][14]*SF[8] + P_[2][14]*SF[7] + P_[3][14]*SF[11] - P_[12][14]*SF[15] + P_[11][14]*SPP[10] - (P_[10][14]*q0)/2) + SPP[5]*(P_[1][15] + P_[0][15]*SF[8] + P_[2][15]*SF[7] + P_[3][15]*SF[11] - P_[12][15]*SF[15] + P_[11][15]*SPP[10] - (P_[10][15]*q0)/2);
    nextP[2][5] = P_[2][5] + P_[0][5]*SF[6] + P_[1][5]*SF[10] + P_[3][5]*SF[8] + P_[12][5]*SF[14] - P_[10][5]*SPP[10] - (P_[11][5]*q0)/2 + SF[4]*(P_[2][0] + P_[0][0]*SF[6] + P_[1][0]*SF[10] + P_[3][0]*SF[8] + P_[12][0]*SF[14] - P_[10][0]*SPP[10] - (P_[11][0]*q0)/2) + SF[3]*(P_[2][2] + P_[0][2]*SF[6] + P_[1][2]*SF[10] + P_[3][2]*SF[8] + P_[12][2]*SF[14] - P_[10][2]*SPP[10] - (P_[11][2]*q0)/2) + SF[5]*(P_[2][3] + P_[0][3]*SF[6] + P_[1][3]*SF[10] + P_[3][3]*SF[8] + P_[12][3]*SF[14] - P_[10][3]*SPP[10] - (P_[11][3]*q0)/2) - SPP[0]*(P_[2][1] + P_[0][1]*SF[6] + P_[1][1]*SF[10] + P_[3][1]*SF[8] + P_[12][1]*SF[14] - P_[10][1]*SPP[10] - (P_[11][1]*q0)/2) - SPP[8]*(P_[2][13] + P_[0][13]*SF[6] + P_[1][13]*SF[10] + P_[3][13]*SF[8] + P_[12][13]*SF[14] - P_[10][13]*SPP[10] - (P_[11][13]*q0)/2) + SPP[2]*(P_[2][14] + P_[0][14]*SF[6] + P_[1][14]*SF[10] + P_[3][14]*SF[8] + P_[12][14]*SF[14] - P_[10][14]*SPP[10] - (P_[11][14]*q0)/2) + SPP[5]*(P_[2][15] + P_[0][15]*SF[6] + P_[1][15]*SF[10] + P_[3][15]*SF[8] + P_[12][15]*SF[14] - P_[10][15]*SPP[10] - (P_[11][15]*q0)/2);
    nextP[3][5] = P_[3][5] + P_[0][5]*SF[7] + P_[1][5]*SF[6] + P_[2][5]*SF[9] + P_[10][5]*SF[15] - P_[11][5]*SF[14] - (P_[12][5]*q0)/2 + SF[4]*(P_[3][0] + P_[0][0]*SF[7] + P_[1][0]*SF[6] + P_[2][0]*SF[9] + P_[10][0]*SF[15] - P_[11][0]*SF[14] - (P_[12][0]*q0)/2) + SF[3]*(P_[3][2] + P_[0][2]*SF[7] + P_[1][2]*SF[6] + P_[2][2]*SF[9] + P_[10][2]*SF[15] - P_[11][2]*SF[14] - (P_[12][2]*q0)/2) + SF[5]*(P_[3][3] + P_[0][3]*SF[7] + P_[1][3]*SF[6] + P_[2][3]*SF[9] + P_[10][3]*SF[15] - P_[11][3]*SF[14] - (P_[12][3]*q0)/2) - SPP[0]*(P_[3][1] + P_[0][1]*SF[7] + P_[1][1]*SF[6] + P_[2][1]*SF[9] + P_[10][1]*SF[15] - P_[11][1]*SF[14] - (P_[12][1]*q0)/2) - SPP[8]*(P_[3][13] + P_[0][13]*SF[7] + P_[1][13]*SF[6] + P_[2][13]*SF[9] + P_[10][13]*SF[15] - P_[11][13]*SF[14] - (P_[12][13]*q0)/2) + SPP[2]*(P_[3][14] + P_[0][14]*SF[7] + P_[1][14]*SF[6] + P_[2][14]*SF[9] + P_[10][14]*SF[15] - P_[11][14]*SF[14] - (P_[12][14]*q0)/2) + SPP[5]*(P_[3][15] + P_[0][15]*SF[7] + P_[1][15]*SF[6] + P_[2][15]*SF[9] + P_[10][15]*SF[15] - P_[11][15]*SF[14] - (P_[12][15]*q0)/2);
    nextP[4][5] = P_[4][5] + SQ[2] + P_[0][5]*SF[5] + P_[1][5]*SF[3] - P_[3][5]*SF[4] + P_[2][5]*SPP[0] + P_[13][5]*SPP[3] + P_[14][5]*SPP[6] - P_[15][5]*SPP[9] + SF[4]*(P_[4][0] + P_[0][0]*SF[5] + P_[1][0]*SF[3] - P_[3][0]*SF[4] + P_[2][0]*SPP[0] + P_[13][0]*SPP[3] + P_[14][0]*SPP[6] - P_[15][0]*SPP[9]) + SF[3]*(P_[4][2] + P_[0][2]*SF[5] + P_[1][2]*SF[3] - P_[3][2]*SF[4] + P_[2][2]*SPP[0] + P_[13][2]*SPP[3] + P_[14][2]*SPP[6] - P_[15][2]*SPP[9]) + SF[5]*(P_[4][3] + P_[0][3]*SF[5] + P_[1][3]*SF[3] - P_[3][3]*SF[4] + P_[2][3]*SPP[0] + P_[13][3]*SPP[3] + P_[14][3]*SPP[6] - P_[15][3]*SPP[9]) - SPP[0]*(P_[4][1] + P_[0][1]*SF[5] + P_[1][1]*SF[3] - P_[3][1]*SF[4] + P_[2][1]*SPP[0] + P_[13][1]*SPP[3] + P_[14][1]*SPP[6] - P_[15][1]*SPP[9]) - SPP[8]*(P_[4][13] + P_[0][13]*SF[5] + P_[1][13]*SF[3] - P_[3][13]*SF[4] + P_[2][13]*SPP[0] + P_[13][13]*SPP[3] + P_[14][13]*SPP[6] - P_[15][13]*SPP[9]) + SPP[2]*(P_[4][14] + P_[0][14]*SF[5] + P_[1][14]*SF[3] - P_[3][14]*SF[4] + P_[2][14]*SPP[0] + P_[13][14]*SPP[3] + P_[14][14]*SPP[6] - P_[15][14]*SPP[9]) + SPP[5]*(P_[4][15] + P_[0][15]*SF[5] + P_[1][15]*SF[3] - P_[3][15]*SF[4] + P_[2][15]*SPP[0] + P_[13][15]*SPP[3] + P_[14][15]*SPP[6] - P_[15][15]*SPP[9]);
    nextP[5][5] = P_[5][5] + P_[0][5]*SF[4] + P_[2][5]*SF[3] + P_[3][5]*SF[5] - P_[1][5]*SPP[0] - P_[13][5]*SPP[8] + P_[14][5]*SPP[2] + P_[15][5]*SPP[5] + dvxVar*sq(SG[7] + 2*q0*q3) + dvzVar*sq(SG[5] - 2*q0*q1) + SF[4]*(P_[5][0] + P_[0][0]*SF[4] + P_[2][0]*SF[3] + P_[3][0]*SF[5] - P_[1][0]*SPP[0] - P_[13][0]*SPP[8] + P_[14][0]*SPP[2] + P_[15][0]*SPP[5]) + SF[3]*(P_[5][2] + P_[0][2]*SF[4] + P_[2][2]*SF[3] + P_[3][2]*SF[5] - P_[1][2]*SPP[0] - P_[13][2]*SPP[8] + P_[14][2]*SPP[2] + P_[15][2]*SPP[5]) + SF[5]*(P_[5][3] + P_[0][3]*SF[4] + P_[2][3]*SF[3] + P_[3][3]*SF[5] - P_[1][3]*SPP[0] - P_[13][3]*SPP[8] + P_[14][3]*SPP[2] + P_[15][3]*SPP[5]) - SPP[0]*(P_[5][1] + P_[0][1]*SF[4] + P_[2][1]*SF[3] + P_[3][1]*SF[5] - P_[1][1]*SPP[0] - P_[13][1]*SPP[8] + P_[14][1]*SPP[2] + P_[15][1]*SPP[5]) - SPP[8]*(P_[5][13] + P_[0][13]*SF[4] + P_[2][13]*SF[3] + P_[3][13]*SF[5] - P_[1][13]*SPP[0] - P_[13][13]*SPP[8] + P_[14][13]*SPP[2] + P_[15][13]*SPP[5]) + SPP[2]*(P_[5][14] + P_[0][14]*SF[4] + P_[2][14]*SF[3] + P_[3][14]*SF[5] - P_[1][14]*SPP[0] - P_[13][14]*SPP[8] + P_[14][14]*SPP[2] + P_[15][14]*SPP[5]) + SPP[5]*(P_[5][15] + P_[0][15]*SF[4] + P_[2][15]*SF[3] + P_[3][15]*SF[5] - P_[1][15]*SPP[0] - P_[13][15]*SPP[8] + P_[14][15]*SPP[2] + P_[15][15]*SPP[5]) + dvyVar*sq(SG[1] - SG[2] + SG[3] - SG[4]);
    nextP[0][6] = P_[0][6] + P_[1][6]*SF[9] + P_[2][6]*SF[11] + P_[3][6]*SF[10] + P_[10][6]*SF[14] + P_[11][6]*SF[15] + P_[12][6]*SPP[10] + SF[4]*(P_[0][1] + P_[1][1]*SF[9] + P_[2][1]*SF[11] + P_[3][1]*SF[10] + P_[10][1]*SF[14] + P_[11][1]*SF[15] + P_[12][1]*SPP[10]) - SF[5]*(P_[0][2] + P_[1][2]*SF[9] + P_[2][2]*SF[11] + P_[3][2]*SF[10] + P_[10][2]*SF[14] + P_[11][2]*SF[15] + P_[12][2]*SPP[10]) + SF[3]*(P_[0][3] + P_[1][3]*SF[9] + P_[2][3]*SF[11] + P_[3][3]*SF[10] + P_[10][3]*SF[14] + P_[11][3]*SF[15] + P_[12][3]*SPP[10]) + SPP[0]*(P_[0][0] + P_[1][0]*SF[9] + P_[2][0]*SF[11] + P_[3][0]*SF[10] + P_[10][0]*SF[14] + P_[11][0]*SF[15] + P_[12][0]*SPP[10]) + SPP[4]*(P_[0][13] + P_[1][13]*SF[9] + P_[2][13]*SF[11] + P_[3][13]*SF[10] + P_[10][13]*SF[14] + P_[11][13]*SF[15] + P_[12][13]*SPP[10]) - SPP[7]*(P_[0][14] + P_[1][14]*SF[9] + P_[2][14]*SF[11] + P_[3][14]*SF[10] + P_[10][14]*SF[14] + P_[11][14]*SF[15] + P_[12][14]*SPP[10]) - SPP[1]*(P_[0][15] + P_[1][15]*SF[9] + P_[2][15]*SF[11] + P_[3][15]*SF[10] + P_[10][15]*SF[14] + P_[11][15]*SF[15] + P_[12][15]*SPP[10]);
    nextP[1][6] = P_[1][6] + P_[0][6]*SF[8] + P_[2][6]*SF[7] + P_[3][6]*SF[11] - P_[12][6]*SF[15] + P_[11][6]*SPP[10] - (P_[10][6]*q0)/2 + SF[4]*(P_[1][1] + P_[0][1]*SF[8] + P_[2][1]*SF[7] + P_[3][1]*SF[11] - P_[12][1]*SF[15] + P_[11][1]*SPP[10] - (P_[10][1]*q0)/2) - SF[5]*(P_[1][2] + P_[0][2]*SF[8] + P_[2][2]*SF[7] + P_[3][2]*SF[11] - P_[12][2]*SF[15] + P_[11][2]*SPP[10] - (P_[10][2]*q0)/2) + SF[3]*(P_[1][3] + P_[0][3]*SF[8] + P_[2][3]*SF[7] + P_[3][3]*SF[11] - P_[12][3]*SF[15] + P_[11][3]*SPP[10] - (P_[10][3]*q0)/2) + SPP[0]*(P_[1][0] + P_[0][0]*SF[8] + P_[2][0]*SF[7] + P_[3][0]*SF[11] - P_[12][0]*SF[15] + P_[11][0]*SPP[10] - (P_[10][0]*q0)/2) + SPP[4]*(P_[1][13] + P_[0][13]*SF[8] + P_[2][13]*SF[7] + P_[3][13]*SF[11] - P_[12][13]*SF[15] + P_[11][13]*SPP[10] - (P_[10][13]*q0)/2) - SPP[7]*(P_[1][14] + P_[0][14]*SF[8] + P_[2][14]*SF[7] + P_[3][14]*SF[11] - P_[12][14]*SF[15] + P_[11][14]*SPP[10] - (P_[10][14]*q0)/2) - SPP[1]*(P_[1][15] + P_[0][15]*SF[8] + P_[2][15]*SF[7] + P_[3][15]*SF[11] - P_[12][15]*SF[15] + P_[11][15]*SPP[10] - (P_[10][15]*q0)/2);
    nextP[2][6] = P_[2][6] + P_[0][6]*SF[6] + P_[1][6]*SF[10] + P_[3][6]*SF[8] + P_[12][6]*SF[14] - P_[10][6]*SPP[10] - (P_[11][6]*q0)/2 + SF[4]*(P_[2][1] + P_[0][1]*SF[6] + P_[1][1]*SF[10] + P_[3][1]*SF[8] + P_[12][1]*SF[14] - P_[10][1]*SPP[10] - (P_[11][1]*q0)/2) - SF[5]*(P_[2][2] + P_[0][2]*SF[6] + P_[1][2]*SF[10] + P_[3][2]*SF[8] + P_[12][2]*SF[14] - P_[10][2]*SPP[10] - (P_[11][2]*q0)/2) + SF[3]*(P_[2][3] + P_[0][3]*SF[6] + P_[1][3]*SF[10] + P_[3][3]*SF[8] + P_[12][3]*SF[14] - P_[10][3]*SPP[10] - (P_[11][3]*q0)/2) + SPP[0]*(P_[2][0] + P_[0][0]*SF[6] + P_[1][0]*SF[10] + P_[3][0]*SF[8] + P_[12][0]*SF[14] - P_[10][0]*SPP[10] - (P_[11][0]*q0)/2) + SPP[4]*(P_[2][13] + P_[0][13]*SF[6] + P_[1][13]*SF[10] + P_[3][13]*SF[8] + P_[12][13]*SF[14] - P_[10][13]*SPP[10] - (P_[11][13]*q0)/2) - SPP[7]*(P_[2][14] + P_[0][14]*SF[6] + P_[1][14]*SF[10] + P_[3][14]*SF[8] + P_[12][14]*SF[14] - P_[10][14]*SPP[10] - (P_[11][14]*q0)/2) - SPP[1]*(P_[2][15] + P_[0][15]*SF[6] + P_[1][15]*SF[10] + P_[3][15]*SF[8] + P_[12][15]*SF[14] - P_[10][15]*SPP[10] - (P_[11][15]*q0)/2);
    nextP[3][6] = P_[3][6] + P_[0][6]*SF[7] + P_[1][6]*SF[6] + P_[2][6]*SF[9] + P_[10][6]*SF[15] - P_[11][6]*SF[14] - (P_[12][6]*q0)/2 + SF[4]*(P_[3][1] + P_[0][1]*SF[7] + P_[1][1]*SF[6] + P_[2][1]*SF[9] + P_[10][1]*SF[15] - P_[11][1]*SF[14] - (P_[12][1]*q0)/2) - SF[5]*(P_[3][2] + P_[0][2]*SF[7] + P_[1][2]*SF[6] + P_[2][2]*SF[9] + P_[10][2]*SF[15] - P_[11][2]*SF[14] - (P_[12][2]*q0)/2) + SF[3]*(P_[3][3] + P_[0][3]*SF[7] + P_[1][3]*SF[6] + P_[2][3]*SF[9] + P_[10][3]*SF[15] - P_[11][3]*SF[14] - (P_[12][3]*q0)/2) + SPP[0]*(P_[3][0] + P_[0][0]*SF[7] + P_[1][0]*SF[6] + P_[2][0]*SF[9] + P_[10][0]*SF[15] - P_[11][0]*SF[14] - (P_[12][0]*q0)/2) + SPP[4]*(P_[3][13] + P_[0][13]*SF[7] + P_[1][13]*SF[6] + P_[2][13]*SF[9] + P_[10][13]*SF[15] - P_[11][13]*SF[14] - (P_[12][13]*q0)/2) - SPP[7]*(P_[3][14] + P_[0][14]*SF[7] + P_[1][14]*SF[6] + P_[2][14]*SF[9] + P_[10][14]*SF[15] - P_[11][14]*SF[14] - (P_[12][14]*q0)/2) - SPP[1]*(P_[3][15] + P_[0][15]*SF[7] + P_[1][15]*SF[6] + P_[2][15]*SF[9] + P_[10][15]*SF[15] - P_[11][15]*SF[14] - (P_[12][15]*q0)/2);
    nextP[4][6] = P_[4][6] + SQ[1] + P_[0][6]*SF[5] + P_[1][6]*SF[3] - P_[3][6]*SF[4] + P_[2][6]*SPP[0] + P_[13][6]*SPP[3] + P_[14][6]*SPP[6] - P_[15][6]*SPP[9] + SF[4]*(P_[4][1] + P_[0][1]*SF[5] + P_[1][1]*SF[3] - P_[3][1]*SF[4] + P_[2][1]*SPP[0] + P_[13][1]*SPP[3] + P_[14][1]*SPP[6] - P_[15][1]*SPP[9]) - SF[5]*(P_[4][2] + P_[0][2]*SF[5] + P_[1][2]*SF[3] - P_[3][2]*SF[4] + P_[2][2]*SPP[0] + P_[13][2]*SPP[3] + P_[14][2]*SPP[6] - P_[15][2]*SPP[9]) + SF[3]*(P_[4][3] + P_[0][3]*SF[5] + P_[1][3]*SF[3] - P_[3][3]*SF[4] + P_[2][3]*SPP[0] + P_[13][3]*SPP[3] + P_[14][3]*SPP[6] - P_[15][3]*SPP[9]) + SPP[0]*(P_[4][0] + P_[0][0]*SF[5] + P_[1][0]*SF[3] - P_[3][0]*SF[4] + P_[2][0]*SPP[0] + P_[13][0]*SPP[3] + P_[14][0]*SPP[6] - P_[15][0]*SPP[9]) + SPP[4]*(P_[4][13] + P_[0][13]*SF[5] + P_[1][13]*SF[3] - P_[3][13]*SF[4] + P_[2][13]*SPP[0] + P_[13][13]*SPP[3] + P_[14][13]*SPP[6] - P_[15][13]*SPP[9]) - SPP[7]*(P_[4][14] + P_[0][14]*SF[5] + P_[1][14]*SF[3] - P_[3][14]*SF[4] + P_[2][14]*SPP[0] + P_[13][14]*SPP[3] + P_[14][14]*SPP[6] - P_[15][14]*SPP[9]) - SPP[1]*(P_[4][15] + P_[0][15]*SF[5] + P_[1][15]*SF[3] - P_[3][15]*SF[4] + P_[2][15]*SPP[0] + P_[13][15]*SPP[3] + P_[14][15]*SPP[6] - P_[15][15]*SPP[9]);
    nextP[5][6] = P_[5][6] + SQ[0] + P_[0][6]*SF[4] + P_[2][6]*SF[3] + P_[3][6]*SF[5] - P_[1][6]*SPP[0] - P_[13][6]*SPP[8] + P_[14][6]*SPP[2] + P_[15][6]*SPP[5] + SF[4]*(P_[5][1] + P_[0][1]*SF[4] + P_[2][1]*SF[3] + P_[3][1]*SF[5] - P_[1][1]*SPP[0] - P_[13][1]*SPP[8] + P_[14][1]*SPP[2] + P_[15][1]*SPP[5]) - SF[5]*(P_[5][2] + P_[0][2]*SF[4] + P_[2][2]*SF[3] + P_[3][2]*SF[5] - P_[1][2]*SPP[0] - P_[13][2]*SPP[8] + P_[14][2]*SPP[2] + P_[15][2]*SPP[5]) + SF[3]*(P_[5][3] + P_[0][3]*SF[4] + P_[2][3]*SF[3] + P_[3][3]*SF[5] - P_[1][3]*SPP[0] - P_[13][3]*SPP[8] + P_[14][3]*SPP[2] + P_[15][3]*SPP[5]) + SPP[0]*(P_[5][0] + P_[0][0]*SF[4] + P_[2][0]*SF[3] + P_[3][0]*SF[5] - P_[1][0]*SPP[0] - P_[13][0]*SPP[8] + P_[14][0]*SPP[2] + P_[15][0]*SPP[5]) + SPP[4]*(P_[5][13] + P_[0][13]*SF[4] + P_[2][13]*SF[3] + P_[3][13]*SF[5] - P_[1][13]*SPP[0] - P_[13][13]*SPP[8] + P_[14][13]*SPP[2] + P_[15][13]*SPP[5]) - SPP[7]*(P_[5][14] + P_[0][14]*SF[4] + P_[2][14]*SF[3] + P_[3][14]*SF[5] - P_[1][14]*SPP[0] - P_[13][14]*SPP[8] + P_[14][14]*SPP[2] + P_[15][14]*SPP[5]) - SPP[1]*(P_[5][15] + P_[0][15]*SF[4] + P_[2][15]*SF[3] + P_[3][15]*SF[5] - P_[1][15]*SPP[0] - P_[13][15]*SPP[8] + P_[14][15]*SPP[2] + P_[15][15]*SPP[5]);
    nextP[6][6] = P_[6][6] + P_[1][6]*SF[4] - P_[2][6]*SF[5] + P_[3][6]*SF[3] + P_[0][6]*SPP[0] + P_[13][6]*SPP[4] - P_[14][6]*SPP[7] - P_[15][6]*SPP[1] + dvxVar*sq(SG[6] - 2*q0*q2) + dvyVar*sq(SG[5] + 2*q0*q1) + SF[4]*(P_[6][1] + P_[1][1]*SF[4] - P_[2][1]*SF[5] + P_[3][1]*SF[3] + P_[0][1]*SPP[0] + P_[13][1]*SPP[4] - P_[14][1]*SPP[7] - P_[15][1]*SPP[1]) - SF[5]*(P_[6][2] + P_[1][2]*SF[4] - P_[2][2]*SF[5] + P_[3][2]*SF[3] + P_[0][2]*SPP[0] + P_[13][2]*SPP[4] - P_[14][2]*SPP[7] - P_[15][2]*SPP[1]) + SF[3]*(P_[6][3] + P_[1][3]*SF[4] - P_[2][3]*SF[5] + P_[3][3]*SF[3] + P_[0][3]*SPP[0] + P_[13][3]*SPP[4] - P_[14][3]*SPP[7] - P_[15][3]*SPP[1]) + SPP[0]*(P_[6][0] + P_[1][0]*SF[4] - P_[2][0]*SF[5] + P_[3][0]*SF[3] + P_[0][0]*SPP[0] + P_[13][0]*SPP[4] - P_[14][0]*SPP[7] - P_[15][0]*SPP[1]) + SPP[4]*(P_[6][13] + P_[1][13]*SF[4] - P_[2][13]*SF[5] + P_[3][13]*SF[3] + P_[0][13]*SPP[0] + P_[13][13]*SPP[4] - P_[14][13]*SPP[7] - P_[15][13]*SPP[1]) - SPP[7]*(P_[6][14] + P_[1][14]*SF[4] - P_[2][14]*SF[5] + P_[3][14]*SF[3] + P_[0][14]*SPP[0] + P_[13][14]*SPP[4] - P_[14][14]*SPP[7] - P_[15][14]*SPP[1]) - SPP[1]*(P_[6][15] + P_[1][15]*SF[4] - P_[2][15]*SF[5] + P_[3][15]*SF[3] + P_[0][15]*SPP[0] + P_[13][15]*SPP[4] - P_[14][15]*SPP[7] - P_[15][15]*SPP[1]) + dvzVar*sq(SG[1] - SG[2] - SG[3] + SG[4]);
    nextP[0][7] = P_[0][7] + P_[1][7]*SF[9] + P_[2][7]*SF[11] + P_[3][7]*SF[10] + P_[10][7]*SF[14] + P_[11][7]*SF[15] + P_[12][7]*SPP[10] + dt*(P_[0][4] + P_[1][4]*SF[9] + P_[2][4]*SF[11] + P_[3][4]*SF[10] + P_[10][4]*SF[14] + P_[11][4]*SF[15] + P_[12][4]*SPP[10]);
    nextP[1][7] = P_[1][7] + P_[0][7]*SF[8] + P_[2][7]*SF[7] + P_[3][7]*SF[11] - P_[12][7]*SF[15] + P_[11][7]*SPP[10] - (P_[10][7]*q0)/2 + dt*(P_[1][4] + P_[0][4]*SF[8] + P_[2][4]*SF[7] + P_[3][4]*SF[11] - P_[12][4]*SF[15] + P_[11][4]*SPP[10] - (P_[10][4]*q0)/2);
    nextP[2][7] = P_[2][7] + P_[0][7]*SF[6] + P_[1][7]*SF[10] + P_[3][7]*SF[8] + P_[12][7]*SF[14] - P_[10][7]*SPP[10] - (P_[11][7]*q0)/2 + dt*(P_[2][4] + P_[0][4]*SF[6] + P_[1][4]*SF[10] + P_[3][4]*SF[8] + P_[12][4]*SF[14] - P_[10][4]*SPP[10] - (P_[11][4]*q0)/2);
    nextP[3][7] = P_[3][7] + P_[0][7]*SF[7] + P_[1][7]*SF[6] + P_[2][7]*SF[9] + P_[10][7]*SF[15] - P_[11][7]*SF[14] - (P_[12][7]*q0)/2 + dt*(P_[3][4] + P_[0][4]*SF[7] + P_[1][4]*SF[6] + P_[2][4]*SF[9] + P_[10][4]*SF[15] - P_[11][4]*SF[14] - (P_[12][4]*q0)/2);
    nextP[4][7] = P_[4][7] + P_[0][7]*SF[5] + P_[1][7]*SF[3] - P_[3][7]*SF[4] + P_[2][7]*SPP[0] + P_[13][7]*SPP[3] + P_[14][7]*SPP[6] - P_[15][7]*SPP[9] + dt*(P_[4][4] + P_[0][4]*SF[5] + P_[1][4]*SF[3] - P_[3][4]*SF[4] + P_[2][4]*SPP[0] + P_[13][4]*SPP[3] + P_[14][4]*SPP[6] - P_[15][4]*SPP[9]);
    nextP[5][7] = P_[5][7] + P_[0][7]*SF[4] + P_[2][7]*SF[3] + P_[3][7]*SF[5] - P_[1][7]*SPP[0] - P_[13][7]*SPP[8] + P_[14][7]*SPP[2] + P_[15][7]*SPP[5] + dt*(P_[5][4] + P_[0][4]*SF[4] + P_[2][4]*SF[3] + P_[3][4]*SF[5] - P_[1][4]*SPP[0] - P_[13][4]*SPP[8] + P_[14][4]*SPP[2] + P_[15][4]*SPP[5]);
    nextP[6][7] = P_[6][7] + P_[1][7]*SF[4] - P_[2][7]*SF[5] + P_[3][7]*SF[3] + P_[0][7]*SPP[0] + P_[13][7]*SPP[4] - P_[14][7]*SPP[7] - P_[15][7]*SPP[1] + dt*(P_[6][4] + P_[1][4]*SF[4] - P_[2][4]*SF[5] + P_[3][4]*SF[3] + P_[0][4]*SPP[0] + P_[13][4]*SPP[4] - P_[14][4]*SPP[7] - P_[15][4]*SPP[1]);
    nextP[7][7] = P_[7][7] + P_[4][7]*dt + dt*(P_[7][4] + P_[4][4]*dt);
    nextP[0][8] = P_[0][8] + P_[1][8]*SF[9] + P_[2][8]*SF[11] + P_[3][8]*SF[10] + P_[10][8]*SF[14] + P_[11][8]*SF[15] + P_[12][8]*SPP[10] + dt*(P_[0][5] + P_[1][5]*SF[9] + P_[2][5]*SF[11] + P_[3][5]*SF[10] + P_[10][5]*SF[14] + P_[11][5]*SF[15] + P_[12][5]*SPP[10]);
    nextP[1][8] = P_[1][8] + P_[0][8]*SF[8] + P_[2][8]*SF[7] + P_[3][8]*SF[11] - P_[12][8]*SF[15] + P_[11][8]*SPP[10] - (P_[10][8]*q0)/2 + dt*(P_[1][5] + P_[0][5]*SF[8] + P_[2][5]*SF[7] + P_[3][5]*SF[11] - P_[12][5]*SF[15] + P_[11][5]*SPP[10] - (P_[10][5]*q0)/2);
    nextP[2][8] = P_[2][8] + P_[0][8]*SF[6] + P_[1][8]*SF[10] + P_[3][8]*SF[8] + P_[12][8]*SF[14] - P_[10][8]*SPP[10] - (P_[11][8]*q0)/2 + dt*(P_[2][5] + P_[0][5]*SF[6] + P_[1][5]*SF[10] + P_[3][5]*SF[8] + P_[12][5]*SF[14] - P_[10][5]*SPP[10] - (P_[11][5]*q0)/2);
    nextP[3][8] = P_[3][8] + P_[0][8]*SF[7] + P_[1][8]*SF[6] + P_[2][8]*SF[9] + P_[10][8]*SF[15] - P_[11][8]*SF[14] - (P_[12][8]*q0)/2 + dt*(P_[3][5] + P_[0][5]*SF[7] + P_[1][5]*SF[6] + P_[2][5]*SF[9] + P_[10][5]*SF[15] - P_[11][5]*SF[14] - (P_[12][5]*q0)/2);
    nextP[4][8] = P_[4][8] + P_[0][8]*SF[5] + P_[1][8]*SF[3] - P_[3][8]*SF[4] + P_[2][8]*SPP[0] + P_[13][8]*SPP[3] + P_[14][8]*SPP[6] - P_[15][8]*SPP[9] + dt*(P_[4][5] + P_[0][5]*SF[5] + P_[1][5]*SF[3] - P_[3][5]*SF[4] + P_[2][5]*SPP[0] + P_[13][5]*SPP[3] + P_[14][5]*SPP[6] - P_[15][5]*SPP[9]);
    nextP[5][8] = P_[5][8] + P_[0][8]*SF[4] + P_[2][8]*SF[3] + P_[3][8]*SF[5] - P_[1][8]*SPP[0] - P_[13][8]*SPP[8] + P_[14][8]*SPP[2] + P_[15][8]*SPP[5] + dt*(P_[5][5] + P_[0][5]*SF[4] + P_[2][5]*SF[3] + P_[3][5]*SF[5] - P_[1][5]*SPP[0] - P_[13][5]*SPP[8] + P_[14][5]*SPP[2] + P_[15][5]*SPP[5]);
    nextP[6][8] = P_[6][8] + P_[1][8]*SF[4] - P_[2][8]*SF[5] + P_[3][8]*SF[3] + P_[0][8]*SPP[0] + P_[13][8]*SPP[4] - P_[14][8]*SPP[7] - P_[15][8]*SPP[1] + dt*(P_[6][5] + P_[1][5]*SF[4] - P_[2][5]*SF[5] + P_[3][5]*SF[3] + P_[0][5]*SPP[0] + P_[13][5]*SPP[4] - P_[14][5]*SPP[7] - P_[15][5]*SPP[1]);
    nextP[7][8] = P_[7][8] + P_[4][8]*dt + dt*(P_[7][5] + P_[4][5]*dt);
    nextP[8][8] = P_[8][8] + P_[5][8]*dt + dt*(P_[8][5] + P_[5][5]*dt);
    nextP[0][9] = P_[0][9] + P_[1][9]*SF[9] + P_[2][9]*SF[11] + P_[3][9]*SF[10] + P_[10][9]*SF[14] + P_[11][9]*SF[15] + P_[12][9]*SPP[10] + dt*(P_[0][6] + P_[1][6]*SF[9] + P_[2][6]*SF[11] + P_[3][6]*SF[10] + P_[10][6]*SF[14] + P_[11][6]*SF[15] + P_[12][6]*SPP[10]);
    nextP[1][9] = P_[1][9] + P_[0][9]*SF[8] + P_[2][9]*SF[7] + P_[3][9]*SF[11] - P_[12][9]*SF[15] + P_[11][9]*SPP[10] - (P_[10][9]*q0)/2 + dt*(P_[1][6] + P_[0][6]*SF[8] + P_[2][6]*SF[7] + P_[3][6]*SF[11] - P_[12][6]*SF[15] + P_[11][6]*SPP[10] - (P_[10][6]*q0)/2);
    nextP[2][9] = P_[2][9] + P_[0][9]*SF[6] + P_[1][9]*SF[10] + P_[3][9]*SF[8] + P_[12][9]*SF[14] - P_[10][9]*SPP[10] - (P_[11][9]*q0)/2 + dt*(P_[2][6] + P_[0][6]*SF[6] + P_[1][6]*SF[10] + P_[3][6]*SF[8] + P_[12][6]*SF[14] - P_[10][6]*SPP[10] - (P_[11][6]*q0)/2);
    nextP[3][9] = P_[3][9] + P_[0][9]*SF[7] + P_[1][9]*SF[6] + P_[2][9]*SF[9] + P_[10][9]*SF[15] - P_[11][9]*SF[14] - (P_[12][9]*q0)/2 + dt*(P_[3][6] + P_[0][6]*SF[7] + P_[1][6]*SF[6] + P_[2][6]*SF[9] + P_[10][6]*SF[15] - P_[11][6]*SF[14] - (P_[12][6]*q0)/2);
    nextP[4][9] = P_[4][9] + P_[0][9]*SF[5] + P_[1][9]*SF[3] - P_[3][9]*SF[4] + P_[2][9]*SPP[0] + P_[13][9]*SPP[3] + P_[14][9]*SPP[6] - P_[15][9]*SPP[9] + dt*(P_[4][6] + P_[0][6]*SF[5] + P_[1][6]*SF[3] - P_[3][6]*SF[4] + P_[2][6]*SPP[0] + P_[13][6]*SPP[3] + P_[14][6]*SPP[6] - P_[15][6]*SPP[9]);
    nextP[5][9] = P_[5][9] + P_[0][9]*SF[4] + P_[2][9]*SF[3] + P_[3][9]*SF[5] - P_[1][9]*SPP[0] - P_[13][9]*SPP[8] + P_[14][9]*SPP[2] + P_[15][9]*SPP[5] + dt*(P_[5][6] + P_[0][6]*SF[4] + P_[2][6]*SF[3] + P_[3][6]*SF[5] - P_[1][6]*SPP[0] - P_[13][6]*SPP[8] + P_[14][6]*SPP[2] + P_[15][6]*SPP[5]);
    nextP[6][9] = P_[6][9] + P_[1][9]*SF[4] - P_[2][9]*SF[5] + P_[3][9]*SF[3] + P_[0][9]*SPP[0] + P_[13][9]*SPP[4] - P_[14][9]*SPP[7] - P_[15][9]*SPP[1] + dt*(P_[6][6] + P_[1][6]*SF[4] - P_[2][6]*SF[5] + P_[3][6]*SF[3] + P_[0][6]*SPP[0] + P_[13][6]*SPP[4] - P_[14][6]*SPP[7] - P_[15][6]*SPP[1]);
    nextP[7][9] = P_[7][9] + P_[4][9]*dt + dt*(P_[7][6] + P_[4][6]*dt);
    nextP[8][9] = P_[8][9] + P_[5][9]*dt + dt*(P_[8][6] + P_[5][6]*dt);
    nextP[9][9] = P_[9][9] + P_[6][9]*dt + dt*(P_[9][6] + P_[6][6]*dt);
    nextP[0][10] = P_[0][10] + P_[1][10]*SF[9] + P_[2][10]*SF[11] + P_[3][10]*SF[10] + P_[10][10]*SF[14] + P_[11][10]*SF[15] + P_[12][10]*SPP[10];
    nextP[1][10] = P_[1][10] + P_[0][10]*SF[8] + P_[2][10]*SF[7] + P_[3][10]*SF[11] - P_[12][10]*SF[15] + P_[11][10]*SPP[10] - (P_[10][10]*q0)/2;
    nextP[2][10] = P_[2][10] + P_[0][10]*SF[6] + P_[1][10]*SF[10] + P_[3][10]*SF[8] + P_[12][10]*SF[14] - P_[10][10]*SPP[10] - (P_[11][10]*q0)/2;
    nextP[3][10] = P_[3][10] + P_[0][10]*SF[7] + P_[1][10]*SF[6] + P_[2][10]*SF[9] + P_[10][10]*SF[15] - P_[11][10]*SF[14] - (P_[12][10]*q0)/2;
    nextP[4][10] = P_[4][10] + P_[0][10]*SF[5] + P_[1][10]*SF[3] - P_[3][10]*SF[4] + P_[2][10]*SPP[0] + P_[13][10]*SPP[3] + P_[14][10]*SPP[6] - P_[15][10]*SPP[9];
    nextP[5][10] = P_[5][10] + P_[0][10]*SF[4] + P_[2][10]*SF[3] + P_[3][10]*SF[5] - P_[1][10]*SPP[0] - P_[13][10]*SPP[8] + P_[14][10]*SPP[2] + P_[15][10]*SPP[5];
    nextP[6][10] = P_[6][10] + P_[1][10]*SF[4] - P_[2][10]*SF[5] + P_[3][10]*SF[3] + P_[0][10]*SPP[0] + P_[13][10]*SPP[4] - P_[14][10]*SPP[7] - P_[15][10]*SPP[1];
    nextP[7][10] = P_[7][10] + P_[4][10]*dt;
    nextP[8][10] = P_[8][10] + P_[5][10]*dt;
    nextP[9][10] = P_[9][10] + P_[6][10]*dt;
    nextP[10][10] = P_[10][10];
    nextP[0][11] = P_[0][11] + P_[1][11]*SF[9] + P_[2][11]*SF[11] + P_[3][11]*SF[10] + P_[10][11]*SF[14] + P_[11][11]*SF[15] + P_[12][11]*SPP[10];
    nextP[1][11] = P_[1][11] + P_[0][11]*SF[8] + P_[2][11]*SF[7] + P_[3][11]*SF[11] - P_[12][11]*SF[15] + P_[11][11]*SPP[10] - (P_[10][11]*q0)/2;
    nextP[2][11] = P_[2][11] + P_[0][11]*SF[6] + P_[1][11]*SF[10] + P_[3][11]*SF[8] + P_[12][11]*SF[14] - P_[10][11]*SPP[10] - (P_[11][11]*q0)/2;
    nextP[3][11] = P_[3][11] + P_[0][11]*SF[7] + P_[1][11]*SF[6] + P_[2][11]*SF[9] + P_[10][11]*SF[15] - P_[11][11]*SF[14] - (P_[12][11]*q0)/2;
    nextP[4][11] = P_[4][11] + P_[0][11]*SF[5] + P_[1][11]*SF[3] - P_[3][11]*SF[4] + P_[2][11]*SPP[0] + P_[13][11]*SPP[3] + P_[14][11]*SPP[6] - P_[15][11]*SPP[9];
    nextP[5][11] = P_[5][11] + P_[0][11]*SF[4] + P_[2][11]*SF[3] + P_[3][11]*SF[5] - P_[1][11]*SPP[0] - P_[13][11]*SPP[8] + P_[14][11]*SPP[2] + P_[15][11]*SPP[5];
    nextP[6][11] = P_[6][11] + P_[1][11]*SF[4] - P_[2][11]*SF[5] + P_[3][11]*SF[3] + P_[0][11]*SPP[0] + P_[13][11]*SPP[4] - P_[14][11]*SPP[7] - P_[15][11]*SPP[1];
    nextP[7][11] = P_[7][11] + P_[4][11]*dt;
    nextP[8][11] = P_[8][11] + P_[5][11]*dt;
    nextP[9][11] = P_[9][11] + P_[6][11]*dt;
    nextP[10][11] = P_[10][11];
    nextP[11][11] = P_[11][11];
    nextP[0][12] = P_[0][12] + P_[1][12]*SF[9] + P_[2][12]*SF[11] + P_[3][12]*SF[10] + P_[10][12]*SF[14] + P_[11][12]*SF[15] + P_[12][12]*SPP[10];
    nextP[1][12] = P_[1][12] + P_[0][12]*SF[8] + P_[2][12]*SF[7] + P_[3][12]*SF[11] - P_[12][12]*SF[15] + P_[11][12]*SPP[10] - (P_[10][12]*q0)/2;
    nextP[2][12] = P_[2][12] + P_[0][12]*SF[6] + P_[1][12]*SF[10] + P_[3][12]*SF[8] + P_[12][12]*SF[14] - P_[10][12]*SPP[10] - (P_[11][12]*q0)/2;
    nextP[3][12] = P_[3][12] + P_[0][12]*SF[7] + P_[1][12]*SF[6] + P_[2][12]*SF[9] + P_[10][12]*SF[15] - P_[11][12]*SF[14] - (P_[12][12]*q0)/2;
    nextP[4][12] = P_[4][12] + P_[0][12]*SF[5] + P_[1][12]*SF[3] - P_[3][12]*SF[4] + P_[2][12]*SPP[0] + P_[13][12]*SPP[3] + P_[14][12]*SPP[6] - P_[15][12]*SPP[9];
    nextP[5][12] = P_[5][12] + P_[0][12]*SF[4] + P_[2][12]*SF[3] + P_[3][12]*SF[5] - P_[1][12]*SPP[0] - P_[13][12]*SPP[8] + P_[14][12]*SPP[2] + P_[15][12]*SPP[5];
    nextP[6][12] = P_[6][12] + P_[1][12]*SF[4] - P_[2][12]*SF[5] + P_[3][12]*SF[3] + P_[0][12]*SPP[0] + P_[13][12]*SPP[4] - P_[14][12]*SPP[7] - P_[15][12]*SPP[1];
    nextP[7][12] = P_[7][12] + P_[4][12]*dt;
    nextP[8][12] = P_[8][12] + P_[5][12]*dt;
    nextP[9][12] = P_[9][12] + P_[6][12]*dt;
    nextP[10][12] = P_[10][12];
    nextP[11][12] = P_[11][12];
    nextP[12][12] = P_[12][12];
    
    // add process noise that is not from the IMU
    for (unsigned i = 0; i <= 12; i++) {
      nextP[i][i] += process_noise[i];
    }
    
    // Inhibit delta velocity bias learning by zeroing the covariance terms
    zeroRows(nextP,13,15);
    zeroCols(nextP,13,15);
        
    // stop position covariance growth if our total position variance reaches 100m
    if ((P_[7][7] + P_[8][8]) > 1.0e4) {
      for (uint8_t i = 7; i <= 8; i++) {
        for (uint8_t j = 0; j < k_num_states_; j++) {
          nextP[i][j] = P_[i][j];
          nextP[j][i] = P_[j][i];
        }
      }
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
    
    fixCovarianceErrors();
  }
  
  void ESKF::controlHeightSensorTimeouts() {
    /*
     * Handle the case where we have not fused height measurements recently and
     * uncertainty exceeds the max allowable. Reset using the best available height
     * measurement source, continue using it after the reset and declare the current
     * source failed if we have switched.
    */

    // check if height has been inertial deadreckoning for too long
    bool hgt_fusion_timeout = ((time_last_imu_ - time_last_hgt_fuse_) > 5e6);

    // reset the vertical position and velocity states
    if ((P_[9][9] > sq(hgt_reset_lim)) && (hgt_fusion_timeout)) {
      // boolean that indicates we will do a height reset
      bool reset_height = false;

      // handle the case where we are using external vision data for height
      if (ev_hgt_) {
        // check if vision data is available
        extVisionSample ev_init = _ext_vision_buffer.get_newest();
        bool ev_data_available = ((time_last_imu_ - ev_init.time_us) < 2 * EV_MAX_INTERVAL);

        // reset to ev data if it is available
        bool reset_to_ev = ev_data_available;

        if (reset_to_ev) {
          // request a reset
          reset_height = true;
        } else {
          // we have nothing to reset to
          reset_height = false;
        }
      }

      // Reset vertical position and velocity states to the last measurement
      if (reset_height) {
        resetHeight();
        // Reset the timout timer
        time_last_hgt_fuse_ = time_last_imu_;
      }
    }
  }
  
  void ESKF::resetHeight() {
    // reset the vertical position
    if (ev_hgt_) {
      // initialize vertical position with newest measurement
      extVisionSample ev_newest = _ext_vision_buffer.get_newest();

      // use the most recent data if it's time offset from the fusion time horizon is smaller
      int32_t dt_newest = ev_newest.time_us - _imu_sample_delayed.time_us;
      int32_t dt_delayed = ev_sample_delayed_.time_us - _imu_sample_delayed.time_us;

      if (std::abs(dt_newest) < std::abs(dt_delayed)) {
        state_.pos(2) = ev_newest.posNED(2);
      } else {
        state_.pos(2) = ev_sample_delayed_.posNED(2);
      }
    }

    // reset the vertical velocity covariance values
    zeroRows(P_, 6, 6);
    zeroCols(P_, 6, 6);

    // we don't know what the vertical velocity is, so set it to zero
    state_.vel(2) = 0.0f;

    // Set the variance to a value large enough to allow the state to converge quickly
    // that does not destabilise the filter
    P_[6][6] = 10.0f;
  }
    
  void ESKF::fuseVelPosHeight() {
    if(!fuse_) return;
    
    bool fuse_map[6] = {}; // map of booleans true when [VN,VE,VD,PN,PE,PD] observations are available
    bool innov_check_pass_map[6] = {}; // true when innovations consistency checks pass for [PN,PE,PD] observations
    scalar_t R[6] = {}; // observation variances for [VN,VE,VD,PN,PE,PD]
    scalar_t gate_size[6] = {}; // innovation consistency check gate sizes for [VN,VE,VD,PN,PE,PD] observations
    scalar_t Kfusion[k_num_states_] = {}; // Kalman gain vector for any single observation - sequential fusion is used
    
    // calculate innovations, innovations gate sizes and observation variances
    if (fuse_hor_vel_) {
      // enable fusion for NE velocity axes
      fuse_map[0] = fuse_map[1] = true;
      _velObsVarNE(1) = _velObsVarNE(0) = sq(fmaxf(gps_sample_delayed_.sacc, gps_vel_noise));
      
      // Set observation noise variance and innovation consistency check gate size for the NE position observations
      R[0] = _velObsVarNE(0);
      R[1] = _velObsVarNE(1);
      
      _hvelInnovGate = fmaxf(vel_innov_gate, 1.0f);
      gate_size[1] = gate_size[0] = _hvelInnovGate;
    }

    if (fuse_vert_vel_) {
      fuse_map[2] = true;
      // observation variance - use receiver reported accuracy with parameter setting the minimum value
      R[2] = fmaxf(gps_vel_noise, 0.01f);
      // use scaled horizontal speed accuracy assuming typical ratio of VDOP/HDOP
      R[2] = 1.5f * fmaxf(R[2], gps_sample_delayed_.sacc);
      R[2] = R[2] * R[2];
      // innovation gate size
      gate_size[2] = fmaxf(vel_innov_gate, 1.0f);
    }
    
    if (fuse_pos_) {
      fuse_map[3] = fuse_map[4] = true;
      
      if(gps_pos_) {
	gps_pos_ = false;
	gps_vel_ = false;
	fuse_ = false;
	// calculate observation process noise
	scalar_t lower_limit = fmaxf(gps_pos_noise, 0.01f);
	scalar_t upper_limit = fmaxf(pos_noaid_noise, lower_limit);
	_posObsNoiseNE = constrain(gps_sample_delayed_.hacc, lower_limit, upper_limit);

	// calculate innovations
	vel_pos_innov_[0] = state_.vel(0) - gps_sample_delayed_.vel(0);
	vel_pos_innov_[1] = state_.vel(1) - gps_sample_delayed_.vel(1);
	vel_pos_innov_[2] = state_.vel(2) - gps_sample_delayed_.vel(2);
	vel_pos_innov_[3] = state_.pos(0) - gps_sample_delayed_.pos(0);
	vel_pos_innov_[4] = state_.pos(1) - gps_sample_delayed_.pos(1);

	// set innovation gate size
	gate_size[3] = fmaxf(posNE_innov_gate, 1.0f);

	R[3] = sq(_posObsNoiseNE);

      } else if (ev_pos_) {
        ev_pos_ = false;
        fuse_ = false;
        // calculate innovations
        // use the absolute position
        vel_pos_innov_[3] = state_.pos(0) - ev_sample_delayed_.posNED(0);
        vel_pos_innov_[4] = state_.pos(1) - ev_sample_delayed_.posNED(1);
        
        // observation 1-STD error
        R[3] = fmaxf(0.05f, 0.01f);

        // innovation gate size
        gate_size[3] = fmaxf(5.0f, 1.0f);

      } else {
        // No observations - use a static position to constrain drift
        if (in_air_) {
          R[3] = fmaxf(10.0f, 0.5f);
        } else {
          R[3] = 0.5f;
        }
        vel_pos_innov_[3] = state_.pos(0) - last_known_posNED_(0);
        vel_pos_innov_[4] = state_.pos(1) - last_known_posNED_(1);

        // glitch protection is not required so set gate to a large value
        gate_size[3] = 100.0f;

        vel_pos_innov_[5] = state_.pos(2) - last_known_posNED_(2);
        fuse_map[5] = true;
        R[5] = 0.5f;
        R[5] = R[5] * R[5];
        // innovation gate size
        gate_size[5] = 100.0f;
      }

      // convert North position noise to variance
      R[3] = R[3] * R[3];

      // copy North axis values to East axis
      R[4] = R[3];
      gate_size[4] = gate_size[3];
    }
    
    if (fuse_height_) {
      fuse_map[5] = true;
      // calculate the innovation assuming the external vision observaton is in local NED frame
      vel_pos_innov_[5] = state_.pos(2) - ev_sample_delayed_.posNED(2);
      // observation variance - defined externally
      R[5] = fmaxf(0.05f, 0.01f);
      R[5] = R[5] * R[5];
      // innovation gate size
      gate_size[5] = fmaxf(5.0f, 1.0f);
    }

    // calculate innovation test ratios
    for (unsigned obs_index = 0; obs_index < 6; obs_index++) {
      if (fuse_map[obs_index]) {
        // compute the innovation variance SK = HPH + R
        unsigned state_index = obs_index + 4;	// we start with vx and this is the 4. state
        vel_pos_innov_var_[obs_index] = P_[state_index][state_index] + R[obs_index];
        // Compute the ratio of innovation to gate size
        vel_pos_test_ratio_[obs_index] = sq(vel_pos_innov_[obs_index]) / (sq(gate_size[obs_index]) * vel_pos_innov_var_[obs_index]);
      }
    }
    
    // check position, velocity and height innovations
    // treat 2D position and height as separate sensors
    bool pos_check_pass = ((vel_pos_test_ratio_[3] <= 1.0f) && (vel_pos_test_ratio_[4] <= 1.0f));
    innov_check_pass_map[3] = innov_check_pass_map[4] = pos_check_pass;
    innov_check_pass_map[5] = (vel_pos_test_ratio_[5] <= 1.0f);
    
    for (unsigned obs_index = 0; obs_index < 6; obs_index++) {
      // skip fusion if not requested or checks have failed
      if (!fuse_map[obs_index] || !innov_check_pass_map[obs_index]) {
        continue;
      }

      unsigned state_index = obs_index + 4;	// we start with vx and this is the 4. state

      // calculate kalman gain K = PHS, where S = 1/innovation variance
      for (int row = 0; row < k_num_states_; row++) {
        Kfusion[row] = P_[row][state_index] / vel_pos_innov_var_[obs_index];
      }

      // update covarinace matrix via Pnew = (I - KH)P
      float KHP[k_num_states_][k_num_states_];
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
        fixCovarianceErrors();

        // apply the state corrections
        fuse(Kfusion, vel_pos_innov_[obs_index]);
      }
    }
  }
  
  void ESKF::updateVision(const quat& q, const vec3& p, uint64_t time_usec, scalar_t dt) {
    // q here is rotation from enu to ros body
    quat q_nb = (q_rb.conjugate() * q.conjugate() * q_ne.conjugate()).conjugate();
    q_nb.normalize();

    if(reset_quat_){
      // reset the yaw angle to the value from the observaton quaternion
      // get the roll, pitch, yaw estimates from the quaternion states
      vec3 euler_init = dcm2vec(quat2dcm(state_.quat_nominal));

      // get initial yaw from the observation quaternion
      vec3 euler_obs = dcm2vec(quat2dcm(q_nb));
      euler_init(2) = euler_obs(2);

      // calculate initial quaternion states for the ekf
      state_.quat_nominal = from_axis_angle(euler_init);

      reset_quat_ = false;
    }

    // convert position from enu to ned
    vec3 pos_nb = q_ne.conjugate().toRotationMatrix() * p;
    
    // limit data rate to prevent data being lost
    if (time_usec - time_last_ext_vision > min_obs_interval_us) {
      extVisionSample ev_sample_new;
      // calculate the system time-stamp for the mid point of the integration period
      // copy required data
      ev_sample_new.angErr = 0.05f;
      ev_sample_new.posErr = 0.05f;
      ev_sample_new.quatNED = q_nb;
      ev_sample_new.posNED = pos_nb;
      ev_sample_new.time_us = time_usec - ev_delay_ms * 1000;
      time_last_ext_vision = time_usec;
      // push to buffer
      _ext_vision_buffer.push(ev_sample_new);
    }
    if(_ext_vision_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &ev_sample_delayed_)) {
      ev_pos_ = true;
      ev_yaw_ = true;
      fuse_ = true;
      ev_hgt_ = true;
    } else {
      ev_pos_ = false;
      ev_yaw_ = false;
      fuse_ = false;
      ev_hgt_ = false;
    }
    if(state_.pos(2) > 1.0) {
      in_air_ = true;
    } else {
      in_air_ = false;
    }
  }
  
  void ESKF::updateGps(const vec3& v, const vec3& p, uint64_t time_us, scalar_t dt) {
    vec3 vel_nb = q_ne.conjugate().toRotationMatrix() * v;
    vec3 pos_nb = q_ne.conjugate().toRotationMatrix() * p;
    
    // check for arrival of new sensor data at the fusion time horizon
    if (time_us - time_last_gps > min_obs_interval_us) {
      gpsSample gps_sample_new;
      gps_sample_new.time_us = time_us - 110 * 1000;

      gps_sample_new.time_us -= FILTER_UPDATE_PERIOD_MS * 1000 / 2;
      time_last_gps = time_us;

      gps_sample_new.time_us = max(gps_sample_new.time_us, _imu_sample_delayed.time_us);
      gps_sample_new.vel = vel_nb;
      gps_sample_new.hacc = 5.0;
      gps_sample_new.vacc = 8.0;
      gps_sample_new.sacc = 1.0;
      
      gps_sample_new.pos(0) = pos_nb(0);
      gps_sample_new.pos(1) = pos_nb(1);
      _gps_buffer.push(gps_sample_new);
    }
    if(_gps_buffer.pop_first_older_than(_imu_sample_delayed.time_us, &_gps_sample_delayed)) {
      gps_pos_ = true;
      gps_vel_ = true;
      fuse_ = true;      
    } else {
      gps_pos_ = false;
      gps_vel_ = false;
      fuse_ = false; 
    }
    if(state_.pos(2) > 1.0) {
      in_air_ = true;
    } else {
      in_air_ = false;
    }
  }
  
  void ESKF::fuseHeading() {
    if(!ev_yaw_) return;
    ev_yaw_ = false;
    
    // assign intermediate state variables
    scalar_t q0 = state_.quat_nominal.w();
    scalar_t q1 = state_.quat_nominal.x();
    scalar_t q2 = state_.quat_nominal.y();
    scalar_t q3 = state_.quat_nominal.z();
    
    scalar_t predicted_hdg, measured_hdg;
    scalar_t H_YAW[4];
    
    // update transformation matrix from body to world frame
    mat3 R_to_earth = quat_to_invrotmat(state_.quat_nominal);
    // determine if a 321 or 312 Euler sequence is best
    if (fabsf(R_to_earth(2, 0)) < fabsf(R_to_earth(2, 1))) {
      // calculate observation jacobian when we are observing the first rotation in a 321 sequence
      scalar_t t9 = q0*q3;
      scalar_t t10 = q1*q2;
      scalar_t t2 = t9+t10;
      scalar_t t3 = q0*q0;
      scalar_t t4 = q1*q1;
      scalar_t t5 = q2*q2;
      scalar_t t6 = q3*q3;
      scalar_t t7 = t3+t4-t5-t6;
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

      H_YAW[0] = t8*t14*(q3*t3-q3*t4+q3*t5+q3*t6+q0*q1*q2*2.0f)*-2.0f;
      H_YAW[1] = t8*t14*(-q2*t3+q2*t4+q2*t5+q2*t6+q0*q1*q3*2.0f)*-2.0f;
      H_YAW[2] = t8*t14*(q1*t3+q1*t4+q1*t5-q1*t6+q0*q2*q3*2.0f)*2.0f;
      H_YAW[3] = t8*t14*(q0*t3+q0*t4-q0*t5+q0*t6+q1*q2*q3*2.0f)*2.0f;

      // rotate the magnetometer measurement into earth frame
      vec3 euler321 = dcm2vec(quat2dcm(state_.quat_nominal));
      predicted_hdg = euler321(2); // we will need the predicted heading to calculate the innovation

      // calculate the yaw angle for a 321 sequence
      // Expressions obtained from yaw_input_321.c produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw321.m
      scalar_t Tbn_1_0 = 2.0f*(ev_sample_delayed_.quatNED.w() * ev_sample_delayed_.quatNED.z() + ev_sample_delayed_.quatNED.x() * ev_sample_delayed_.quatNED.y());
      scalar_t Tbn_0_0 = sq(ev_sample_delayed_.quatNED.w()) + sq(ev_sample_delayed_.quatNED.x()) - sq(ev_sample_delayed_.quatNED.y()) - sq(ev_sample_delayed_.quatNED.z());
      measured_hdg = atan2f(Tbn_1_0,Tbn_0_0);
    }
    else {
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

      scalar_t yaw = atan2f(-R_to_earth(0, 1), R_to_earth(1, 1)); // first rotation (yaw)

      predicted_hdg = yaw; // we will need the predicted heading to calculate the innovation
      // calculate the yaw angle for a 312 sequence
      // Values from yaw_input_312.c file produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
      scalar_t Tbn_0_1_neg = 2.0f * (ev_sample_delayed_.quatNED.w() * ev_sample_delayed_.quatNED.z() - ev_sample_delayed_.quatNED.x() * ev_sample_delayed_.quatNED.y());
      scalar_t Tbn_1_1 = sq(ev_sample_delayed_.quatNED.w()) - sq(ev_sample_delayed_.quatNED.x()) + sq(ev_sample_delayed_.quatNED.y()) - sq(ev_sample_delayed_.quatNED.z());
      measured_hdg = atan2f(Tbn_0_1_neg, Tbn_1_1);
    }
    
    scalar_t R_YAW = sq(fmaxf(ev_sample_delayed_.angErr, 1.0e-2f));
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
      initialiseCovariance();
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
      if(in_air_)
        return;
      else {
        // constrain the innovation to the maximum set by the gate
        scalar_t gate_limit = sqrtf(sq(heading_innov_gate) * heading_innov_var);
        heading_innov = constrain(heading_innov, -gate_limit, gate_limit);
      }
    }
    
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
      fixCovarianceErrors();

      // apply the state corrections
      fuse(Kfusion, heading_innov);
    }
  }
  
  // calculate the inverse rotation matrix from a quaternion rotation
  mat3 ESKF::quat_to_invrotmat(const quat &q) {
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
  
  void ESKF::fixCovarianceErrors() {
    scalar_t P_lim[4] = {};
    P_lim[0] = 1.0f;		// quaternion max var
    P_lim[1] = 1e6f;		// velocity max var
    P_lim[2] = 1e6f;		// positiion max var
    P_lim[3] = 1.0f;		// gyro bias max var
        
    for (int i = 0; i <= 3; i++) {
      // quaternion states
      P_[i][i] = constrain(P_[i][i], 0.0f, P_lim[0]);
    }
    
    for (int i = 4; i <= 6; i++) {
      // NED velocity states
      P_[i][i] = constrain(P_[i][i], 0.0f, P_lim[1]);
    }

    for (int i = 7; i <= 9; i++) {
      // NED position states
      P_[i][i] = constrain(P_[i][i], 0.0f, P_lim[2]);
    }
    
    for (int i = 10; i <= 12; i++) {
      // gyro bias states
      P_[i][i] = constrain(P_[i][i], 0.0f, P_lim[3]);
    }
    
    // force symmetry on the quaternion, velocity, positon and gyro bias state covariances
    makeSymmetrical(P_,0,12);
    
    // accelerometer bias states zeros (inhibited)
    zeroRows(P_,13,15);
    zeroCols(P_,13,15);
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
      state_.vel(i) = state_.vel(i) - K[i + 4] * innovation;
    }

    for (unsigned i = 0; i < 3; i++) {
      state_.pos(i) = state_.pos(i) - K[i + 7] * innovation;
    }

    for (unsigned i = 0; i < 3; i++) {
      state_.gyro_bias(i) = state_.gyro_bias(i) - K[i + 10] * innovation;
    }

    for (unsigned i = 0; i < 3; i++) {
      state_.accel_bias(i) = state_.accel_bias(i) - K[i + 13] * innovation;
    }
  }
  
  vec3 ESKF::getRPY(const mat3& R) {
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
    return q_ne.toRotationMatrix() * rpy;
  }

  quat ESKF::getQuat() { 
    quat q_eb = (q_rb.conjugate() * state_.quat_nominal.conjugate() * q_ne.conjugate()).conjugate();
    q_eb.normalize();
    return q_eb; 
  }

  vec3 ESKF::getXYZ() {
    return q_ne.toRotationMatrix() * state_.pos;
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
	
  void ESKF::constrainStates() {
    state_.quat_nominal.w() = constrain(state_.quat_nominal.w(), -1.0f, 1.0f);
    state_.quat_nominal.x() = constrain(state_.quat_nominal.x(), -1.0f, 1.0f);
    state_.quat_nominal.y() = constrain(state_.quat_nominal.y(), -1.0f, 1.0f);
    state_.quat_nominal.z() = constrain(state_.quat_nominal.z(), -1.0f, 1.0f);
	  
    for (int i = 0; i < 3; i++) {
      state_.vel(i) = constrain(state_.vel(i), -1000.0f, 1000.0f);
    }

    for (int i = 0; i < 3; i++) {
      state_.pos(i) = constrain(state_.pos(i), -1.e6f, 1.e6f);
    }

    for (int i = 0; i < 3; i++) {
      state_.gyro_bias(i) = constrain(state_.gyro_bias(i), -0.349066f * _dt_ekf_avg, 0.349066f * _dt_ekf_avg);
    }

    for (int i = 0; i < 3; i++) {
      state_.accel_bias(i) = constrain(state_.accel_bias(i), -acc_bias_lim * _dt_ekf_avg, acc_bias_lim * _dt_ekf_avg);
    }
  }
} //  namespace eskf
