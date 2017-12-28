#ifndef ESKF_H_
#define ESKF_H_

#include <Eigen/Core>
#include <Eigen/Dense>

namespace eskf {

  class ESKF {
  public:
    static constexpr int k_num_states_ = 7;
	  typedef double scalar_t;
    typedef Eigen::Matrix<scalar_t, 3, 1> vec3; /// Vector in R3
    typedef Eigen::Matrix<scalar_t, k_num_states_, 1> vec7; /// Vector in R3
    typedef Eigen::Matrix<scalar_t, 3, 3> mat3; /// Matrix in R3
    typedef Eigen::Quaternion<scalar_t> quat;   /// Member of S4
    
    ESKF();

    void predict(const vec3& w, scalar_t dt);

    void update(const quat& q, scalar_t dt);
    
    bool initialize(scalar_t dt);
    
    const quat& getQuat() const;

    vec3 getRPY(const mat3& R);

  private:
    
    void constrainStates(scalar_t dt);
    void initialiseQuatCovariances(const vec3& rot_vec_var);
	  void zeroRows(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
    void zeroCols(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
    void makeSymmetrical(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
	  void fuse(scalar_t *K, scalar_t innovation);
    void fixCovarianceErrors(scalar_t dt);
    void initialiseCovariance(scalar_t dt);
    void updateYaw(const quat& q, scalar_t dt);
    mat3 quat_to_invrotmat(const quat &q);
    quat from_axis_angle(vec3 vec);
    quat from_axis_angle(const vec3 &axis, scalar_t theta);
    /* State vector:
     * Attitude quaternion
     * Delta Angle bias - rad (X,Y,Z)
    */
    struct state {
		  quat	  quat_nominal;	  ///< quaternion defining the rotaton from NED to XYZ frame
		  vec3    gyro_bias;	    ///< delta angle bias estimate in rad
    } state_;
    
    scalar_t P_[k_num_states_][k_num_states_]; /// System covariance matrix
    
    vec7 dx_; /// State vector increment
    
    // process noise
    scalar_t gyro_bias_p_noise{1.0e-3};		///< process noise for IMU rate gyro bias prediction (rad/sec**2)
	
	  // input noise
	  scalar_t gyro_noise{1.5e-2};		///< IMU angular rate noise used for covariance prediction (rad/sec)
    
	  // initialization errors
	  scalar_t switch_on_gyro_bias{0.1f};		///< 1-sigma gyro bias uncertainty at switch on (rad/sec)
	  scalar_t initial_tilt_err{0.1f};		///< 1-sigma tilt error after initial alignment using gravity vector (rad)
    
    scalar_t yaw_err{0.05f};
    scalar_t heading_innov_gate{2.6f};		///< heading fusion innovation consistency gate size (STD)
    
    mat3 rosb2px4b; ///< rotation from ROS body to PX4 body frame. need conversion for IMU
    
    bool firstPredict = false;
    quat q_ne; ///< rotation from NED to ENU
    quat q_rb; ///< rotation from ROS body to PX4 body
    quat q_nb; ///< resulting rotation from ned to body
    
    scalar_t last_known_yaw, last_known_pitch, last_known_roll;
  };
}

#endif /* defined(ESKF_H_) */
