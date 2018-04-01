#ifndef ESKF_H_
#define ESKF_H_

#include <Eigen/Core>
#include <Eigen/Dense>

namespace eskf {

  template <typename data_type>
  class RingBuffer {
  public:
	  RingBuffer() {
		_buffer = NULL;
		_head = _tail = _size = 0;
		_first_write = true;
	}
	~RingBuffer() { delete[] _buffer; }

	bool allocate(int size) {
		if (size <= 0) {
			return false;
		}

		if (_buffer != NULL) {
			delete[] _buffer;
		}

		_buffer = new data_type[size];

		if (_buffer == NULL) {
			return false;
		}

		_size = size;
		// set the time elements to zero so that bad data is not retrieved from the buffers
		//for (unsigned index = 0; index < _size; index++) {
	  //  _buffer[index].time_us = 0;
		//}
		_first_write = true;
		return true;
	}

	void unallocate()	{
		if (_buffer != NULL) {
			delete[] _buffer;
		}
	}

	inline void push(data_type sample) {
		int head_new = _head;

		if (_first_write) {
			head_new = _head;

		} else {
			head_new = (_head + 1) % _size;
		}

		_buffer[head_new] = sample;
		_head = head_new;

		// move tail if we overwrite it
		if (_head == _tail && !_first_write) {
			_tail = (_tail + 1) % _size;

		} else {
			_first_write = false;
		}
	}

	inline data_type get_oldest() {
		return _buffer[_tail];
	}

	unsigned get_oldest_index() {
		return _tail;
	}

	inline data_type get_newest() {
		return _buffer[_head];
	}

	data_type &operator[](unsigned index) {
		return _buffer[index];
	}

	// return data at the specified index
	inline data_type get_from_index(unsigned index) {
		if (index >= _size) {
			index = _size-1;
		}
		return _buffer[index];
	}

	// push data to the specified index
	inline void push_to_index(unsigned index, data_type sample) {
		if (index >= _size) {
			index = _size-1;
		}
		_buffer[index] = sample;
	}

	// return the length of the buffer
	unsigned get_length() {
		return _size;
	}

  private:
    data_type *_buffer;
    unsigned _head, _tail, _size;
    bool _first_write;
  };

  class ESKF {
  public:
    static constexpr int k_num_states_ = 16;
	  typedef double scalar_t;
    typedef Eigen::Matrix<scalar_t, 3, 1> vec3; /// Vector in R3
    typedef Eigen::Matrix<scalar_t, k_num_states_, 1> vec7; /// Vector in R3
    typedef Eigen::Matrix<scalar_t, 3, 3> mat3; /// Matrix in R3
    typedef Eigen::Quaternion<scalar_t> quat;   /// Member of S4
    static const unsigned FILTER_UPDATE_PERIOD_MS = 12;	// ekf prediction period in milliseconds - this should ideally be an integer multiple of the IMU time delta
    
    ESKF();

    void predict(const vec3& w, const vec3& a, scalar_t dt);

    void update(const quat& q, scalar_t dt);
    
    const quat& getQuat() const;

    vec3 getRPY(const mat3& R);

  private:
    
    void constrainStates();
    void initialiseQuatCovariances(const vec3& rot_vec_var);
	  void zeroRows(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
    void zeroCols(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
    void makeSymmetrical(scalar_t (&cov_mat)[k_num_states_][k_num_states_], uint8_t first, uint8_t last);
	  void fuse(scalar_t *K, scalar_t innovation);
    void fixCovarianceErrors();
    void initialiseCovariance();
    void predictCovariance();
    //void updatePos(const vec3& p, scalar_t dt);
    //void updateYaw(const quat& q, scalar_t dt);
    mat3 quat_to_invrotmat(const quat &q);
    quat from_axis_angle(vec3 vec);
    quat from_axis_angle(const vec3 &axis, scalar_t theta);
    vec3 to_axis_angle(const quat& q);
    mat3 quat2dcm(const quat& q);

    /* State vector:
     * Attitude quaternion
     * Delta Angle bias - rad (X,Y,Z)
    */
    struct state {
		  quat	  quat_nominal;	  ///< quaternion defining the rotaton from NED to XYZ frame
		  vec3    vel;	          ///< NED velocity in earth frame in m/s
		  vec3    pos;	          ///< NED position in earth frame in m
		  vec3    gyro_bias;	    ///< delta angle bias estimate in rad
		  vec3    accel_bias;	    ///< delta velocity bias estimate in m/s
    } state_;
    
    struct imuSample {
      vec3    delta_ang;		///< delta angle in body frame (integrated gyro measurements) (rad)
      vec3    delta_vel;		///< delta velocity in body frame (integrated accelerometer measurements) (m/sec)
      scalar_t  delta_ang_dt;	///< delta angle integration period (sec)
      scalar_t  delta_vel_dt;	///< delta velocity integration period (sec)
    };

    imuSample _imu_sample_new{};		///< imu sample capturing the newest imu data
    imuSample _imu_down_sampled{};  ///< down sampled imu data (sensor rate -> filter update rate)
    RingBuffer<imuSample> _imu_buffer;

    quat _q_down_sampled;
    float _imu_collection_time_adj{0.0f};	///< the amount of time the IMU collection needs to be advanced to meet the target set by FILTER_UPDATE_PERIOD_MS (sec)
    imuSample _imu_sample_delayed{};	// captures the imu sample on the delayed time horizon

    scalar_t _dt_ekf_avg{0.001f * FILTER_UPDATE_PERIOD_MS}; ///< average update rate of the ekf
    
    bool collect_imu(imuSample& imu);
    
    scalar_t P_[k_num_states_][k_num_states_]; /// System covariance matrix

    static constexpr scalar_t kOneG = 9.80665;  /// Earth gravity (m/s^2)
    static constexpr scalar_t acc_bias_lim = 0.4; ///< maximum accel bias magnitude (m/sec**2)

    // process noise
    scalar_t gyro_bias_p_noise{1.0e-3};		///< process noise for IMU rate gyro bias prediction (rad/sec**2)
    scalar_t accel_bias_p_noise{3.0e-3};	///< process noise for IMU accelerometer bias prediction (m/sec**3)

    // input noise
    scalar_t gyro_noise{1.5e-2};		///< IMU angular rate noise used for covariance prediction (rad/sec)
    scalar_t accel_noise{3.5e-1};		///< IMU acceleration noise use for covariance prediction (m/sec**2)
    
	  // initialization errors
	  scalar_t switch_on_gyro_bias{0.1f};		///< 1-sigma gyro bias uncertainty at switch on (rad/sec)
	  scalar_t switch_on_accel_bias{0.2f};	///< 1-sigma accelerometer bias uncertainty at switch on (m/sec**2)
	  scalar_t initial_tilt_err{0.1f};		///< 1-sigma tilt error after initial alignment using gravity vector (rad)
    
    scalar_t vel_noise{5.0e-1f};	///< minimum allowed observation noise for velocity fusion (m/sec)
	  scalar_t pos_noise{0.5f};		///< minimum allowed observation noise for position fusion (m)
	  scalar_t baro_noise{2.0f};			///< observation noise for barometric height fusion (m)
    
    scalar_t yaw_err{0.05f};
    scalar_t heading_innov_gate{2.6f};		///< heading fusion innovation consistency gate size (STD)
    
    mat3 rosb2px4b; ///< rotation from ROS body to PX4 body frame. need conversion for IMU
    
    quat q_ne; ///< rotation from NED to ENU
    quat q_rb; ///< rotation from ROS body to PX4 body
    quat q_nb; ///< resulting rotation from ned to body
    
  };
}

#endif /* defined(ESKF_H_) */
