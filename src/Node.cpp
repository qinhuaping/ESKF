#include <eskf/Node.hpp>
#include <geometry_msgs/Vector3Stamped.h>
#include <assert.h>

namespace eskf {

  Node::Node(const ros::NodeHandle &nh, const ros::NodeHandle &pnh) : nh_(pnh), init_(false) {
    //  subscribe
    ROS_INFO("Subscribing to ~imu.");
    subImu_ = nh_.subscribe<sensor_msgs::Imu>("imu", 100, boost::bind(&Node::inputCallback, this, _1));
    ROS_INFO("Subscribing to ~pose.");
    subPOSE_ = nh_.subscribe("pose", 1, &Node::measurementCallback, this);
    
    pubRPY_ = nh_.advertise<geometry_msgs::Vector3Stamped>("rpy", 1);
  }

  void Node::inputCallback(const sensor_msgs::ImuConstPtr& imuMsg) {
    vec3 wm;        //  measured angular rate
            
    wm = vec3(imuMsg->angular_velocity.x, imuMsg->angular_velocity.y, imuMsg->angular_velocity.z);
         
    if (prevStampIMU_.sec != 0) {
      
      const double delta = (imuMsg->header.stamp - prevStampIMU_).toSec();
      
      if (!init_) {
	      eskf_.initialize(delta);
        init_ = true;
        ROS_INFO("Initialized ESKF");
      }
      
      //  run kalman filter
      eskf_.predict(wm, delta);
      
      const eskf::ESKF::quat n2b = eskf_.getQuat();
      const vec3 orientation = eskf_.getRPY(n2b.matrix());
      
      sensor_msgs::Imu imu = *imuMsg;
      imu.header.seq = 0;
            
      // roll-pitch-yaw (rad)
      geometry_msgs::Vector3Stamped rpy;
      rpy.header = imu.header;
      rpy.vector.x = orientation[0];
      rpy.vector.y = orientation[1];
      rpy.vector.z = orientation[2];
      // publish our topics
      pubRPY_.publish(rpy);
    }
    prevStampIMU_ = imuMsg->header.stamp;
  }
  
  void Node::measurementCallback(const geometry_msgs::PoseStampedConstPtr& poseMsg) {
    /*
    if(prevStampPOSE_.sec != 0) {
      const double delta = (poseMsg->header.stamp - prevStampPOSE_).toSec();
      // get measurements
      quat z_q = quat(poseMsg->pose.orientation.w, poseMsg->pose.orientation.x, poseMsg->pose.orientation.y, poseMsg->pose.orientation.z);
      eskf_.update(z_q, delta);
    }
    prevStampPOSE_ = poseMsg->header.stamp;
    */
  }
}
