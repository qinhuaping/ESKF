#include <eskf/Node.hpp>
#include <geometry_msgs/Vector3Stamped.h>

namespace eskf
{

Node::Node(const ros::NodeHandle &nh, const ros::NodeHandle &pnh) : nh_(pnh), init_(false) {
  //  subscribe
  ROS_INFO("Subscribing to ~imu.");
  subImu_ = nh_.subscribe<sensor_msgs::Imu>("imu", 1000, &Node::inputCallback, this, ros::TransportHints().tcpNoDelay(true));
  ROS_INFO("Subscribing to ~vision_pose.");
  subVisionPose_ = nh_.subscribe("vision_pose", 1, &Node::visionCallback, this);
  ROS_INFO("Subscribing to ~gps_pose.");
  subGpsPose_ = nh_.subscribe("gps_pose", 1, &Node::gpsCallback, this);
  
  pubRPY_ = nh_.advertise<geometry_msgs::Vector3Stamped>("rpy", 1);
  pubXYZ_ = nh_.advertise<geometry_msgs::Vector3Stamped>("xyz", 1);
  pubPose_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("pose", 1);
}

void Node::inputCallback(const sensor_msgs::ImuConstPtr &imuMsg) {

  vec3 wm = vec3(imuMsg->angular_velocity.x, imuMsg->angular_velocity.y, imuMsg->angular_velocity.z); //  measured angular rate
  vec3 am = vec3(imuMsg->linear_acceleration.x, imuMsg->linear_acceleration.y, imuMsg->linear_acceleration.z); //  measured linear acceleration

  if (prevStampImu_.sec != 0) {
    const double delta = (imuMsg->header.stamp - prevStampImu_).toSec();

    if (!init_) {
      init_ = true;
      ROS_INFO("Initialized ESKF");
    }

    //  run kalman filter
    eskf_.predict(wm, am, static_cast<uint64_t>(imuMsg->header.stamp.toSec()*1e6f), delta);

    // get kalman filter result
    const quat n2b = eskf_.getQuat();
    const vec3 orientation = eskf_.getRPY(n2b.matrix());
    const vec3 position = eskf_.getXYZ();

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

    // x-y-z (m)
    geometry_msgs::Vector3Stamped xyz;
    xyz.header = imu.header;
    xyz.vector.x = position[0];
    xyz.vector.y = position[1];
    xyz.vector.z = position[2];
    // publish our topics
    pubXYZ_.publish(xyz);

    geometry_msgs::PoseWithCovarianceStamped pose;
    pose.header = imu.header;
    pose.pose.pose.position.x = position[0];
    pose.pose.pose.position.y = position[1];
    pose.pose.pose.position.z = position[2];
    pose.pose.pose.orientation.x = n2b.x();
    pose.pose.pose.orientation.y = n2b.y();
    pose.pose.pose.orientation.z = n2b.z();
    pose.pose.pose.orientation.w = n2b.w();

    //px4 doesn't use covariance for vision so set it up to zero 
    for(size_t i = 0; i < 36; ++i)
      pose.pose.covariance[i] = 0;
    pubPose_.publish(pose);
  }
  prevStampImu_ = imuMsg->header.stamp;
}
  
void Node::visionCallback(const geometry_msgs::PoseWithCovarianceStampedConstPtr& poseMsg) {
  if(prevStampVisionPose_.sec != 0) {
    const double delta = (poseMsg->header.stamp - prevStampVisionPose_).toSec();
    // get measurements
    quat z_q = quat(poseMsg->pose.pose.orientation.w, poseMsg->pose.pose.orientation.x, poseMsg->pose.pose.orientation.y, poseMsg->pose.pose.orientation.z);
    vec3 z_p = vec3(poseMsg->pose.pose.position.x, poseMsg->pose.pose.position.y, poseMsg->pose.pose.position.z);
    eskf_.updateVision(z_q, z_p, static_cast<uint64_t>(poseMsg->header.stamp.toSec()*1e6f), delta);
  }
  prevStampVisionPose_ = poseMsg->header.stamp;
}

void Node::gpsCallback(const nav_msgs::OdometryConstPtr &odomMsg) {
  if (prevStampGpsPose_.sec != 0) {
    const double delta = (odomMsg->header.stamp - prevStampGpsPose_).toSec();
    // get gps measurements
    vec3 z_v = vec3(odomMsg->twist.twist.linear.x, odomMsg->twist.twist.linear.y, odomMsg->twist.twist.linear.z);
    vec3 z_p = vec3(odomMsg->pose.pose.position.x, odomMsg->pose.pose.position.y, odomMsg->pose.pose.position.z);
    eskf_.updateGps(z_v, z_p, static_cast<uint64_t>(odomMsg->header.stamp.toSec() * 1e6f), delta);
  }
  prevStampGpsPose_ = odomMsg->header.stamp;
}

}
