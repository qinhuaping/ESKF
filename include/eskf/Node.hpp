#ifndef ESKF_NODE_HPP_
#define ESKF_NODE_HPP_

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PoseStamped.h>
#include <message_filters/subscriber.h>
#include <eskf/ESKF.hpp>

namespace eskf {

  class Node {
  public:
    Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh);
    
  private:
    
    ros::NodeHandle nh_;
    
    // publishers
    ros::Publisher pubPose_;
    
    //  subsribers
    ros::Subscriber subImu_;
    ros::Subscriber subVisionPose_;
    ros::Subscriber subGpsPose_;
    
    // implementation
    eskf::ESKF eskf_;
    ros::Time prevStampImu_;
    ros::Time prevStampVisionPose_;
    ros::Time prevStampGpsPose_;
    bool init_;
    
    //  callbacks
    void inputCallback(const sensor_msgs::ImuConstPtr&);
    void visionCallback(const geometry_msgs::PoseStampedConstPtr&);
    void gpsCallback(const nav_msgs::OdometryConstPtr&);
  };
} //  namespace eskf

#endif // ESKF_NODE_HPP_
