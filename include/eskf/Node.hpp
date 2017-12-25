#ifndef ESKF_NODE_HPP_
#define ESKF_NODE_HPP_

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/PoseStamped.h>
#include <message_filters/subscriber.h>

#include <eskf/ESKF.hpp>

namespace eskf {

  class Node {
  public:
    Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh);
    
  private:
    typedef eskf::ESKF::vec3 vec3;
    typedef eskf::ESKF::quat quat;
    
    ros::NodeHandle nh_;
    
    // publishers
    ros::Publisher pubRPY_;
    
    //  subsribers
    ros::Subscriber subImu_;
    ros::Subscriber subPOSE_;
    
    // implementation
    eskf::ESKF eskf_;
    ros::Time prevStampIMU_;
    ros::Time prevStampPOSE_;
    bool init_;
    
    //  callbacks
    void inputCallback(const sensor_msgs::ImuConstPtr&);
    void measurementCallback(const geometry_msgs::PoseStampedConstPtr&);
  };
} //  namespace eskf

#endif // ESKF_NODE_HPP_
