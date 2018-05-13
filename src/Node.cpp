#include <eskf/Node.hpp>
#include <geometry_msgs/PoseWithCovarianceStamped.h>

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
		//  advertise
		pubPose_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("pose", 1);
	}

	void Node::inputCallback(const sensor_msgs::ImuConstPtr &imuMsg) {
		vec3 wm;        //  measured angular rate
		vec3 am;        //  measured linear acceleration

		wm = vec3(imuMsg->angular_velocity.x, imuMsg->angular_velocity.y, imuMsg->angular_velocity.z);
		am = vec3(imuMsg->linear_acceleration.x, imuMsg->linear_acceleration.y, imuMsg->linear_acceleration.z);

		if (prevStampImu_.sec != 0) {

			const double delta = (imuMsg->header.stamp - prevStampImu_).toSec();

			if (!init_) {
				init_ = true;
				ROS_INFO("Initialized ESKF");
			}

			//  run kalman filter
			eskf_.run(wm, am, static_cast<uint64_t>(imuMsg->header.stamp.toSec() * 1e6f), delta);

			//	get kalman output
			const quat n2b = eskf_.getQuat();
			const vec3 position = eskf_.getXYZ();

			sensor_msgs::Imu imu = *imuMsg;
			imu.header.seq = 0;
			
			//	compose output topic
			geometry_msgs::PoseWithCovarianceStamped msg_pose;
			msg_pose.header = imu.header;
			msg_pose.pose.pose.position.x = position[0];
			msg_pose.pose.pose.position.y = position[1];
			msg_pose.pose.pose.position.z = position[2];
			msg_pose.pose.pose.orientation.x = n2b.x();
			msg_pose.pose.pose.orientation.y = n2b.y();
			msg_pose.pose.pose.orientation.z = n2b.z();
			msg_pose.pose.pose.orientation.w = n2b.w();
			//px4 doesn't use covariance for vision so set it up to zero 
			for(size_t i=0; i<36; ++i)
				msg_pose.pose.covariance[i] = 0;
			
			// 	publish output topic 
			pubPose_.publish(msg_pose);
		}

		prevStampImu_ = imuMsg->header.stamp;
	}

	void Node::visionCallback(const geometry_msgs::PoseStampedConstPtr &poseMsg) {
		if (prevStampVisionPose_.sec != 0) {
			const double delta = (poseMsg->header.stamp - prevStampVisionPose_).toSec();
			// get vision measurements
			quat z_q = quat(poseMsg->pose.orientation.w, poseMsg->pose.orientation.x, poseMsg->pose.orientation.y, poseMsg->pose.orientation.z);
			vec3 z_p = vec3(poseMsg->pose.position.x, poseMsg->pose.position.y, poseMsg->pose.position.z);
			eskf_.updateVision(z_q, z_p, static_cast<uint64_t>(poseMsg->header.stamp.toSec() * 1e6f), delta);
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
