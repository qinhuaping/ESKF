#include <eskf/Node.hpp>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <assert.h>

namespace eskf
{

Node::Node(const ros::NodeHandle &nh, const ros::NodeHandle &pnh) : nh_(pnh), init_(false)
{
	//  subscribe
	ROS_INFO("Subscribing to ~imu.");
	subImu_ = nh_.subscribe<sensor_msgs::Imu>("imu", 1000, &Node::inputCallback, this,
			ros::TransportHints().tcpNoDelay(true));
	ROS_INFO("Subscribing to ~pose.");
	subPOSE_ = nh_.subscribe("pose", 1, &Node::measurementCallback, this);

	pubPose_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("pose", 1);
}

void Node::inputCallback(const sensor_msgs::ImuConstPtr &imuMsg)
{
	vec3 wm;        //  measured angular rate
	vec3 am;        //  measured linear acceleration

	wm = vec3(imuMsg->angular_velocity.x, imuMsg->angular_velocity.y, imuMsg->angular_velocity.z);
	am = vec3(imuMsg->linear_acceleration.x, imuMsg->linear_acceleration.y, imuMsg->linear_acceleration.z);

	if (prevStampIMU_.sec != 0) {

		const double delta = (imuMsg->header.stamp - prevStampIMU_).toSec();

		if (!init_) {
			init_ = true;
			ROS_INFO("Initialized ESKF");
		}

		//  run kalman filter
		eskf_.predict(wm, am, static_cast<uint64_t>(imuMsg->header.stamp.toSec() * 1e6f), delta);

		const quat n2b = eskf_.getQuat();
		const vec3 position = eskf_.getXYZ();

		sensor_msgs::Imu imu = *imuMsg;
		imu.header.seq = 0;
		
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
		pubPose_.publish(msg_pose);
	}

	prevStampIMU_ = imuMsg->header.stamp;
}

void Node::measurementCallback(const geometry_msgs::PoseStampedConstPtr &poseMsg)
{
	if (prevStampPOSE_.sec != 0) {
		const double delta = (poseMsg->header.stamp - prevStampPOSE_).toSec();
		// get measurements
		quat z_q = quat(poseMsg->pose.orientation.w, poseMsg->pose.orientation.x, poseMsg->pose.orientation.y,
				poseMsg->pose.orientation.z);
		vec3 z_p = vec3(poseMsg->pose.position.x, poseMsg->pose.position.y, poseMsg->pose.position.z);
		eskf_.update(z_q, z_p, static_cast<uint64_t>(poseMsg->header.stamp.toSec() * 1e6f), delta);
	}

	prevStampPOSE_ = poseMsg->header.stamp;
}
}
