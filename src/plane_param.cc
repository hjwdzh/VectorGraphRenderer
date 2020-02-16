#include "plane_param.h"

#include <Eigen/Dense>

PlaneParam::PlaneParam()
{}

PlaneParam::PlaneParam(
		const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2) {
	Eigen::Vector3d norm = (v1 - v0).cross(v2 - v0);
	if (norm.norm() < 1e-10 || norm[2] == 0) {
		invalid_ = true;
	}
	norm.normalize();
	d_ = v0.dot(norm);
	n1_ = norm[0];
	n2_ = norm[1];
	n3_ = norm[2];
	invalid_ = false;
}
