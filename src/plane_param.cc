#include "plane_param.h"

#include <Eigen/Dense>

PlaneParam::PlaneParam()
{}

PlaneParam::PlaneParam(
		const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2) {

	Eigen::Vector3d V0(v0[0] * v0[2], v0[1] * v0[2], v0[2]);
	Eigen::Vector3d V1(v1[0] * v1[2], v1[1] * v1[2], v1[2]);
	Eigen::Vector3d V2(v2[0] * v2[2], v2[1] * v2[2], v2[2]);
	/*
	Eigen::Vector3d V0(v0);
	Eigen::Vector3d V1(v1);
	Eigen::Vector3d V2(v2);
	*/
	Eigen::Vector3d norm = (V1 - V0).cross(V2 - V0);
	if (norm.norm() < 1e-10 || norm[2] == 0) {
		invalid_ = true;
	}
	norm.normalize();
	d_ = -V0.dot(norm);
	n1_ = norm[0];
	n2_ = norm[1];
	n3_ = norm[2];
	invalid_ = false;
	if (d_ == 0)
		invalid_ = true;
}

bool PlaneParam::ProjectiveIntersection(const PlaneParam& other, LineSegment* l) const {
	if (invalid_ || other.invalid_)
		return false;
	//figure out the intersection line direction
	l->n1 = n2_ * other.n3_ - n3_ * other.n2_;
	l->n2 = n3_ * other.n1_ - n1_ * other.n3_;
	l->n3 = n1_ * other.n2_ - n2_ * other.n1_;

	//figure out one intersection point
	if (l->n1 != 0) {
		l->x = 0;
		l->z = (other.n2_ * d_ - n2_ * other.d_) / l->n1;
		l->y = (-other.n3_ * d_ + n3_ * other.d_) / l->n1;
	}
	else if (l->n2 != 0) {
		l->y = 0;
		l->x = (other.n3_ * d_ - n3_ * other.d_) / l->n2;
		l->z = (-other.n1_ * d_ + n1_ * other.d_) / l->n2;
	}
	else {
		return false;
	}

	if (l->x == 0 && l->y == 0)
		return false;

	if (l->z == 0) {
		l->x += l->n1;
		l->y += l->n2;
		l->z += l->n3;
	}
	K nx, ny;
	if (l->z + l->n3 != 0) {
		nx = (l->x+l->n1) / (l->z+l->n3);
		ny = (l->y+l->n2) / (l->z+l->n3);
	} else if (l->z - l->n3 != 0) {
		nx = (l->x-l->n1) / (l->z-l->n3);
		ny = (l->y-l->n2) / (l->z-l->n3);
	} else {
		return false;
	}

	l->x /= l->z;
	l->y /= l->z;
	l->n1 = nx - l->x;
	l->n2 = ny - l->y;
	return true;
}
