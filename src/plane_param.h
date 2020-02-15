#ifndef VECTORGRAPH_RENDERER_PLANE_PARAM_H_
#define VECTORGRAPH_RENDERER_PLANE_PARAM_H_

#include <Eigen/Core>

#include "types.h"

struct LineSegment {
	bool Intersect(const Point_2& p0, const Point_2& p1, Point_2* p) {
		K a = p1.x() - p0.x();
		K b = -n1;
		K c = x - p0.x();
		K d = p1.y() - p0.y();
		K e = -n2;
		K f = y - p0.y();
		K n = a * e - b * d;
		if (n == 0)
			return false;
		K xx = (c * e - f * b) / n;
		K yy = (-c * d + f * a) / n;

		if (xx < 0 || xx >= 1)
			return false;
		*p = Point_2(xx * a + p0.x(), xx * d + p0.y());

		return true;
	}
	K x, y, z;
	K n1, n2, n3;
};

struct PlaneParam
{
public:
	PlaneParam();
	PlaneParam(
		const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2);

	// make this an inline funcion for fast computation
	/*
	K ComputeDepth(const Point_2& p) const {
		if (invalid_)
			return K(-1e30);
		return (d_ - p.x() * n1_ - p.y() * n2_) / n3_;
	}
	*/
	bool ProjectiveIntersection(const PlaneParam& other, LineSegment* pl) const;

	K ComputeDepth(const Point_2& p) const {
		if (invalid_ )
			return K(-1e30);
		K dx = p.x();
		K dy = p.y();
		K dz = K(1);
		K dis = dx * n1_ + dy * n2_ + dz * n3_;

		K z = -d_ / dis;
		if (z < 0)
			printf("Inverse Depth!\n");
		return z;
	}

	bool InsidePlane(K x, K y, K z) const {
		return (!invalid_) && (n1_ * x + n2_ * y + n3_ * z + d_ == 0);
	}

	bool invalid_;
	K n1_, n2_, n3_, d_;
};

#endif