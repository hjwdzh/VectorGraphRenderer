#ifndef VECTORGRAPH_RENDERER_PLANE_PARAM_H_
#define VECTORGRAPH_RENDERER_PLANE_PARAM_H_

#include <Eigen/Core>

#include "types.h"

class PlaneParam
{
public:
	PlaneParam();
	PlaneParam(
		const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2);

	// make this an inline funcion for fast computation
	K compute_z(const Point_2& p) const {
		if (invalid_)
			return K(-1e30);
		return (d_ - p.x() * n1_ - p.y() * n2_) / n3_;
	}

private:
	bool invalid_;
	K n1_, n2_, n3_, d_;
};

#endif