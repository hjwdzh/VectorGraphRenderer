#ifndef VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_
#define VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_

#include <vector>

#include <Eigen/Core>

#include "types.h"
#include "plane_param.h"

void MergeArrangement(const Arrangement_2& arr1, const Arrangement_2& arr2,
	const std::vector<Eigen::Vector3d>& vertices,
	const std::vector<Eigen::Vector3i>& faces,
	const std::vector<PlaneParam>& plane_params,
	Arrangement_2* pout);

#endif