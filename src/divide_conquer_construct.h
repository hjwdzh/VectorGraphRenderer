#ifndef VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_
#define VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_

#include <vector>

#include <Eigen/Core>

#include "types.h"
#include "plane_param.h"
#include "split_data.h"

void ConstructArrangement(const std::vector<Eigen::Vector3d>& vertices,
	const std::vector<Eigen::Vector3i>& faces,
	const std::vector<PlaneParam>& plane_params,
	int start, int end, Arrangement_2* overlay);

void MergeArrangement(const Arrangement_2& arr1, const Arrangement_2& arr2,
	const std::vector<Eigen::Vector3d>& vertices,
	const std::vector<Eigen::Vector3i>& faces,
	const std::vector<PlaneParam>& plane_params,
	Arrangement_2* pout);

void RelabelFaceFromArrangement(const std::vector<PlaneParam>& plane_params, Arrangement_2* pout);

#endif