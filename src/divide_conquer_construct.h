#ifndef VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_
#define VECTORGRAPH_RENDER_DIVIDE_CONQUER_CONSTRUCT_H_

#include <vector>

#include <Eigen/Core>

#include "types.h"
#include "mesh.h"
#include "plane_param.h"
#include "split_data.h"

void ConstructArrangement(const Mesh& mesh,
	int start, int end, Arrangement_2* overlay);

void MergeArrangement(const Arrangement_2& arr1, const Arrangement_2& arr2,
	const Mesh& mesh,
	Arrangement_2* pout);

void RelabelFaceFromArrangement(const std::vector<PlaneParam>& plane_params, Arrangement_2* pout);

#endif