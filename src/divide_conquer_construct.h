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

K ComputeDepth(const PlaneParam& plane, Arrangement_2::Halfedge_handle u, int source = 1);

void CollectSelfIntersection(Arrangement_2::Ccb_halfedge_circulator curr,
	int id1, int id2,
	const std::unordered_set<Arrangement_2::Halfedge_handle>& halfedges,
	const std::vector<PlaneParam>& plane_params,
	std::vector<SplitData>* p_split_points);

void SortSplitPoint(std::vector<SplitData>* p_split_points);

#endif