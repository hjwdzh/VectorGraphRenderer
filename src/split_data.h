#ifndef VECTORGRAPH_RENDER_SPLIT_DATA_H_
#define VECTORGRAPH_RENDER_SPLIT_DATA_H_

#include "plane_param.h"
#include "types.h"

struct SplitData
{
public:
	SplitData();
	SplitData(Arrangement_2::Halfedge_handle _e, K dis1, K dis2);

	Arrangement_2::Halfedge_handle e;
	Arrangement_2::Vertex_handle v;
	long long fid;
	Point_2 split_point;
	bool trivial;

	// make this an inline funcion for fast computation
	K signature() const {
		return K(0.3245) * split_point.x() + K(0.6842) * split_point.y();
	}

	// make this an inline funcion for fast computation
	bool operator<(const SplitData& n) const {
		return signature() < n.signature();
	}
};

K ComputeDepth(const PlaneParam& plane, Arrangement_2::Halfedge_handle u, int source = 1);

void CollectSelfIntersection(Arrangement_2::Ccb_halfedge_circulator curr,
	int id1, int id2,
	const std::unordered_set<Arrangement_2::Halfedge_handle>& halfedges,
	const std::vector<PlaneParam>& plane_params,
	std::vector<SplitData>* p_split_points);

void SortSplitPoint(std::vector<SplitData>* p_split_points);

void SubdivideAtIntersection(std::vector<std::vector<SplitData> >& splits, Arrangement_2* pout);

#endif