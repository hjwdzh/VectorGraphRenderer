#ifndef VECTORGRAPH_RENDER_SPLIT_DATA_H_
#define VECTORGRAPH_RENDER_SPLIT_DATA_H_

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

#endif