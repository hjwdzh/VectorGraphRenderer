#include "split_data.h"

SplitData::SplitData()
{}

SplitData::SplitData(Arrangement_2::Halfedge_handle _e, K dis1, K dis2) {
	e = _e;

	K ratio1 = dis2 / (dis2 - dis1);
	K ratio2 = -dis1 / (dis2 - dis1);

	split_point = Point_2(e->source()->point().x() * ratio1 + e->target()->point().x() * ratio2,
		e->source()->point().y() * ratio1 + e->target()->point().y() * ratio2);
}