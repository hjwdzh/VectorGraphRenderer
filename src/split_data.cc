#include "split_data.h"

SplitData::SplitData()
{}

/*
SplitData::SplitData(Arrangement_2::Halfedge_handle _e, K dis1, K dis2) {
	e = _e;

	K ratio1 = dis2 / (dis2 - dis1);
	K ratio2 = -dis1 / (dis2 - dis1);

	split_point = Point_2(e->source()->point().x() * ratio1 + e->target()->point().x() * ratio2,
		e->source()->point().y() * ratio1 + e->target()->point().y() * ratio2);
}
*/
SplitData::SplitData(Arrangement_2::Halfedge_handle _e, const Point_2& p) {
	e = _e;
	split_point = p;
}

K ComputeDepth(const PlaneParam& param, Arrangement_2::Halfedge_handle u, int source) {
	auto v_src = u->source();
	auto v_tar = u->target();
	auto& p = (source == 1) ? v_src->point() : v_tar->point();
	return param.ComputeDepth(p);
}

void CollectSelfIntersection(Arrangement_2::Ccb_halfedge_circulator curr,
	int id1, int id2,
	const std::unordered_set<Arrangement_2::Halfedge_handle>& halfedges,
	const std::vector<PlaneParam>& plane_params,
	std::vector<SplitData>* p_split_points) {
	auto& split_points = *p_split_points;
	int current_size = split_points.size();
	bool valid = true;
	auto origin = curr;

	auto f = curr->face();
	LineSegment l;
	
	plane_params[id1].ProjectiveIntersection(plane_params[id2], &l);

	do {
		Arrangement_2::Halfedge_handle he = curr;
		if (halfedges.count(he) == 0) {
			valid = false;
			break;
		}
		//auto pp = he->source()->point();
		//vs.push_back(Eigen::Vector3d(pp.x().convert_to<double>(), pp.y().convert_to<double>(), 0));
		Point_2 p;
		if (l.Intersect(he->source()->point(), he->target()->point(), &p)) {
			split_points.push_back(SplitData(he, p));
		}
	} while (++curr != origin);
	if (!valid) {
		split_points.resize(current_size);
		return;
	}
	//if (valid && split_points.size() % 2 == 1) {
		//printf("HHHH %d %d\n", id1, id2);
		//exit(0);
	//}
	if (!valid || split_points.size() % 2 == 1) {
		split_points.resize(current_size);
	}
}

void SortSplitPoint(std::vector<SplitData>* p_split_points) {
	auto& split_points = *p_split_points;

	std::sort(split_points.begin(), split_points.end());

	for (int j = 0; j < split_points.size(); j += 2) {
		if (split_points[j].e->face() == split_points[j+1].e->face() ||
			split_points[j].e->face() == split_points[j+1].e->twin()->face()) {
			split_points[j].e = split_points[j].e->twin();
		}
		if (split_points[j].e->twin()->face() == split_points[j+1].e->twin()->face()) {
			split_points[j+1].e = split_points[j+1].e->twin();
		}
	}

	for (int j = 0; j < split_points.size(); j += 2)
		split_points[j + 1].fid = split_points[j + 1].e->face()->data();

}

void SubdivideAtIntersection(std::vector<std::vector<SplitData> > & splits, Arrangement_2* pout) {
	auto& out = *pout;
	std::unordered_map<Arrangement_2::Halfedge_handle, std::vector<std::pair<int, int> > > edge_handles;
	for (int i = 0; i < splits.size(); ++i) {
		for (int j = 0; j < splits[i].size(); ++j) {
			auto it = edge_handles.find(splits[i][j].e);
			if (it != edge_handles.end()) {
				it->second.push_back(std::make_pair(i, j));
			}
			else {
				it = edge_handles.find(splits[i][j].e->twin());
				if (it != edge_handles.end()) {
					it->second.push_back(std::make_pair(i, j));
				} else {
					std::vector<std::pair<int, int> > v;
					v.push_back(std::make_pair(i, j));
					edge_handles[splits[i][j].e] = v;
				}
			}
		}
	}


	for (auto& info : edge_handles) {
		auto p = info.first->source()->point();
		std::vector<std::pair<K, int> > dis(info.second.size());
		for (int i = 0; i < info.second.size(); ++i) {
			auto e = info.second[i];
			auto diff_p = splits[e.first][e.second].split_point - p;
			dis[i] = std::make_pair(diff_p.x() * diff_p.x() + diff_p.y() * diff_p.y(), i);
		}

		std::sort(dis.begin(), dis.end());

		std::vector<int> dis_id;
		for (int i = 0; i < dis.size(); ++i) {
			if (i > 0) {
				if (dis[i].first != dis[i - 1].first) {
					dis_id.push_back(i);
				}
			} else {
				dis_id.push_back(0);
			}
		}
		auto e_handle = info.first;
		for (int ii = dis_id.size() - 1; ii >= 0; --ii) {
			int i = dis_id[ii];
			auto e = info.second[dis[i].second];
			auto& split_point = splits[e.first][e.second].split_point;
			
			if (split_point == e_handle->source()->point()) {
				for (int j = 0; j <= ii; ++j) {
					int i = dis_id[j];
					auto e = info.second[dis[i].second];
					splits[e.first][e.second].v = e_handle->source();
				}
				break;
			}
			else if (split_point == e_handle->target()->point()) {
				splits[e.first][e.second].v = e_handle->target();
				continue;
			}
			Segment_2 s1(e_handle->source()->point(), split_point);
			Segment_2 s2(split_point, e_handle->target()->point());

			e_handle = out.split_edge(e_handle, s1, s2);
			splits[e.first][e.second].v = e_handle->target();
		}
		for (int i = 1; i < dis.size(); ++i) {
			if ((dis[i].first == dis[i - 1].first)) {
				auto e1 = info.second[dis[i - 1].second];
				auto e = info.second[dis[i].second];
				splits[e.first][e.second].v = splits[e1.first][e1.second].v;
			}
		}
	}

	for (int i = 0; i < splits.size(); ++i) {
		auto& points = splits[i];
		for (int j = 0; j < points.size(); j += 2) {
			if (points[j].v->point() != points[j + 1].v->point()) {
				Segment_2 s(points[j].v->point(), points[j+1].v->point());
				auto new_e = out.insert_at_vertices(s, points[j].v, points[j + 1].v);		
				new_e->face()->set_data(points[j + 1].fid);
				new_e->twin()->face()->set_data(points[j + 1].fid);
			}
		}
	}	
}