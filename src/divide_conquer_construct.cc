#include "divide_conquer_construct.h"

#include <unordered_map>

#include "arrangement.h"
#include "split_data.h"

void ConstructArrangement(const std::vector<Eigen::Vector3d>& vertices,
	const std::vector<Eigen::Vector3i>& faces,
	const std::vector<PlaneParam>& plane_params,
	int start, int end, Arrangement_2* overlay) {
	printf("<%d %d>\n", start, end);
	if (start == end) {
		ComputeTriangleArrangement(vertices, faces[start], start, overlay);
	}
	else {
		//printf("<%d %d>\n", start, end);
		Arrangement_2 overlay1, overlay2;
		int m = (start + end) / 2;
		ConstructArrangement(vertices, faces, plane_params, start, m, &overlay1);
		ConstructArrangement(vertices, faces, plane_params, m + 1, end, &overlay2);
		MergeArrangement(overlay1, overlay2, vertices, faces, plane_params, overlay);
	}
	printf("finish <%d %d>\n", start, end);
}

void MergeArrangement(const Arrangement_2& arr1, const Arrangement_2& arr2,
	const std::vector<Eigen::Vector3d>& vertices,
	const std::vector<Eigen::Vector3i>& faces,
	const std::vector<PlaneParam>& plane_params,
	Arrangement_2* pout) {

	Overlay_traits         overlay_traits;
	overlay (arr1, arr2, *pout, overlay_traits);

	auto& out = *pout;

	auto compute_z_dz = [&](int id1, Arrangement_2::Halfedge_handle u, int source = 1) {
		auto& param = plane_params[id1];
		auto v_src = u->source();
		auto v_tar = u->target();
		auto& p = (source == 1) ? v_src->point() : v_tar->point();
		return param.compute_z(p);
	};

	std::vector<std::vector<SplitData> > splits;
	std::unordered_set<Arrangement_2::Halfedge_handle> halfedges;
	for (auto e = out.halfedges_begin(); e != out.halfedges_end(); ++e) {
		halfedges.insert(e);
	}

	for (auto fit = out.faces_begin(); fit != out.faces_end(); ++fit) {
		if (fit == out.unbounded_face())
			continue;
		long long seg = 1 << 16;
		seg *= seg;
		int id1 = fit->data() / seg - 1;
		int id2 = fit->data() % seg - 1;
		if (id1 != -1 && id2 != -1) {
			std::vector<SplitData> split_points;
			Arrangement_2::Ccb_halfedge_circulator curr = fit->outer_ccb();
			int current_size = split_points.size();
			bool valid = true;

			do {
				Arrangement_2::Halfedge_handle he = curr;
				if (halfedges.count(he) == 0) {
					valid = false;
					break;
				}
				auto z_src1 = compute_z_dz(id1, he, 1);
				auto z_src2 = compute_z_dz(id2, he, 1);
				auto z_tar1 = compute_z_dz(id1, he, 0);
				auto z_tar2 = compute_z_dz(id2, he, 0);

				if ((z_src1 - z_src2) < K(0) && (z_tar1 - z_tar2) > K(0) ||
					(z_src1 - z_src2) > K(0) && (z_tar1 - z_tar2) < K(0)) {
					split_points.push_back(SplitData(he, z_src1 - z_src2, z_tar1 - z_tar2));
				}
			} while (++curr != fit->outer_ccb());
			if (!valid) {
				split_points.resize(current_size);
			} else {
				std::sort(split_points.data() + current_size, split_points.data() + split_points.size());
			}

			for (auto hi = fit->holes_begin(); hi != fit->holes_end(); ++hi) {
				auto circ = *hi;
				Arrangement_2::Halfedge_handle curr = circ;
				int current_size = split_points.size();
				bool valid = true;
				do {
					Arrangement_2::Halfedge_handle he = curr;
					if (halfedges.count(he) == 0) {
						valid = false;
						break;
					}
					auto z_src1 = compute_z_dz(id1, he, 1);
					auto z_src2 = compute_z_dz(id2, he, 1);
					auto z_tar1 = compute_z_dz(id1, he, 0);
					auto z_tar2 = compute_z_dz(id2, he, 0);

					if ((z_src1 - z_src2) < K(0) && (z_tar1 - z_tar2) > K(0) ||
						(z_src1 - z_src2) > K(0) && (z_tar1 - z_tar2) < K(0)) {
							split_points.push_back(SplitData(he, z_src1 - z_src2, z_tar1 - z_tar2));
					}
				} while (++curr != circ);
				if (!valid) {
					split_points.resize(current_size);
				} else {
					std::sort(split_points.data() + current_size, split_points.data() + split_points.size());
				}
			}

			if (split_points.size() != 0) {
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

				splits.push_back(split_points);
			}			
		}
	}

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

	for (auto fit = out.faces_begin(); fit != out.faces_end(); ++fit) {
		if (fit == out.unbounded_face())
			continue;
		long long seg = 1 << 16;
		seg *= seg;
		int id1 = fit->data() / seg - 1;
		int id2 = fit->data() % seg - 1;
		if (id1 == -1) {
			fit->set_data((id2 + 1));
		}
		else if (id2 == -1) {
			fit->set_data((id1 + 1));
		}
		else {
			Arrangement_2::Ccb_halfedge_circulator curr = fit->outer_ccb();
			int selected = -1;
			do {
				Arrangement_2::Halfedge_handle he = curr;
				auto z_src1 = compute_z_dz(id1, he, 1);
				auto z_src2 = compute_z_dz(id2, he, 1);
				z_src1 -= z_src2;
				if (z_src1 < 0) {
					selected = 0;
					break;
				}
				if (z_src1 > 0) {
					selected = 1;
					break;
				}
			} while (++curr != fit->outer_ccb());

			if (selected == 0)
				fit->set_data(id1 + 1);
			else
				fit->set_data(id2 + 1);
		}
	}

	MergeFaceInsideArrangement(out);
}
