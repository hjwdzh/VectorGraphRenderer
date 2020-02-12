#include "arrangement.h"

#include <unordered_set>

void ComputeTriangleArrangement(const std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3i& f, int fid, Arrangement_2* pout) {
	*pout = Arrangement_2();
	auto& out = *pout;
	
	Eigen::Vector3d v1 = vertices[f[0]];
	Eigen::Vector3d v2 = vertices[f[1]];
	Eigen::Vector3d v3 = vertices[f[2]];

	auto off1 = v2 - v1;
	auto off2 = v3 - v1;
	double t = off1.x() * off2.y() - off1.y() * off2.x();
	if (std::abs(t) < 1e-15) {
		for (auto fit = out.faces_begin(); fit != out.faces_end(); ++fit)
			fit->set_data(0);
		return;
	}

	for (int i = 0; i < 3; ++i) {
		auto v1 = vertices[f[i]];
		auto v2 = vertices[f[(i + 1) % 3]];
		Segment_2      s1 (Point_2(v1.x(), v1.y()), Point_2(v2.x(), v2.y()));	
		insert_non_intersecting_curve (out, s1);
	}

	for (auto fit = out.faces_begin(); fit != out.faces_end(); ++fit)
		fit->set_data ((fid + 1) * (fit != out.unbounded_face()));
}


void MergeFaceInsideArrangement(Arrangement_2& out) {

	// Detect half edges with same original face ID and Remove them
	std::unordered_set<Arrangement_2::Halfedge_handle> halfedges;
	for (auto e = out.halfedges_begin(); e != out.halfedges_end(); ++e) {
		if (halfedges.count(e) || halfedges.count(e->twin()))
			continue;
		if (e->face() == out.unbounded_face() || e->twin()->face() == out.unbounded_face())
			continue;
		if (e->face()->data() == e->twin()->face()->data()) {
			halfedges.insert(e);
		}
	}

	for (auto e : halfedges) {
		out.remove_edge(e);
	}

	// Remove isolated or redundant vertices
	for (auto it = out.vertices_begin(); it != out.vertices_end(); ++it) {
		if (it->is_isolated()) {
			remove_vertex(out, it);	
		}
		else {
			auto p1 = it->point();

			auto e1 = it->incident_halfedges();
			auto e2 = e1;
			int incidents = 0;
			do {
				incidents += 1;
				if (incidents > 2)
					break;
				e2++;
			} while (e1 != e2);
			if (incidents == 2) {
				e2++;
				auto& p1 = it->point();
				auto& p2 = e1->source()->point();
				auto& p3 = e2->source()->point();

				auto diff1 = p2 - p1;
				auto diff2 = p3 - p2;
				if (diff1.x() * diff2.y() - diff2.x() * diff1.y() == 0) {
					remove_vertex(out, it);
				}
			}
		}
	}
}