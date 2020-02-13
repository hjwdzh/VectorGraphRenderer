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

			CollectSelfIntersection(fit->outer_ccb(), id1, id2, halfedges, plane_params, &split_points);
			for (auto hi = fit->holes_begin(); hi != fit->holes_end(); ++hi) {
				auto circ = *hi;
				CollectSelfIntersection(circ, id1, id2, halfedges, plane_params, &split_points);
			}

			if (split_points.size() != 0) {
				SortSplitPoint(&split_points);
				splits.push_back(split_points);
			}			
		}
	}

	SubdivideAtIntersection(splits, &out);

	RelabelFaceFromArrangement(plane_params, &out);

	MergeFaceInsideArrangement(out);
}

void RelabelFaceFromArrangement(const std::vector<PlaneParam>& plane_params, Arrangement_2* pout) {
	auto& out = *pout;
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
				auto z_src1 = ComputeDepth(plane_params[id1], he, 1);
				auto z_src2 = ComputeDepth(plane_params[id2], he, 1);
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
}
