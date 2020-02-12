#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <strstream>
#include <set>
#include <unordered_set>
#include <chrono>
struct combine_indices
{
long long operator()( const long long& lhs, const long long& rhs ) const {
	return lhs * ((long long)(1 << 16)) * ((long long)(1 << 16)) + rhs;
}
};

typedef CGAL::Exact_rational									K;
typedef CGAL::Cartesian<K>                   Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::X_monotone_curve_2                            Segment_2;
typedef CGAL::Arr_face_extended_dcel<Traits_2, long long>            Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>                     Arrangement_2;
typedef CGAL::Arr_face_overlay_traits<Arrangement_2,
                                      Arrangement_2,
                                      Arrangement_2,
                                      combine_indices >  Overlay_traits;

#define THRES 1e-10
inline bool equal(double a, double b) {
	return std::abs(a - b) < THRES;
}
inline bool smaller(double a, double b) {
	return a + THRES < b;
}

struct PlaneParam
{
	PlaneParam(){}
	PlaneParam(const Eigen::Vector3d& v0,
		const Eigen::Vector3d& v1,
		const Eigen::Vector3d& v2) {
		Eigen::Vector3d norm = (v1 - v0).cross(v2 - v0);
		if (norm.norm() < THRES || norm[2] == 0) {
			invalid = true;
		}
		norm.normalize();
		d = v0.dot(norm);
		n1 = norm[0];
		n2 = norm[1];
		n3 = norm[2];
		invalid = false;
	}

	K compute_z(const Point_2& p) const {
		if (invalid)
			return K(-1e30);
		return (d - p.x() * n1 - p.y() * n2) / n3;
	}
	bool invalid;
	K n1, n2, n3, d;
};

struct SplitData
{
	SplitData()
	{
		//printf("Interesting....\n");
		//exit(0);
	}
	SplitData(Arrangement_2::Halfedge_handle _e, K dis1, K dis2) {
		e = _e;

		K ratio1 = dis2 / (dis2 - dis1);
		K ratio2 = -dis1 / (dis2 - dis1);

		split_point = Point_2(e->source()->point().x() * ratio1 + e->target()->point().x() * ratio2,
			e->source()->point().y() * ratio1 + e->target()->point().y() * ratio2);
		
		auto p1 = e->source()->point();
		auto p2 = e->target()->point();
		auto split_point = e->source()->point();

		/*
		printf("<%f %f> <%f %f> <%f %f>\n",
			p1.x().convert_to<double>(),p1.y().convert_to<double>(),
			p2.x().convert_to<double>(),p2.y().convert_to<double>(),
			split_point.x().convert_to<double>(),split_point.y().convert_to<double>());
		*/
		if (split_point == e->source()->point() || split_point == e->target()->point())
			trivial = true;
		else
			trivial = false;
		v1 = e->source();
		v2 = e->target();
		v = Arrangement_2::Vertex_handle();
	}
	Arrangement_2::Halfedge_handle e;
	Arrangement_2::Vertex_handle v, v1, v2;
	long long fid;
	Point_2 split_point;
	bool trivial;

	K signature() const {
		return K(0.3245) * split_point.x() + K(0.6842) * split_point.y();
	}
	bool operator<(const SplitData& n) const {
		return signature() < n.signature();
	}
};

void ComputeTriangleOverlay(const std::vector<Eigen::Vector3d>& vertices, const Eigen::Vector3i& f, int fid, Arrangement_2& out) {
	out = Arrangement_2();
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

	//Count(out);
}

void MergeFace(Arrangement_2& out, bool final_step) {

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

void MergeOverlay(Arrangement_2& arr1, Arrangement_2& arr2, std::vector<Eigen::Vector3d>& vertices, std::vector<Eigen::Vector3i>& faces,
	std::vector<PlaneParam>& plane_params, Arrangement_2& out, bool final_step, int start, int end) {
	//printf("Merge Overlay!\n");
	Overlay_traits         overlay_traits;
	overlay (arr1, arr2, out, overlay_traits);

	/*
	auto compute_z_dz = [&](int id1, Arrangement_2::Halfedge_handle u, int source = 1) {
		Eigen::Vector3d v0 = vertices[faces[id1][0]];
		Eigen::Vector3d v1 = vertices[faces[id1][1]];
		Eigen::Vector3d v2 = vertices[faces[id1][2]];

		Eigen::Vector3d norm = (v1 - v0).cross(v2 - v0);
		if (norm.norm() < THRES || norm[2] == 0) {
			return K(-1e30);
		}
		norm.normalize();
		double d = v0.dot(norm);
		auto v_src = u->source();
		auto v_tar = u->target();
		auto p1 = v_src->point();
		auto p2 = v_tar->point();
		auto p = (source == 1) ? p1 : p2;

		K z = (K(d) - p.x() * K(norm[0]) - p.y() * K(norm[1])) / K(norm[2]);

		auto& param = plane_params[id1];
		K z1 = param.compute_z(p);
		if (z != z1) {
			printf("<%f %f>\n", z.convert_to<double>(), z1.convert_to<double>());
			printf("<%f %f>\n", param.n1.convert_to<double>(), norm[0]);
			printf("<%f %f>\n", param.n2.convert_to<double>(), norm[1]);
			printf("<%f %f>\n", param.n3.convert_to<double>(), norm[2]);
			printf("<%f %f>\n", param.d.convert_to<double>(), d);
			exit(0);
		}
		return z;
	};
	*/

	auto compute_z_dz = [&](int id1, Arrangement_2::Halfedge_handle u, int source = 1) {
		auto& param = plane_params[id1];
		auto v_src = u->source();
		auto v_tar = u->target();
		auto& p = (source == 1) ? v_src->point() : v_tar->point();
		return param.compute_z(p);
	};

	auto check_edges = [&]() {
		int invalid_count = 0;
		for (auto e = out.halfedges_begin(); e != out.halfedges_end(); ++e) {
			if (e->face() == out.unbounded_face())
				continue;
			if (e->face() == e->twin()->face()) {
				invalid_count += 1;
			}
		}
		return invalid_count;
	};

	//printf("Merge Overlay 1!\n");
	std::vector<std::vector<SplitData> > splits;
	std::unordered_set<Arrangement_2::Halfedge_handle> halfedges;
	for (auto e = out.halfedges_begin(); e != out.halfedges_end(); ++e) {
		halfedges.insert(e);
	}
	auto t1 = std::chrono::steady_clock::now();
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
					if (split_points[j].e->twin()->face() != split_points[j+1].e->face()) {
						printf("Unable to pair!\n");
						exit(0);
					}
				}

				for (int j = 0; j < split_points.size(); j += 2)
					split_points[j + 1].fid = split_points[j + 1].e->face()->data();

				splits.push_back(split_points);
			}			
			if (split_points.size() % 2 != 0) {
				printf("...... %d\n", split_points.size());
				exit(0);
			}
		}
	}

	/*
	for (int i = 0; i < splits.size(); ++i) {
		int top = 0;
		for (int j = 0; j < splits[i].size(); j += 2) {
			if (!splits[i][j].trivial || !splits[i][j + 1].trivial) {
				splits[i][top++] = splits[i][j];
				splits[i][top++] = splits[i][j + 1];
			}
		}
		splits[i].resize(top);
	}
	*/
	auto t2 = std::chrono::steady_clock::now();

	std::map<Arrangement_2::Halfedge_handle, std::vector<std::pair<int, int> > > edge_handles;
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

	auto t3 = std::chrono::steady_clock::now();


	//printf("================================\n");
	//std::ofstream os;
	//os.open("split_point.obj");
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

			//os << "v " << split_point.x().convert_to<double>() << " " << split_point.y().convert_to<double>() << " 0\n";
			e_handle = out.split_edge(e_handle, s1, s2);
			splits[e.first][e.second].v = e_handle->target();
			/*
			if (vid.count(e_handle->target()) == 0) {
				int n = vid.size();
				vid[e_handle->target()] = n;
				printf("New splits %d %d\n", vid.size(), i);
			}
			*/
		}
		for (int i = 1; i < dis.size(); ++i) {
			if ((dis[i].first == dis[i - 1].first)) {
				auto e1 = info.second[dis[i - 1].second];
				auto e = info.second[dis[i].second];
				//printf("Skip %d\n", i);
				splits[e.first][e.second].v = splits[e1.first][e1.second].v;
			}
		}
	}

	auto t4 = std::chrono::steady_clock::now();
	/*
	os.close();

	int invalid_count = check_edges();
	if (invalid_count > 0) {
		exit(0);
	}

	for (int i = 0; i < splits.size(); ++i) {
		for (int j = 0; j < splits[i].size(); ++j) {
			auto& p = splits[i][j];
			printf("<%d %d %d>\n", vid[splits[i][j].v1], vid[splits[i][j].v2], vid[splits[i][j].v]);
		}
	}
	//Count(out);
	//printf("Merge Overlay 2!\n");
	os.open("split_new_edge.obj");
	*/
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
	auto t5 = std::chrono::steady_clock::now();
	/*
	for (int i = 0; i < count; ++i) {
		os << "l " << i * 2 + 1 << " " << i * 2 + 2 << "\n";
	}
	os.close();
	*/
	//Count(out);
	//printf("Merge Overlay 3!\n");
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
			K min_k = K(100), max_k = K(-100);
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

	auto t6 = std::chrono::steady_clock::now();
	//printf("Merge Overlay 4!\n");
	//printf("OK1!\n");

	//Count(out);
	MergeFace(out, (start == 0) && (end == 1000));
	auto t7 = std::chrono::steady_clock::now();

	//printf("OK2!\n");
	//printf("Merge Overlay Done!\n");
	//Count(out);

	auto diff1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
	auto diff2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count();
	auto diff3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3).count();
	auto diff4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t4).count();
	auto diff5 = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5).count();
	auto diff6 = std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count();
	std::cout << diff1 << " " << diff2 << " " << diff3 << " " << diff4 << " " << diff5 << " " << diff6 << "\n";
}

void ComputeOverlay(std::vector<Eigen::Vector3d>& vertices,
	std::vector<Eigen::Vector3i>& faces,
	std::vector<PlaneParam>& plane_params,
	int start, int end, Arrangement_2& overlay) {
	printf("<%d %d>\n", start, end);
	if (start == end) {
		ComputeTriangleOverlay(vertices, faces[start], start, overlay);
	}
	else {
		//printf("<%d %d>\n", start, end);
		Arrangement_2 overlay1, overlay2;
		int m = (start + end) / 2;
		ComputeOverlay(vertices, faces, plane_params, start, m, overlay1);
		ComputeOverlay(vertices, faces, plane_params, m + 1, end, overlay2);
		MergeOverlay(overlay1, overlay2, vertices, faces, plane_params, overlay, false, start, end);//(start == 0 && end == faces.size() - 1));
	}
	printf("finish <%d %d>\n", start, end);
}

void Clamp(std::vector<Eigen::Vector3d>& vertices, std::vector<Eigen::Vector3i>& faces, std::vector<Eigen::Vector3d>& face_normals, int dim, double thres, int comp) {
	int i = 0;
	int faces_num = faces.size();
	while (i < faces_num) {
		Eigen::Vector3d v[3];
		for (int j = 0; j < 3; ++j)
			v[j] = vertices[faces[i][j]];
		int valids[3];
		for (int j = 0; j < 3; ++j) {
			if (comp == 0 && v[j][dim] < thres + -1e-6 || comp == 1 && v[j][dim] > thres + 1e-6) {
				valids[j] = 0;
			} else {
				valids[j] = 1;
			}
		}
		int sum = valids[0] + valids[1] + valids[2];
		if (sum == 0) {
			faces[i] = faces.back();
			face_normals[i] = face_normals.back();
			faces.pop_back();
			face_normals.pop_back();
			faces_num = faces.size();
			continue;
		}
		else if (sum == 3) {
		}
		else if (sum == 1) {
			int j = 0;
			while (valids[j] == 0)
				j += 1;
			double step0 = 0, step1 = 0;
			if (comp == 0) {
				step0 = (thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			else {
				step0 = (thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			int nx = vertices.size();
			int nx2 = vertices.size() + 1;
			vertices.push_back(v[j] + step0 * (v[(j+2) % 3] - v[j]));
			vertices.push_back(v[j] + step1 * (v[(j+1) % 3] - v[j]));
			faces[i] = Eigen::Vector3i(nx, faces[i][j], nx2);
		}
		else {
			int j = 0;
			while (valids[j] == 1)
				j += 1;
			double step0 = 0, step1 = 0;
			if (comp == 0) {
				step0 = (thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			else {
				step0 = (thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			int v0 = faces[i][(j + 2) % 3];
			int nx = vertices.size();
			int nx2 = vertices.size() + 1;
			int v1 = faces[i][(j + 1) % 3];
			vertices.push_back(v[j] + step0 * (v[(j+2) % 3] - v[j]));
			vertices.push_back(v[j] + step1 * (v[(j+1) % 3] - v[j]));
			faces[i] = Eigen::Vector3i(v0, nx, nx2);
			faces.push_back(Eigen::Vector3i(v0, nx2, v1));
			face_normals.push_back(face_normals[i]);
		}
		i += 1;
	}	
}

int main (int argc, char** argv)
{
	std::ifstream is(argv[1]);
	char buffer[1024];
	std::vector<Eigen::Vector3d> vertices;
	std::vector<Eigen::Vector3d> face_normals;
	std::vector<Eigen::Vector3i> faces;

	while (is >> buffer) {
		if (strcmp(buffer, "v") == 0) {
			Eigen::Vector3d v;
			is >> v[0] >> v[1] >> v[2];
			vertices.push_back(v);
		}
		else if (strcmp(buffer, "f") == 0) {
			Eigen::Vector3i f;
			for (int i = 0; i < 3; ++i) {
				is >> buffer;
				int t = 0;
				int p = 0;
				int l = strlen(buffer);
				while (p != l && buffer[p] >= '0' && buffer[p] <= '9') {
					t = t * 10 + (buffer[p] - '0');
					p += 1;
				}
				f[i] = t - 1;
			}
			faces.push_back(f);
		}
	}
	is.close();

	is.open(argv[2]);

	Eigen::Matrix4d world2cam;

	int width, height;
	double fx, fy, cx, cy;
	double angle;
	is >> height >> width;
	is >> fx >> fy >> cx >> cy;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			is >> world2cam(i, j);
		}
	}
	is >> angle;


	for (int i = 0; i < vertices.size(); ++i) {
		Eigen::Vector4d v(vertices[i][0], vertices[i][1], vertices[i][2], 1);
		v = world2cam * v;
		vertices[i] = Eigen::Vector3d(v[0],v[1],v[2]);
	}

	face_normals.resize(faces.size());
	for (int i = 0; i < faces.size(); ++i) {
		Eigen::Vector3d d1 = vertices[faces[i][1]] - vertices[faces[i][0]];
		Eigen::Vector3d d2 = vertices[faces[i][2]] - vertices[faces[i][0]];
		Eigen::Vector3d n = d1.cross(d2);
		if (n.norm() > 0)
			n = n.normalized();
		face_normals[i] = n;
	}

	Clamp(vertices, faces, face_normals, 2, 1e-2, 0);
	std::vector<Eigen::Vector3d> vertices0 = vertices;

	for (int i = 0; i < vertices.size(); ++i) {
		Eigen::Vector3d v = vertices[i];
		if (v[2] != 0) {
			vertices[i] = Eigen::Vector3d((v[0]/v[2]*fx+cx)/width,
				(v[1]/v[2]*fy+cy)/height,v[2]);
		} else {
			vertices[i] = Eigen::Vector3d(0, 0, 0);
		}
	}

	Clamp(vertices, faces, face_normals, 0, 0, 0);
	Clamp(vertices, faces, face_normals, 1, 0, 0);
	Clamp(vertices, faces, face_normals, 0, 1, 1);
	Clamp(vertices, faces, face_normals, 1, 1, 1);

	std::vector<int> source_face_indices(faces.size());
	for (int i = 0; i < faces.size(); ++i)
		source_face_indices[i] = i;
	//SelfIntersection(vertices, faces, source_face_indices);

	face_normals.resize(faces.size());
	for (int i = 0; i < faces.size(); ++i) {
		Eigen::Vector3d d1 = vertices[faces[i][1]] - vertices[faces[i][0]];
		Eigen::Vector3d d2 = vertices[faces[i][2]] - vertices[faces[i][0]];
		Eigen::Vector3d n = d1.cross(d2);
		if (n.norm() > 0)
			n = n.normalized();
		face_normals[i] = n;
	}

	Arrangement_2 overlay;

	printf("Overlay...\n");
	std::vector<PlaneParam> params(faces.size());
	for (int i = 0; i < faces.size(); ++i) {
		params[i] = PlaneParam(vertices[faces[i][0]],vertices[faces[i][1]],vertices[faces[i][2]]);
	}

	ComputeOverlay(vertices, faces, params, 0, faces.size() - 1, overlay);
	printf("Done...\n");

	auto compute_z = [&](int id1, double x, double y) {
		Eigen::Vector3d v0 = vertices[faces[id1][0]];
		Eigen::Vector3d v1 = vertices[faces[id1][1]];
		Eigen::Vector3d v2 = vertices[faces[id1][2]];
		Eigen::Vector3d norm = (v1 - v0).cross(v2 - v0);
		if (norm.norm() < THRES || norm[2] == 0) {
			return (v0.z() + v1.z() + v2.z()) / 3.0;
		}
		norm.normalize();
		double d = v0.dot(norm);
		//x * norm.x + y * norm.y + z * norm.z = d
		double z = (d - x * norm[0] - y * norm[1]) / norm[2];
		return z;
	};
	
	int detected_deletion = 0;
	for (auto v = overlay.vertices_begin(); v != overlay.vertices_end(); ++v) {
		auto e1 = v->incident_halfedges();
		auto e2 = e1;
		int c = 0;
		do {
			c += 1;
		} while (e2 != e1);
		if (c == 2) {
			detected_deletion += 1;
		}
	}

	struct EdgeList {
		std::vector<int> outer_indices;
		std::vector<std::vector<int> > inner_indices;
	};
	std::vector<EdgeList> facets;
	std::map<std::pair<Arrangement_2::Vertex_const_handle,int>, int> vertexID;


	std::vector<int> facets_id;
	std::vector<std::pair<Arrangement_2::Vertex_const_handle,int> > points;
	for (auto fit = overlay.faces_begin(); fit != overlay.faces_end(); ++fit) {
		int fid = fit->data() - 1;
		if (fit == overlay.unbounded_face() || fid < 0)
			continue;

		EdgeList e;
		auto e_handle = fit->outer_ccb();
		Arrangement_2::Ccb_halfedge_const_circulator curr = e_handle;
		do {
			Arrangement_2::Halfedge_const_handle he = curr;
			auto key = std::make_pair(he->source(), fid);
			auto it = vertexID.find(key);
			if (it == vertexID.end()) {
				vertexID[key] = points.size();
				e.outer_indices.push_back(points.size());

				points.push_back(key);
			} else {
				e.outer_indices.push_back(it->second);
			}
		} while (++curr != e_handle);

		for (auto hi = fit->holes_begin(); hi != fit->holes_end(); ++hi) {
			auto circ = *hi;
			Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
			e.inner_indices.push_back(std::vector<int>());
			do {
				Arrangement_2::Halfedge_const_handle he = curr;
				auto key = std::make_pair(he->source(), fid);
				auto it = vertexID.find(key);
				if (it == vertexID.end()) {
					vertexID[key] = points.size();
					e.inner_indices.back().push_back(points.size());
					points.push_back(key);
				} else {
					e.inner_indices.back().push_back(it->second);
				}
			} while (++curr != circ);
		}
		facets.push_back(e);
		facets_id.push_back(fid);
	}

	// remove 2-degree vertex on edge
	std::vector<int> degrees(points.size());

	auto compute_degree = [&]() {
		for (auto& d :degrees)
			d = 0;
		for (auto& f : facets) {
			for (auto& e : f.outer_indices) {
				degrees[e] += 2;
			}
			for (auto& es : f.inner_indices) {
				for (auto& e : es) {
					degrees[e] += 2;
				} 
			}
		}
	};
	compute_degree();

	auto remove_middle = [&](std::vector<int>& indices) {
		std::vector<int> mask(indices.size(), 0);
		for (int i = 0; i < indices.size(); ++i) {
			if (degrees[indices[i]] > 2) {
				mask[i] = 1;
				continue;
			}
			int prev_id = indices[(i + indices.size() - 1) % indices.size()];
			int next_id = indices[(i + 1) % indices.size()];
			typedef std::pair<Arrangement_2::Vertex_const_handle,int> key;
			const key& p1 = points[prev_id];
			const key& p2 = points[indices[i]];
			const key& p3 = points[next_id];

			auto diff1 = p2.first->point() - p1.first->point();
			auto diff2 = p3.first->point() - p2.first->point();

			K a = diff1.x() * diff2.y() - diff1.y() * diff2.x();
			if (a != K(0)) {
				mask[i] = 1;
			}
		}
		int top = 0;
		for (int i = 0; i < indices.size(); ++i) {
			if (mask[i] == 1) {
				indices[top++] = indices[i];
			}
		}
		indices.resize(top);

	};

	for (auto& f : facets) {
		remove_middle(f.outer_indices);
		for (auto& es : f.inner_indices) {
			remove_middle(es);
		}
	}

	compute_degree();
	std::vector<int> compressed_vertexID(degrees.size());
	compressed_vertexID[0] = 0;
	for (int i = 1; i < compressed_vertexID.size(); ++i) {
		compressed_vertexID[i] = compressed_vertexID[i - 1] + (degrees[i - 1] > 0 ? 1 : 0);
	}

	int top = 0;
	for (int i = 0; i < points.size(); ++i) {
		if (degrees[i]) {
			points[top++] = points[i];
		}
	}

	points.resize(top);
	for (auto& i : facets) {
		for (auto& e : i.outer_indices) {
			e = compressed_vertexID[e];
		}
		for (auto& es : i.inner_indices) {
			for (auto& e : es)
				e = compressed_vertexID[e];
		}
	}

	/*
	// merge duplex
	std::map<std::pair<int,std::pair<int,int> >, int> vID;
	top = 0;
	compressed_vertexID.resize(points.size());
	std::vector<Eigen::Vector3d> points_buf;
	for (int i = 0; i < points.size(); ++i) {
		auto key = std::make_pair(int(points[i][0] * 1e5), std::make_pair(int(points[i][1] * 1e5), int(points[i][2] * 1e5)));
		auto it = vID.find(key);
		if (it == vID.end()) {
			compressed_vertexID[i] = top;
			points_buf.push_back(points[i]);
			vID[key] = top++;
		} else {
			compressed_vertexID[i] = it->second;
		}
	}
	points = points_buf;

	auto shrink = [&](std::vector<int>& v) {
		std::vector<int> m(v.size(), 1);
		for (int i = 0; i < v.size(); ++i) {
			int curr = v[i];
			int next = v[(i + 1) % v.size()];
			if (curr == next) {
				m[i] = 0;
			}
		}
		int top = 0;
		for (int i = 0; i < v.size(); ++i) {
			if (m[i] == 1)
				v[top++] = v[i];
		}
		v.resize(top);
	};

	for (auto& i : facets) {
		for (auto& e : i.outer_indices) {
			e = compressed_vertexID[e];
		}
		shrink(i.outer_indices);
		for (auto& es : i.inner_indices) {
			for (auto& e : es) {
				e = compressed_vertexID[e];
			}
			shrink(es);
		}
	}

	*/
	//recognize faces
	std::map<std::pair<int,int>, std::set<int> > edge2face;
	for (int i = 0; i < facets.size(); ++i) {
		for (int j = 0; j < facets[i].outer_indices.size(); ++j) {
			int v0 = facets[i].outer_indices[j];
			int v1 = facets[i].outer_indices[(j + 1) % facets[i].outer_indices.size()];
			auto key = (v0 < v1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
			edge2face[key].insert(facets_id[i]);
		}
		for (auto& es : facets[i].inner_indices) {
			for (int j = 0; j < es.size(); ++j) {
				int v0 = es[j];
				int v1 = es[(j + 1) % es.size()];
				auto key = (v0 < v1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
				edge2face[key].insert(facets_id[i]);
			}			
		}
	}
	std::ofstream os;

	os.open(argv[3]);

	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = compute_z(points[i].second, x, y);

		//x = (x * width - cx) / fx * z;
		//y = (y * height - cy) / fy * z;

		Eigen::Vector4d v(x - 0.5, y - 0.5, z, 1);
		//v = world2cam.inverse() * v;
		os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
	}

	for (auto& i : facets) {
		os << "f";
		for (auto& e : i.outer_indices) {
			os << " " << e + 1;
		}
		os << "\n";
		for (auto& es : i.inner_indices) {
			os << "###holes### f";
			for (auto& e : es) {
				os << " " << e + 1;
			}
			os << "\n";
		}
	}

	os << "### occlusion boundaries\n";
	for (auto& p : edge2face) {
		if (p.second.size() < 2) {
			os << "l " << p.first.first + 1  << " " << p.first.second + 1 << "\n";
		}
	}

	os << "### sharp edges\n";
	for (auto& p : edge2face) {
		if (p.second.size() >= 2) {
			auto it = p.second.begin();
			int f0 = *it;
			it++;
			int f1 = *it;
			double t = std::abs(face_normals[f0].dot(face_normals[f1]));
			if (t < cos(angle / 180.0 * 3.141592654)) {
				os << "l " << p.first.first + 1 << " " << p.first.second + 1 << "\n";
			}

		}
	}
	os.close();
}
