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

#include "arrangement.h"
#include "divide_conquer_construct.h"
#include "plane_param.h"

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

	ConstructArrangement(vertices, faces, params, 0, faces.size() - 1, &overlay);
	printf("Done...\n");

	auto compute_z = [&](int id1, double x, double y) {
		Eigen::Vector3d v0 = vertices[faces[id1][0]];
		Eigen::Vector3d v1 = vertices[faces[id1][1]];
		Eigen::Vector3d v2 = vertices[faces[id1][2]];
		Eigen::Vector3d norm = (v1 - v0).cross(v2 - v0);
		if (norm.norm() < 1e-10 || norm[2] == 0) {
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
