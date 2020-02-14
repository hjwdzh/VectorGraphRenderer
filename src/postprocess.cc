#include "postprocess.h"

#include <fstream>
#include <map>

#include <Eigen/Dense>

#include "simple_svg.h"

using namespace svg;


void PostProcess::CollectFaceAndVertices(const Arrangement_2& overlay) {
	auto& points = points_;
	auto& facets_id = facets_id_;
	auto& facets = facets_;

	std::map<std::pair<Arrangement_2::Vertex_const_handle,int>, int> vertexID;

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
}

void PostProcess::RemoveRedundantVertices()
{
	auto& facets = facets_;
	auto& points = points_;
	auto& degrees = degrees_;

	ComputeDegree();

	for (auto& f : facets) {
		RemoveRedundantVerticesFromLoop(f.outer_indices);
		for (auto& es : f.inner_indices) {
			RemoveRedundantVerticesFromLoop(es);
		}
	}

	ComputeDegree();
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
}

void PostProcess::MergeDuplex(const Mesh& mesh)
{
	auto & vertices = mesh.GetVertices();
	auto & faces = mesh.GetFaces();
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

	auto& facets = facets_;
	auto& points = points_;
	// merge duplex
	std::map<std::pair<int,std::pair<int,int> >, int> vID;
	std::vector<int> compressed_vertexID;

	int top = 0;
	compressed_vertexID.resize(points.size());
	std::vector<VertexSignature> points_buf;
	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = compute_z(points[i].second, x, y);

		auto key = std::make_pair(int(x * 1e5), std::make_pair(int(y * 1e5), int(z * 1e5)));
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
}


void PostProcess::CollectEdges(const Mesh& mesh) {
	auto & vertices = mesh.GetVertices();
	auto & faces = mesh.GetFaces();
	auto & face_normals = mesh.GetFaceNormals();
	auto & facets_id = facets_id_;
	auto & edge2face = edge2face_;

	auto& facets = facets_;
	auto& points = points_;

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
}

void PostProcess::SaveToFile(const Mesh& mesh, const Camera& camera, const char* filename) {
	auto & vertices = mesh.GetVertices();
	auto & faces = mesh.GetFaces();
	auto & facets = facets_;
	auto & face_normals = mesh.GetFaceNormals();
	auto & facets_id = facets_id_;
	auto & edge2face = edge2face_;
	auto & points = points_;

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

	std::ofstream os;

	os.open(filename);

	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = compute_z(points[i].second, x, y);

		Eigen::Vector4d v(x - 0.5, y - 0.5, z, 1);
		os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
	}

	for (auto& i : facets) {
		if (i.outer_indices.size() < 3)
			continue;
		os << "f";
		for (auto& e : i.outer_indices) {
			os << " " << e + 1;
		}
		os << "\n";
		for (auto& es : i.inner_indices) {
			if (es.size() < 3)
				continue;
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
			if (t < cos(angle_thres / 180.0 * 3.141592654)) {
				os << "l " << p.first.first + 1 << " " << p.first.second + 1 << "\n";
			}

		}
	}
	os.close();	
}

void PostProcess::SaveToSVG(const Mesh& mesh, double angle_thres, const char* filename) {
    auto & facets = facets_;
    auto & points = points_;
    auto & face_normals = mesh.GetFaceNormals();
    auto & edge2face = edge2face_;
   	
   	Dimensions dimensions(1000, 1000);
    Document doc(filename, Layout(dimensions, Layout::BottomLeft));

    // Red image border.
    Polygon border(Stroke(1, Color::Red));
    border << Point(0, 0) << Point(dimensions.width, 0)
        << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    doc << border;
    // Render polygon

	for (auto& i : facets) {
		if (i.outer_indices.size() < 3)
			continue;

	    Polygon poly(Color(rand() % 256, rand() % 256, rand() % 256), Stroke(.0, Color(150, 160, 200)));

	    for (auto& e : i.outer_indices) {
	    	auto& p = points[e].first->point();
	    	int x = p.x().convert_to<double>() * 1000;
	    	int y = p.y().convert_to<double>() * 1000;
	    	poly << Point(x, 1000-y);
	    }
	    doc << poly;

		for (auto& es : i.inner_indices) {
			if (es.size() < 3)
				continue;
		    Polygon poly(Color(255,255,255), Stroke(.0, Color(150, 160, 200)));
			for (auto& e : es) {
				auto& p = points[e].first->point();
		    	int x = p.x().convert_to<double>() * 1000;
		    	int y = p.y().convert_to<double>() * 1000;
		    	poly << Point(x, 1000-y);
			}
			doc << poly;
		}
	}

	for (auto& p : edge2face) {
		if (p.second.size() < 2) {
		    Polyline poly(Stroke(.5, Color::Blue));
		    {
			    auto& p1 = points[p.first.first].first->point();
		    	int x = p1.x().convert_to<double>() * 1000;
		    	int y = p1.y().convert_to<double>() * 1000;
		    	poly << Point(x, 1000-y);
		    }
		    {
			    auto& p1 = points[p.first.second].first->point();
		    	int x = p1.x().convert_to<double>() * 1000;
		    	int y = p1.y().convert_to<double>() * 1000;
		    	poly << Point(x, 1000-y);
		    }
		    doc << poly;
		}
	}

	for (auto& p : edge2face) {
		if (p.second.size() >= 2) {
			auto it = p.second.begin();
			int f0 = *it;
			it++;
			int f1 = *it;
			double t = std::abs(face_normals[f0].dot(face_normals[f1]));
			if (t < cos(angle_thres / 180.0 * 3.141592654)) {
			    Polyline poly(Stroke(.5, Color::Green));
			    {
				    auto& p1 = points[p.first.first].first->point();
			    	int x = p1.x().convert_to<double>() * 1000;
			    	int y = p1.y().convert_to<double>() * 1000;
			    	poly << Point(x, 1000-y);
			    }
			    {
				    auto& p1 = points[p.first.second].first->point();
			    	int x = p1.x().convert_to<double>() * 1000;
			    	int y = p1.y().convert_to<double>() * 1000;
			    	poly << Point(x, 1000-y);
			    }
			    doc << poly;
			}

		}
	}


    doc.save();
}

void PostProcess::ComputeDegree()
{
	degrees_.resize(points_.size());

	for (auto& d :degrees_)
		d = 0;
	for (auto& f : facets_) {
		for (auto& e : f.outer_indices) {
			degrees_[e] += 2;
		}
		for (auto& es : f.inner_indices) {
			for (auto& e : es) {
				degrees_[e] += 2;
			} 
		}
	}
}

void PostProcess::RemoveRedundantVerticesFromLoop(std::vector<int>& indices) {
	auto& points = points_;
	auto& degrees = degrees_;

	std::vector<int> mask(indices.size(), 0);
	for (int i = 0; i < indices.size(); ++i) {
		if (degrees[indices[i]] > 2) {
			mask[i] = 1;
			continue;
		}
		int prev_id = indices[(i + indices.size() - 1) % indices.size()];
		int next_id = indices[(i + 1) % indices.size()];

		auto& p1 = points[prev_id];
		auto& p2 = points[indices[i]];
		auto& p3 = points[next_id];

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

}