#include "postprocess.h"

#include <fstream>
#include <map>

#include <Eigen/Dense>

#include "simple_svg.h"

using namespace svg;

int EvaluateEdgeType(Arrangement_2::Halfedge_const_handle he, const Mesh& mesh, const Arrangement_2& overlay, double angle_thres) {
	int fid = he->face()->data() - 1;
	auto& params = mesh.GetPlanes();
	if (he->twin()->face() == overlay.unbounded_face()) {
		return 2;
	}
	else {
		int fid1 = he->twin()->face()->data() - 1;
		if (fid1 == -1)
			return 2;
		if (fid1 != fid) {
			K z = params[fid].ComputeDepth(he->source()->point());
			K n_z = params[fid1].ComputeDepth(he->source()->point());
			if (std::abs((z-n_z).convert_to<double>()) > 1e-6) {
				return 2;
			} else {
				K z = params[fid].ComputeDepth(he->target()->point());
				K n_z = params[fid1].ComputeDepth(he->target()->point());
				if (std::abs((z-n_z).convert_to<double>()) > 1e-6)
					return 2;

				auto& p1 = mesh.GetPlanes()[fid];
				auto& p2 = mesh.GetPlanes()[fid1];
				K dot = p1.n1_ * p2.n1_ + p1.n2_ * p2.n2_ + p1.n3_ * p2.n3_;
				double t1 = dot.convert_to<double>();
				double t2 = cos(angle_thres * 3.141592654 / 180.0);

				if (std::abs(t1) < t2) {
					return 1;
				}
			}
		}
	}

	return 0;
}

void PostProcess::CollectFaceAndVertices(const Mesh& mesh, const Arrangement_2& overlay, double angle_thres) {
	auto& points = points_;
	auto& facets_id = facets_id_;
	auto& facets = facets_;
	auto& params = mesh.GetPlanes();
	auto& depths = depths_;
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
				K z = params[fid].ComputeDepth(he->source()->point());
				points.push_back(key);
				depths.push_back(z);
			} else {
				e.outer_indices.push_back(it->second);
			}
			e.outer_type.push_back(EvaluateEdgeType(he, mesh, overlay, angle_thres));
		} while (++curr != e_handle);

		for (auto hi = fit->holes_begin(); hi != fit->holes_end(); ++hi) {
			auto circ = *hi;
			Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
			e.inner_indices.push_back(std::vector<int>());
			e.inner_type.push_back(std::vector<int>());
			do {
				Arrangement_2::Halfedge_const_handle he = curr;
				auto key = std::make_pair(he->source(), fid);
				auto it = vertexID.find(key);
				if (it == vertexID.end()) {
					vertexID[key] = points.size();
					e.inner_indices.back().push_back(points.size());
					K z = params[fid].ComputeDepth(he->source()->point());
					points.push_back(key);
					depths.push_back(z);
	
				} else {
					e.inner_indices.back().push_back(it->second);
				}
				e.inner_type.back().push_back(EvaluateEdgeType(he, mesh, overlay, angle_thres));
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
	auto& depths = depths_;
	auto& degrees = degrees_;

	ComputeDegree();

	for (auto& f : facets) {
		RemoveRedundantVerticesFromLoop(f.outer_indices, f.outer_type);
		for (int i = 0; i < f.inner_indices.size(); ++i) {
			RemoveRedundantVerticesFromLoop(f.inner_indices[i], f.inner_type[i]);
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
			points[top] = points[i];
			depths[top] = depths[i];
			top += 1;
		}
	}

	points.resize(top);
	depths.resize(top);
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

	auto& facets = facets_;
	auto& points = points_;
	auto& depths = depths_;
	// merge duplex
	std::map<std::pair<int,std::pair<int,int> >, int> vID;
	std::vector<int> compressed_vertexID;

	int top = 0;
	compressed_vertexID.resize(points.size());
	std::vector<VertexSignature> points_buf;
	std::vector<K> depths_buf;
	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = depths[i].convert_to<double>();

		auto key = std::make_pair(int(x * 1e5), std::make_pair(int(y * 1e5), int(z * 1e5)));
		auto it = vID.find(key);
		if (it == vID.end()) {
			compressed_vertexID[i] = top;
			points_buf.push_back(points[i]);
			depths_buf.push_back(depths[i]);
			vID[key] = top++;
		} else {
			compressed_vertexID[i] = it->second;
		}
	}
	points = points_buf;
	depths = depths_buf;

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

void PostProcess::SaveToFile(const Mesh& mesh, const char* filename) {
	auto & vertices = mesh.GetVertices();
	auto & faces = mesh.GetFaces();
	auto & facets = facets_;
	auto & face_normals = mesh.GetFaceNormals();
	auto & facets_id = facets_id_;
	auto & edge2face = edge2face_;
	auto & points = points_;
	auto & depths = depths_;

	std::ofstream os;

	os.open(filename);

	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = depths[i].convert_to<double>();

		Eigen::Vector3d v(x * z, y * z, z);
		//Eigen::Vector3d v(x, y, z);
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

	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = depths[i].convert_to<double>();

		Eigen::Vector3d v(x * z, y * z, z);
		//Eigen::Vector3d v(x, y, z);
		os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
	}
	os << "### occlude boundaries\n";
	for (auto& i : facets) {
		if (i.outer_indices.size() < 3)
			continue;
		for (int j = 0; j < i.outer_indices.size(); ++j) {
			if (i.outer_type[j] == 2) {
				int next_j = (j + 1) % i.outer_indices.size();
				os << "l " << i.outer_indices[j] + 1
				   << " " << i.outer_indices[next_j] + 1 << "\n";
			}
		}
	}

	for (int i = 0; i < points.size(); ++i) {
		double x = points[i].first->point().x().convert_to<double>();
		double y = points[i].first->point().y().convert_to<double>();
		double z = depths[i].convert_to<double>();

		Eigen::Vector3d v(x * z, y * z, z);
		//Eigen::Vector3d v(x, y, z);
		os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
	}
	os << "### sharp edges\n";
	for (auto& i : facets) {
		if (i.outer_indices.size() < 3)
			continue;
		for (int j = 0; j < i.outer_indices.size(); ++j) {
			if (i.outer_type[j] == 1) {
				int next_j = (j + 1) % i.outer_indices.size();
				os << "l " << i.outer_indices[j] + 1
				   << " " << i.outer_indices[next_j] + 1 << "\n";
			}
		}
	}
	os.close();	
}

void PostProcess::SaveToSVG(const Mesh& mesh, const char* filename) {
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

	    Polygon poly(Color(128,128,128), Stroke(.0, Color(150, 160, 200)));

	    for (auto& e : i.outer_indices) {
	    	auto& p = points[e].first->point();
	    	int x = (p.x().convert_to<double>()+0.5) * 1000;
	    	int y = (p.y().convert_to<double>()+0.5) * 1000;
	    	poly << Point(x, 1000-y);
	    }
	    doc << poly;
	}

	for (int seq = 1; seq <= 2; ++seq) {
		for (auto& i : facets) {
			if (i.outer_indices.size() < 3)
				continue;
			for (int j = 0; j < i.outer_indices.size(); ++j) {
				if (i.outer_type[j] > 0 && i.outer_type[j] == seq) {
			    	Polyline poly(Stroke(2, Color::Blue));
			    	if (i.outer_type[j] == 2)
			    		poly = Polyline(Stroke(2, Color::Red));
				    {
					    auto& p1 = points[i.outer_indices[j]].first->point();
				    	int x = (p1.x().convert_to<double>()+0.5) * 1000;
				    	int y = (p1.y().convert_to<double>()+0.5) * 1000;
				    	poly << Point(x, 1000-y);
				    }
				    {
				    	int next_j = (j + 1) % i.outer_indices.size();
					    auto& p1 = points[i.outer_indices[next_j]].first->point();
				    	int x = (p1.x().convert_to<double>()+0.5) * 1000;
				    	int y = (p1.y().convert_to<double>()+0.5) * 1000;
				    	poly << Point(x, 1000-y);
				    }
				    doc << poly;
			    }
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

void PostProcess::RemoveRedundantVerticesFromLoop(std::vector<int>& indices, std::vector<int>& types) {
	auto& points = points_;
	auto& degrees = degrees_;

	std::vector<int> mask(indices.size(), 0);
	for (int i = 0; i < indices.size(); ++i) {
		if (degrees[indices[i]] > 2) {
			mask[i] = 1;
			continue;
		}
		if (types[i] != types[(i + indices.size() - 1) % indices.size()]) {
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
			indices[top] = indices[i];
			types[top] = types[i];
			top += 1;
		}
	}
	indices.resize(top);
	types.resize(top);

}