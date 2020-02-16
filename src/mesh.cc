#include "mesh.h"

#include <fstream>

#include <Eigen/Dense>

Mesh::Mesh()
{}

void Mesh::SaveOBJ(const char* filename, const Camera& camera, bool readjust) {
	auto& vertices = vertices_;
	auto& faces = faces_;
	std::ofstream os(filename);
	for (int i = 0; i < vertices.size(); ++i) {
		auto v = vertices[i];
		if (readjust)
			camera.Undistort(v);
		os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
	}
	for (int i = 0; i < faces.size(); ++i) {
		auto& f = faces[i];
		os << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << "\n";
	}
	os.close();
}

Eigen::Vector3d PerspectiveInterp(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double step, bool perspective) {
	if (!perspective)
		return a * (1 - step) + b * step;
	Eigen::Vector3d v0(a[0] * a[2],a[1] * a[2],a[2]);
	Eigen::Vector3d v1(b[0] * b[2],b[1] * b[2],b[2]);
	Eigen::Vector3d n1 = (v1 - v0);
	if (n1.norm() == 0)
		return a;
	n1 /= n1.norm();

	Eigen::Vector3d dir1(a[0] * (1 - step) + b[0] * step, a[1] * (1 - step) + b[1] * step, 1);
	
	Eigen::Vector3d dir = dir1 - dir1.dot(n1) * n1;
	
	double d_norm = dir.squaredNorm();
	double z = v0.dot(dir) / d_norm;

	dir1[2] = z;
	return dir1;
}

void Mesh::LoadFromFile(const char* filename) {
	auto& vertices = vertices_;
	auto& faces = faces_;
	vertices.clear();
	faces.clear();
	char buffer[1024];	
	std::ifstream is(filename);
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
}

void Mesh::BoundaryClip(int dim, double clamp_thres, int comp, bool perspective) {
	auto& vertices = vertices_;
	auto& faces = faces_;

	int i = 0;
	int faces_num = faces.size();
	while (i < faces_num) {
		Eigen::Vector3d v[3];
		for (int j = 0; j < 3; ++j)
			v[j] = vertices[faces[i][j]];
		int valids[3];
		for (int j = 0; j < 3; ++j) {
			if ((comp == 0 && v[j][dim] < clamp_thres + -1e-6)
			 || (comp == 1 && v[j][dim] > clamp_thres + 1e-6)) {
				valids[j] = 0;
			} else {
				valids[j] = 1;
			}
		}
		int sum = valids[0] + valids[1] + valids[2];
		if (sum == 0) {
			faces[i] = faces.back();
			faces.pop_back();
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
				step0 = (clamp_thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (clamp_thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			else {
				step0 = (clamp_thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (clamp_thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			int nx = vertices.size();
			int nx2 = vertices.size() + 1;

			vertices.push_back(PerspectiveInterp(v[j], v[(j+2)%3], step0, perspective));
			vertices.push_back(PerspectiveInterp(v[j], v[(j+1)%3], step1, perspective));

			faces[i] = Eigen::Vector3i(nx, faces[i][j], nx2);
		}
		else {
			int j = 0;
			while (valids[j] == 1)
				j += 1;
			double step0 = 0, step1 = 0;
			if (comp == 0) {
				step0 = (clamp_thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (clamp_thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			else {
				step0 = (clamp_thres - v[j][dim]) / (v[(j + 2) % 3][dim] - v[j][dim]);
				step1 = (clamp_thres - v[j][dim]) / (v[(j + 1) % 3][dim] - v[j][dim]);
			}
			int v0 = faces[i][(j + 2) % 3];
			int nx = vertices.size();
			int nx2 = vertices.size() + 1;
			int v1 = faces[i][(j + 1) % 3];

			vertices.push_back(PerspectiveInterp(v[j], v[(j+2)%3], step0, perspective));
			vertices.push_back(PerspectiveInterp(v[j], v[(j+1)%3], step1, perspective));

			faces[i] = Eigen::Vector3i(v0, nx, nx2);
			faces.push_back(Eigen::Vector3i(v0, nx2, v1));
		}
		i += 1;
	}	
}

void Mesh::Recenter() {
	auto& vertices = vertices_;
	for (int i = 0; i < vertices.size(); ++i) {
		vertices[i][0] -= 0.5;
		vertices[i][1] -= 0.5;
	}
}

void Mesh::ComputeNormals()
{
	auto& face_normals = face_normals_;
	auto& vertices = vertices_;
	auto& faces = faces_;

	face_normals.resize(faces.size());
	for (int i = 0; i < faces.size(); ++i) {
		auto v0 = vertices[faces[i][0]];
		auto v1 = vertices[faces[i][1]];
		auto v2 = vertices[faces[i][2]];
		v0[0] *= v0[2];
		v0[1] *= v0[2];
		v1[0] *= v1[2];
		v1[1] *= v1[2];
		v2[0] *= v2[2];
		v2[1] *= v2[2];
		Eigen::Vector3d d1 = v1 - v0;
		Eigen::Vector3d d2 = v2 - v0;
		Eigen::Vector3d n = d1.cross(d2);
		if (n.norm() > 0)
			n = n.normalized();
		face_normals[i] = n;
	}
}

void Mesh::ComputePlaneParameters() {
	auto& faces = faces_;
	auto& vertices = vertices_;
	auto& params = params_;

	params_.resize(faces.size());
	for (int i = 0; i < faces.size(); ++i) {
		params[i] = PlaneParam(vertices[faces[i][0]],vertices[faces[i][1]],vertices[faces[i][2]]);
	}
}
