#include "camera.h"

#include <fstream>

Camera::Camera()
{}

void Camera::LoadFromFile(const char* filename) {
	std::ifstream is(filename);
	is >> height_ >> width_;
	is >> fx_ >> fy_ >> cx_ >> cy_;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			is >> world2cam_(i, j);
		}
	}
	is >> angle_;
}

void Camera::ApplyExtrinsic(Mesh& mesh) {
	auto& vertices = mesh.GetVertices();

	// apply extrinsic transformation
	for (int i = 0; i < vertices.size(); ++i) {
		Eigen::Vector4d v(vertices[i][0], vertices[i][1], vertices[i][2], 1);
		v = world2cam_ * v;
		vertices[i] = Eigen::Vector3d(v[0],v[1],v[2]);
	}
}

void Camera::ApplyIntrinsic(Mesh& mesh) {
	auto& vertices = mesh.GetVertices();

	// apply intrinsic projection
	for (int i = 0; i < vertices.size(); ++i) {
		Eigen::Vector3d v = vertices[i];
		if (v[2] != 0) {
			vertices[i] = Eigen::Vector3d((v[0]/v[2]*fx_+cx_)/width_,
				(v[1]/v[2]*fy_+cy_)/height_,v[2]);
		} else {
			vertices[i] = Eigen::Vector3d(0, 0, 0);
		}
	}
}