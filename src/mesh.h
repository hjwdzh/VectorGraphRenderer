#ifndef VECTORGRAPH_RENDERER_SCENE_H_
#define VECTORGRAPH_RENDERER_SCENE_H_

#include "plane_param.h"
#include "camera.h"
class Mesh
{
public:
	Mesh();

	void LoadFromFile(const char* filename);

	void BoundaryClip(int dim, double clamp_thres, int comp, bool perspective);

	void Recenter();

	void ComputeNormals();

	void ComputePlaneParameters();

	void SaveOBJ(const char* filename, const Camera& camera, bool readjust);

	int FaceNum() const {
		return faces_.size();
	}

	const std::vector<Eigen::Vector3d>& GetVertices() const {
		return vertices_;
	}
	const std::vector<Eigen::Vector3d>& GetFaceNormals() const {
		return face_normals_;
	}
	const std::vector<Eigen::Vector3i>& GetFaces() const {
		return faces_;
	}
	const std::vector<PlaneParam>& GetPlanes() const {
		return params_;
	}

	std::vector<Eigen::Vector3d>& GetVertices() {
		return vertices_;
	}
	std::vector<Eigen::Vector3d>& GetFaceNormals() {
		return face_normals_;
	}
	std::vector<Eigen::Vector3i>& GetFaces() {
		return faces_;
	}
	std::vector<PlaneParam>& GetPlanes() {
		return params_;
	}
private:
	std::vector<Eigen::Vector3d> vertices_;
	std::vector<Eigen::Vector3d> face_normals_;
	std::vector<Eigen::Vector3i> faces_;
	std::vector<PlaneParam> params_;
};

#endif
